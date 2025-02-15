import multigenomic_api
from datetime import datetime
import re


class TFGene:

    @property
    def objects(self):
        ri_objects = multigenomic_api.regulatory_interactions.get_all()
        for ri_object in ri_objects:
            # print(ri_object.id)
            ri_row = TFGene.TFGeneDatamart(ri_object)
            yield ri_row
        del ri_objects

    class TFGeneDatamart:
        def __init__(self, ri):
            self.ri = ri
            self.trans_factor = ri.regulator
            self.genes = self.trans_factor
            self.regulated_genes = ri.regulated_entity
            self.function = ri.function
            self.ri_ev_category = ri

        @property
        def trans_factor(self):
            return self._trans_factor

        @trans_factor.setter
        def trans_factor(self, regulator):
            synonyms = []
            products_ids = []
            tf = multigenomic_api.transcription_factors.find_tf_id_by_conformation_id(regulator.id)
            if tf is None:
                tf = multigenomic_api.transcription_factors.find_by_name(regulator.name)
            if tf is None:
                if regulator.type == "regulatoryComplex":
                    tf = multigenomic_api.transcription_factors.find_by_name(regulator.name)
                    if tf:
                        synonyms.extend(tf.synonyms)
                        synonyms.append(tf.name)
                        self._trans_factor = {
                            "id": tf.id,
                            "name": tf.abbreviated_name,
                            "products_ids": tf.products_ids,
                            "synonyms": synonyms,
                            "type": "TF"
                        }
                    else:
                        complx = multigenomic_api.regulatory_complexes.find_by_id(regulator.id)
                        for product in complx.products:
                            if product.products_id not in products_ids:
                                products_ids.append(product.products_id)
                        synonyms.extend(complx.synonyms)
                        synonyms.append(regulator.name)
                        self._trans_factor = {
                            "id": regulator.id,
                            "name": complx.abbreviated_name or regulator.name,
                            "products_ids": products_ids,
                            "synonyms": synonyms,
                            "type": "Regulator"
                        }
                elif regulator.type == "product":
                    products_ids.append(regulator.id)
                    product = multigenomic_api.products.find_by_id(regulator.id)
                    synonyms.extend(product.synonyms)
                    synonyms.append(product.name)
                    self._trans_factor = {
                        "id": regulator.id,
                        "name": product.abbreviated_name or regulator.name,
                        "products_ids": products_ids,
                        "synonyms": synonyms,
                        "type": "sRNA"
                    }
                else:
                    self._trans_factor = {
                        "id": regulator.id,
                        "name": regulator.name,
                        "products_ids": products_ids,
                        "synonyms": [],
                        "type": "compound"
                    }
            else:
                synonyms.extend(tf.synonyms)
                synonyms.append(tf.name)
                self._trans_factor = {
                    "id": tf.id,
                    "name": tf.abbreviated_name,
                    "products_ids": tf.products_ids,
                    "synonyms": synonyms,
                    "type": "TF"
                }

        @property
        def genes(self):
            return self._genes

        @genes.setter
        def genes(self, tf):
            self._genes = ""
            for product_id in tf["products_ids"]:
                product = multigenomic_api.products.find_by_id(product_id)
                gene = multigenomic_api.genes.find_by_id(product.genes_id)
                self._genes += f"{gene.name};"
            if len(self._genes) > 0:
                self._genes = self._genes[:-1]

        @property
        def function(self):
            return self._function

        @function.setter
        def function(self, function):
            self._function = ""
            if function == "repressor":
                self._function = '-'
            elif function == "activator":
                self._function = '+'

        @property
        def regulated_genes(self):
            return self._regulated_genes

        @regulated_genes.setter
        def regulated_genes(self, regulated_entity):
            self._regulated_genes = []
            if regulated_entity.type == "gene":
                gene = multigenomic_api.genes.find_by_id(regulated_entity.id)
                self._regulated_genes.append(gene)
            elif regulated_entity.type == "promoter":
                trans_units = multigenomic_api.transcription_units.find_by_promoter_id(regulated_entity.id)
                genes_ids = []
                for tu in trans_units:
                    genes_ids.extend(tu.genes_ids)
                for gene_id in genes_ids:
                    gene = multigenomic_api.genes.find_by_id(gene_id)
                    if gene not in self._regulated_genes:
                        self._regulated_genes.append(gene)
            elif regulated_entity.type == "transcriptionUnit":
                tu = multigenomic_api.transcription_units.find_by_id(regulated_entity.id)
                for gene_id in tu.genes_ids:
                    gene = multigenomic_api.genes.find_by_id(gene_id)
                    if gene not in self._regulated_genes:
                        self._regulated_genes.append(gene)

        @property
        def ri_ev_category(self):
            return self._ri_ev_category

        @ri_ev_category.setter
        def ri_ev_category(self, ri):
            self._ri_ev_category = ""
            ev_categories = []
            citations = ri.citations
            if ri.regulatory_sites_id:
                reg_site = multigenomic_api.regulatory_sites.find_by_id(ri.regulatory_sites_id)
                citations.extend(reg_site.citations)
            for citation in citations:
                if citation.evidences_id:
                    citation_dict = multigenomic_api.evidences.find_by_id(citation.evidences_id)
                    if citation_dict.category:
                        if citation_dict.category not in ev_categories:
                            ev_categories.append(citation_dict.category)
            self._ri_ev_category = "|".join(ev_categories)

        def to_row(self):
            response = []
            for gene in self._regulated_genes:
                row = f"{self.trans_factor['id']}" \
                      f"\t{self.trans_factor['name']}" \
                      f"\t{self.trans_factor['type']}" \
                      f"\t{gene['id']}" \
                      f"\t{gene['name']}" \
                      f"\t{','.join(gene['synonyms'])}" \
                      f"\t{','.join(self.trans_factor['synonyms'])}" \
                      f"\t{self.function}" \
                      f"\t{self.ri.confidence_level or '?'}" \
                      f"\t{self.ri_ev_category}"
                response.append(row)
            return response


def remove_similar_items(lista):
    resultado = []
    omitir = set()
    for cadena in lista:
        cad_inicial = cadena.rsplit('\t')
        inicio = '\t'.join(cad_inicial[:-1])[:-1]
        if inicio not in omitir:
            resultado.append(cadena)
            omitir.add(inicio)
    return resultado


def find_existent_items_without_function(list):
    resultado = []
    for i in range(0, len(list) - 1):
        current_item = list[i]
        next_item = list[i + 1]

        cad_inicial_curr = current_item.rsplit('\t')
        cad_inicial_next = next_item.rsplit('\t')
        if cad_inicial_curr[:-3] == cad_inicial_next[:-3]:
            if cad_inicial_curr[-3] != "":
                resultado.append(current_item)
        else:
            resultado.append(current_item)
    return resultado


def find_dual_items(list):
    new_list = []
    for i in range(0, len(list) - 1):
        current_item = list[i]
        next_item = list[i + 1]

        pattern_rep = r'\t-\t'
        pattern_act = r'\t+\t'

        current_item_cat = '\t'.join(current_item.rsplit('\t')[-1:])
        current_item_cl = '\t'.join(current_item.rsplit('\t')[-2:-1])

        current_item_rep = re.sub(pattern_act, "\t-\t", current_item)
        next_item_rep = re.sub(pattern_act, "\t-\t", next_item)
        current_item_act = re.sub(pattern_rep, "\t+\t", current_item)
        next_item_act = re.sub(pattern_rep, "\t+\t", next_item)

        current_item_rep = '\t'.join(current_item_rep.rsplit('\t')[:-2])
        next_item_rep = '\t'.join(next_item_rep.rsplit('\t')[:-2])
        current_item_act = '\t'.join(current_item_act.rsplit('\t')[:-2])
        next_item_act = '\t'.join(next_item_act.rsplit('\t')[:-2])

        if current_item_rep == next_item_rep or current_item_act == next_item_act:
            new_list.append(current_item_act.replace("+", "-+") + f"\t{current_item_cl}\t{current_item_cat}")
            list[i+1] = next_item_act.replace("+", "-+") + f"\t{current_item_cl}\t{current_item_cat}"
        else:
            new_list.append(current_item)
    return remove_similar_items(new_list)


def get_all_rows(rdb_version, citation):
    trans_factors = TFGene()
    tfs_content = [
        "1)regulatorId\t2)regulatorName\t3)regulatorType\t4)geneId\t5)GeneName\t6)genesSynonyms\t7)regulatorSynonyms\t8)function\t9)confidenceLevel\t10)evCategory"]
    for tf in trans_factors.objects:
        tfs_content.extend(tf.to_row())
    tfs_content = list(set(tfs_content))
    tfs_content = sorted(remove_similar_items(tfs_content))
    tfs_content = find_dual_items(find_existent_items_without_function(tfs_content))
    creation_date = datetime.now()
    tfs_doc = {
        "_id": "RDBECOLIDLF00006",
        "fileName": "NetworkRegulatorGene_internal",
        "title": "Complete Regulator-Gene Network Set release 4",
        "fileFormat": "rif-version 1",
        "license": "# RegulonDB is free for academic/noncommercial use\n# User is not entitled to change or erase data sets of the RegulonDB\n# database or to eliminate copyright notices from RegulonDB. Furthermore,\n# User is not entitled to expand RegulonDB or to integrate RegulonDB partly\n# or as a whole into other databank systems, without prior written consent\n# from CCG-UNAM.\n# Please check the license at https://regulondb.ccg.unam.mx/manual/aboutUs/terms-conditions",
        "citation": citation,
        "contact": {
            "person": "RegulonDB Team",
            "webPage": None,
            "email": "regulondb@ccg.unam.mx"
        },
        "version": "1.0",
        "creationDate": f"{creation_date.strftime('%m-%d-%Y')}",
        "columnsDetails": "# Columns:\n# (1) regulatorId. Regulator identifier\n# (2) regulatorName. Regulator Name\n# (3) regulatorType. Gene(s) coding for the TF\n# (4) geneId. Id of the regulated gene\n# (5) geneName. Name of the regulated gene\n# (6) geneSynonyms. Regulated Gene synonyms\n# (7) regulatorSynonyms. Regulator synonyms\n# (8) function. Regulatory Function of the Regulator on the regulated Gene (+ activator, - repressor, -+ dual, ? unknown)\n# (9) confidenceLevel. RI confidence level based on its evidence (Values: Confirmed[C], Strong[S], Weak[W], Unknown[?])\n# (10)evCategory. Category of the RI evidences",
        "content": " \n".join(tfs_content),
        "rdbVersion": rdb_version,
        "description": "Regulatory Network Interactions between Regulators and their regulated genes (Internal File)",
        "group": "REGULATORY NETWORK"
    }
    return tfs_doc
