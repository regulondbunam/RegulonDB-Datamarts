import multigenomic_api
from datetime import datetime


class SigmaGene:

    @property
    def objects(self):
        ri_objects = multigenomic_api.regulatory_interactions.get_all()
        for ri_object in ri_objects:
            # print(ri_object.id)
            ri_row = SigmaGene.SigmaGeneDatamart(ri_object)
            yield ri_row
        del ri_objects

    class SigmaGeneDatamart:
        def __init__(self, ri):
            self.ri = ri
            self.regulated_genes = ri.regulated_entity
            self.sigma = ri.regulated_entity
            self.function = ri.function
            self.ri_evidences = ri.citations

        @property
        def sigma(self):
            return self._sigma

        @sigma.setter
        def sigma(self, regulated_entity):
            self._sigma = ""
            if regulated_entity.type == "promoter":
                promoter = multigenomic_api.promoters.find_by_id(regulated_entity.id)
                if promoter.binds_sigma_factor:
                    sigma = multigenomic_api.sigma_factors.find_by_id(promoter.binds_sigma_factor.sigma_factors_id)
                    self._sigma = sigma

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
            if regulated_entity.type == "promoter":
                trans_units = multigenomic_api.transcription_units.find_by_promoter_id(regulated_entity.id)
                genes_ids = []
                for tu in trans_units:
                    genes_ids.extend(tu.genes_ids)
                for gene_id in genes_ids:
                    gene = multigenomic_api.genes.find_by_id(gene_id)
                    if gene not in self._regulated_genes:
                        self._regulated_genes.append(gene)

        @property
        def ri_evidences(self):
            return self._ri_evidences

        @ri_evidences.setter
        def ri_evidences(self, citations):
            self._ri_evidences = []
            for citation in citations:
                if citation.evidences_id:
                    citation_dict = multigenomic_api.evidences.find_by_id(citation.evidences_id)
                    if citation_dict.code not in self._ri_evidences:
                        self._ri_evidences.append(citation_dict.code)
            self._ri_evidences = ",".join(self._ri_evidences)
            self._ri_evidences = f"[{self._ri_evidences}]"

        def to_row(self):
            response = []
            if self.sigma != "":
                for gene in self._regulated_genes:
                    if self.function != "":
                        row = f"{self.sigma.abbreviated_name}" \
                              f"\t{gene['name']}" \
                              f"\t{self.function}" \
                              f"\t{self.ri_evidences}" \
                              f"\t{self.ri.confidence_level or '?'}"
                        response.append(row)
            return response


def remove_similar_items(lista):
    resultado = []
    omitir = set()
    for cadena in lista:
        inicio = cadena[:-1]
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
        if cad_inicial_curr[:-2] == cad_inicial_next[:-2]:
            if cad_inicial_curr[-2] != "":
                resultado.append(current_item)
        else:
            resultado.append(current_item)
    return resultado


def find_dual_items(list):
    new_list = []
    for i in range(0, len(list) - 1):
        current_item = list[i]
        next_item = list[i + 1]

        current_item_cl = current_item[-1:]

        current_item_rep = current_item.replace("+", "-")[:-1]
        next_item_rep = next_item.replace("+", "-")[:-1]
        current_item_act = current_item.replace("-", "+")[:-1]
        next_item_act = next_item.replace("-", "+")[:-1]

        if current_item_rep == next_item_rep or current_item_act == next_item_act:
            new_list.append(current_item_act.replace("+", "-+") + current_item_cl)
            list[i+1] = next_item_act.replace("+", "-+") + current_item_cl
        else:
            new_list.append(current_item)
    return remove_similar_items(new_list)


def get_all_rows():
    trans_factors = SigmaGene()
    tfs_content = [
        "1)sigmaName\t2)regulatedName\t3)function\t4)evidences\t5)confidenceLevel"]
    for tf in trans_factors.objects:
        tfs_content.extend(tf.to_row())
    tfs_content = list(set(tfs_content))
    tfs_content = sorted(remove_similar_items(tfs_content))
    tfs_content = find_dual_items(find_existent_items_without_function(tfs_content))
    creation_date = datetime.now()
    tfs_doc = {
        "_id": "RDBECOLIDLF00018",
        "fileName": "NetworkSigmaGene",
        "title": "Complete Sigma-Gene Network Set",
        "fileFormat": "rif-version 1",
        "license": "RegulonDB is free for academic/noncommercial use\nUser is not entitled to change or erase data sets of the RegulonDB\ndatabase or to eliminate copyright notices from RegulonDB. Furthermore,\nUser is not entitled to expand RegulonDB or to integrate RegulonDB partly\nor as a whole into other databank systems, without prior written consent\nfrom CCG-UNAM.\nPlease check the license at https://regulondb.ccg.unam.mx/manual/aboutUs/terms-conditions",
        "citation": "Salgado H., Gama-Castro S. et al (2023). RegulonDB 12.0: A Comprehensive resource of transcriptional regulation in E. coli K-12",
        "contact": {
            "person": "RegulonDB Team",
            "webPage": None,
            "email": "regulondb@ccg.unam.mx"
        },
        "version": "1.0",
        "creationDate": f"{creation_date.strftime('%m-%d-%Y')}",
        "columnsDetails": "Columns:\n(1) sigmaName. Sigma Name\n(2) regulatedName. Gene regulated by the Sigma Factor (regulated Gene)\n(3) function. Regulatory Function of the Regulator on the regulated Gene (+ activator, - repressor, -+ dual, ? unknown)\n(4)riEvidences. Evidence that supports the existence of the regulatory interaction\n(5) confidenceLevel. RI confidence level based on its evidence (Values: Confirmed[C], Strong[S], Weak[W], Unknown[?])",
        "content": " \n".join(tfs_content),
        "rdbVersion": "12.0"
    }
    return tfs_doc
