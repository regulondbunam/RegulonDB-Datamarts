import multigenomic_api
from datetime import datetime


class TFTU:

    @property
    def objects(self):
        ri_objects = multigenomic_api.regulatory_interactions.get_all()
        for ri_object in ri_objects:
            # print(ri_object.id)
            ri_row = TFTU.TFTUDatamart(ri_object)
            yield ri_row
        del ri_objects

    class TFTUDatamart:
        def __init__(self, ri):
            self.ri = ri
            self.trans_factor = ri.regulator
            self.genes = self.trans_factor
            self.regulated_tus = ri.regulated_entity
            self.function = ri.function

        @property
        def trans_factor(self):
            return self._trans_factor

        @trans_factor.setter
        def trans_factor(self, regulator):
            tf = multigenomic_api.transcription_factors.find_tf_id_by_conformation_id(regulator.id)
            if tf is None:
                tf = multigenomic_api.transcription_factors.find_by_name(regulator.name)
            if tf:
                self._trans_factor = {
                    "id": tf.id,
                    "name": tf.abbreviated_name,
                    "products_ids": tf.products_ids
                }
            else:
                products_ids = []
                if regulator.type == "regulatoryComplex":
                    complex = multigenomic_api.regulatory_complexes.find_by_id(regulator.id)
                    for product in complex.products:
                        if product.products_id not in products_ids:
                            products_ids.append(product.products_id)
                    self._trans_factor = {
                        "id": regulator.id,
                        "name": complex.abbreviated_name or regulator.name,
                        "products_ids": products_ids
                    }
                elif regulator.type == "product":
                    products_ids.append(regulator.id)
                    product = multigenomic_api.products.find_by_id(regulator.id)
                    self._trans_factor = {
                        "id": regulator.id,
                        "name": product.abbreviated_name or regulator.name,
                        "products_ids": products_ids
                    }
                else:
                    self._trans_factor = {
                        "id": regulator.id,
                        "name": regulator.name,
                        "products_ids": products_ids
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
        def regulated_tus(self):
            return self._regulated_tus

        @regulated_tus.setter
        def regulated_tus(self, regulated_entity):
            self._regulated_tus = []
            if regulated_entity.type == "promoter":
                trans_units = multigenomic_api.transcription_units.find_by_promoter_id(regulated_entity.id)
                for tu in trans_units:
                    if tu not in self._regulated_tus:
                        self._regulated_tus.append(tu)
            elif regulated_entity.type == "transcriptionUnit":
                tu = multigenomic_api.transcription_units.find_by_id(regulated_entity.id)
                self._regulated_tus.append(tu)

        def to_row(self):
            response = []
            for tu in self._regulated_tus:
                promoter_name = ""
                if tu["promoters_id"]:
                    promoter = multigenomic_api.promoters.find_by_id(tu["promoters_id"])
                    promoter_name = promoter.name
                row = f"{self.trans_factor['id']}" \
                      f"\t{self.trans_factor['name']}" \
                      f"\t{self.genes}" \
                      f"\t{tu['id']}" \
                      f"\t{tu['name']}[{promoter_name}]" \
                      f"\t{self.function}" \
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


def get_all_rows(rdb_version, citation):
    trans_factors = TFTU()
    tfs_content = [
        "1)regulatorId	2)regulatorName	3)RegulatorGeneName	4)regulatedId	5)regulatedName	6)function	7)confidenceLevel"]
    for tf in trans_factors.objects:
        tfs_content.extend(tf.to_row())
    tfs_content = list(set(tfs_content))
    tfs_content = sorted(remove_similar_items(tfs_content))
    tfs_content = find_dual_items(find_existent_items_without_function(tfs_content))
    creation_date = datetime.now()
    tfs_doc = {
        "_id": "RDBECOLIDLF00015",
        "fileName": "NetworkRegulatorTU",
        "title": "Complete Regulator-TU Network Set",
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
        "columnsDetails": "# Columns:\n# (1) regulatorId. Regulator identifier\n# (2) regulatorName. Regulator Name\n# (3) regulatorGeneName. Gene(s) coding for the TF\n# (4) regulatedId. Transcription Unit ID regulated by the TF (regulated Transcription Unit)\n# (5) regulatedName. Transcription Unit regulated by the TF (regulated Transcription Unit)\n# (6) function. Regulatory Function of the TF on the regulated Transcription Unit (+ activator, - repressor, +- dual, ? unknown)\n# (7) confidenceLevel. RI confidence level based on its evidence (Values: Confirmed[C], Strong[S], Weak[W], Unknown[?])",
        "content": " \n".join(tfs_content),
        "rdbVersion": rdb_version,
        "description": "Regulatory Network Interactions between Regulators and their regulated Transcription Units",
        "group": "REGULATORY NETWORK"
    }
    return tfs_doc
