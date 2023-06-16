import multigenomic_api
from src.datamarts.domain.general.remove_items import remove_empty_items


class ListRegulonDM:

    @property
    def objects(self):
        regulator_objects = get_all_regulators(multigenomic_api.regulatory_interactions.get_all())
        tf_objects = multigenomic_api.transcription_factors.get_all()
        for tf_obj in tf_objects:
            print(tf_obj.id)
            tf_obj["regulator_type"] = "transcriptionFactor"
            regulon_datamart = ListRegulonDM.ListRegulon(tf_obj)
            yield regulon_datamart
        for reg_obj in regulator_objects:
            print(reg_obj.id)
            regulator = reg_obj
            if reg_obj.type == "srna":
                regulator = multigenomic_api.products.find_by_id(reg_obj.id)
            if reg_obj.type == "compound":
                regulator = multigenomic_api.regulatory_continuants.find_by_id(reg_obj.id)
            regulator["regulator_type"] = reg_obj.type
            regulon_datamart = ListRegulonDM.ListRegulon(regulator)
            yield regulon_datamart
        del regulator_objects
        del tf_objects

    class ListRegulon:

        def __init__(self, regulator):
            self.id = regulator.id
            self.trans_factor = regulator
            self.genes = regulator
            self.products = regulator

        @property
        def genes(self):
            return self._genes

        @genes.setter
        def genes(self, regulator):
            self._genes = []
            if regulator.regulator_type == "transcriptionFactor":
                for product_id in regulator.products_ids:
                    prod = multigenomic_api.products.find_by_id(product_id)
                    gene = multigenomic_api.genes.find_by_id(prod.genes_id)
                    self.genes.append(gene.name)
            elif regulator.regulator_type == "srna":
                if regulator.genes_id:
                    gene = multigenomic_api.genes.find_by_id(regulator.genes_id)
                    self.genes.append(gene.name)

        @property
        def products(self):
            return self._products

        @products.setter
        def products(self, regulator):
            self._products = []
            if regulator.regulator_type == "transcriptionFactor":
                for product_id in regulator.products_ids:
                    prod = multigenomic_api.products.find_by_id(product_id)
                    if prod.name not in self._products:
                        self._products.append(prod.name)

        def to_dict(self):
            regulon_datamart = {
                "_id": self.id,
                "name": self.trans_factor.name,
                "synonyms": self.trans_factor.synonyms,
                "encodedGenes": self.genes,
                "productsName": self.products,
                "datamartType": "regulon"
            }
            return regulon_datamart


def list_all_regulon_datamarts():
    list_regulons = ListRegulonDM()
    json_regulons = []
    for regulon in list_regulons.objects:
        regulon_dict = remove_empty_items(regulon.to_dict().copy())
        json_regulons.append(remove_empty_items(regulon_dict))
    return json_regulons


def get_all_regulators(reg_ints):
    regulators = []
    for ri in reg_ints:
        if ri.regulator.type == "regulatoryContinuant":
            regulator = ri.regulator
            regulator["type"] = "compound"
            if ri.regulator not in regulators:
                regulators.append(ri.regulator)
        if ri.regulator.type == "product":
            product = multigenomic_api.products.find_by_id(ri.regulator.id)
            if product.type == "small RNA":
                regulator = ri.regulator
                regulator["type"] = "srna"
                if ri.regulator not in regulators:
                    regulators.append(ri.regulator)
    return regulators