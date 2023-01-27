import multigenomic_api
from src.datamarts.domain.general.remove_items import remove_empty_items


class ListRegulonDM:

    @property
    def objects(self):
        regulator_objects = multigenomic_api.transcription_factors.get_all()
        for regulator_object in regulator_objects:
            print(regulator_object.id)
            regulon_datamart = ListRegulonDM.ListRegulon(regulator_object)
            yield regulon_datamart
        del regulator_objects

    class ListRegulon:

        def __init__(self, trans_factor):
            self.id = trans_factor.id
            self.trans_factor = trans_factor
            self.genes = trans_factor.products_ids
            self.products = trans_factor.products_ids

        @property
        def genes(self):
            return self._genes

        @genes.setter
        def genes(self, products_ids):
            self._genes = []
            for product_id in products_ids:
                prod = multigenomic_api.products.find_by_id(product_id)
                gene = multigenomic_api.genes.find_by_id(prod.genes_id)
                self.genes.append(gene.name)

        @property
        def products(self):
            return self._products

        @products.setter
        def products(self, products_ids):
            self._products = []
            for product_id in products_ids:
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
