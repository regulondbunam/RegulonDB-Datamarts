import multigenomic_api
from src.datamarts.domain.general.remove_items import remove_empty_items


class ListGeneDM:

    @property
    def objects(self):
        gene_objects = multigenomic_api.genes.get_all()
        for gene_obj in gene_objects:
            print(gene_obj.id)
            regulon_datamart = ListGeneDM.ListGene(gene_obj)
            yield regulon_datamart
        del gene_objects

    class ListGene:

        def __init__(self, gene):
            self.id = gene.id
            self.gene = gene
            self.products = gene.id

        @property
        def products(self):
            return self._products

        @products.setter
        def products(self, gene_id):
            self._products = []
            products = multigenomic_api.products.find_by_gene_id(gene_id)
            for prod in products:
                if prod.name not in self._products:
                    self._products.append(prod.name)

        def to_dict(self):
            regulon_datamart = {
                "_id": self.id,
                "name": self.gene.name,
                "synonyms": self.gene.synonyms,
                "productsName": self.products,
                "datamartType": "gene"
            }
            return regulon_datamart


def list_all_gene_datamarts():
    list_genes = ListGeneDM()
    json_genes = []
    for gene in list_genes.objects:
        gene_dict = remove_empty_items(gene.to_dict().copy())
        json_genes.append(remove_empty_items(gene_dict))
    return json_genes
