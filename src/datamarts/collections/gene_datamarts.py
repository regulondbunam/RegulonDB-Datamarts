import multigenomic_api
from src.datamarts.domain.gene_datamart.gene import Gene
from src.datamarts.domain.gene_datamart.product import Product
from src.datamarts.domain.gene_datamart.regulation import Regulation
from src.datamarts.domain.general.biological_base import BiologicalBase
from src.datamarts.domain.general.remove_items import remove_empty_items
from src.datamarts.domain.general.print_progress import print_progress


class GeneDatamarts:

    @property
    def objects(self):
        gene_objects = multigenomic_api.genes.get_all()
        for n, gene_object in enumerate(gene_objects, start=0):
            print_progress(n + 1, len(gene_objects), "Genes")
            yield GeneDatamarts.GeneDatamart(gene_object)
        del gene_objects

    class GeneDatamart:

        def __init__(self, gene):
            self.id = gene.id
            self.gene = gene
            self.products = gene.id
            self.regulation = gene.id
            self.organism = gene.organisms_id

        @property
        def gene(self):
            return self._gene

        @gene.setter
        def gene(self, gene):
            gene = Gene(gene)
            self._gene = gene.to_dict()

        @property
        def products(self):
            return self._products

        @products.setter
        def products(self, gene_id):
            gene_products = multigenomic_api.products.find_by_gene_id(gene_id)
            self._products = []
            for related_product in gene_products:
                related_product = Product(related_product)
                self._products.append(related_product.to_dict().copy())

        @property
        def regulation(self):
            return self._regulation

        @regulation.setter
        def regulation(self, gene_id):
            try:
                operon = multigenomic_api.operons.find_by_gene_id(gene_id)
            except ValueError:
                operon = None
            if operon:
                operon = multigenomic_api.operons.find_by_id(operon.id)
                regulation = Regulation(operon, gene_id)
                self._regulation = regulation.to_dict()
            else:
                self._regulation = None

        @property
        def organism(self):
            return self._organism

        @organism.setter
        def organism(self, organism_id):
            organism = multigenomic_api.organisms.find_by_id(organism_id)
            if organism:
                self._organism = {
                    "_id": organism.id,
                    "name": organism.name
                }
            else:
                self._organism = {
                    "_id": organism_id
                }

        def to_dict(self):
            gene_datamart = {
                "_id": self.id,
                "gene": self.gene,
                "products": self.products,
                # TODO: functions extraction for shine dalgarnos are missing, not exists in MG
                "shineDalgarnos": [],
                "regulation": self.regulation,
                # TODO: functions extraction for growth conditions are missing, not exists in MG
                "growthConditions": [],
                # TODO: organism docs aren't in collection, instead ID is used
                "organism": self.organism,
                "allCitations": BiologicalBase.get_all_citations()
            }
            return gene_datamart


def all_genes_datamarts():
    genes = GeneDatamarts()
    json_genes = []
    for gene in genes.objects:
        gene_dict = remove_empty_items(gene.to_dict().copy())
        json_genes.append(remove_empty_items(gene_dict))
    return json_genes
