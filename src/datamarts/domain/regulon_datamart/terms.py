import multigenomic_api
from src.datamarts.domain.general.biological_base import BiologicalBase

class Terms():
    def __init__(self, products_ids):
        super().__init__()
        self.gene_ontology = products_ids
        self.multifun = products_ids

    def to_dict(self):
        terms = {
            "multifun": self.multifun,
            "geneOntology": self.gene_ontology
        }
        return terms

    @property
    def multifun(self):
        return self._multifun

    @multifun.setter
    def multifun(self, products_ids):
        self._multifun = []
        gene_ids = []

        for product_id in products_ids:
            product = multigenomic_api.products.find_by_id(product_id)
            gene = multigenomic_api.genes.find_by_id(product.genes_id)
            terms = gene.terms
            for term in terms:
                term = {
                    'term_id': term.terms_id,
                    'name': term.terms_name,
                    #TODO a que se refieren los gene_ids?
                    'gene_ids': []
                }
                self._multifun.append(term.copy())

    @property
    def gene_ontology(self):
        return self._gene_ontology

    @gene_ontology.setter
    def gene_ontology(self, products_ids):
        self._gene_ontology = {
            'biologicalProcess': [],
            'cellularComponent': [],
            'molecularFunction': []
        }
        for product_id in products_ids:
            product = multigenomic_api.products.find_by_id(product_id)
            terms = product.terms
            for term in terms.biological_process:
                term = Term(term)
                self._gene_ontology['biologicalProcess'].append(term.to_dict())

            for term in terms.cellular_component:
                term = Term(term)
                self._gene_ontology['cellularComponent'].append(term.to_dict())

            for term in terms.molecular_function:
                term = Term(term)
                self._gene_ontology['molecularFunction'].append(term.to_dict())

class Term(BiologicalBase):

    def __init__(self, term):
        super().__init__([], term.citations, [])
        self.term = term

    def to_dict(self):
        term = {
            'id': self.term.terms_id,
            'name': self.term.terms_name,
            #'genes_ids': self.product_ids,
        }
        return term