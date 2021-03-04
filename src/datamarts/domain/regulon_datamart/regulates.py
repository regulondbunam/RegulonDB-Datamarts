import multigenomic_api
from src.datamarts.domain.general.biological_base import BiologicalBase


class Regulates(BiologicalBase):
    def __init__(self, regulator):
        super().__init__(regulator.external_cross_references, regulator.citations, regulator.note)
        self.genes = regulator.products_ids
        self.regulators = regulator
        self.transcriptionUnits = regulator.products_ids
        self.operons = regulator.products_ids
        self.sigmaFactors = regulator

    def to_dict(self):
        regulates = {
            "genes": self.genes,
            # "regulators": self.regulators,
            # "transcriptionUnits": self.transcriptionUnits,
            # "operons": self.operons,
            # "sigmaFactors": self.sigmaFactors
        }
        return regulates

    @property
    def genes(self):
        return self._genes

    @genes.setter
    def genes(self, products_ids):
        self._genes = []
        for product_id in products_ids:
            product = multigenomic_api.products.find_by_id(product_id)
            gene = multigenomic_api.genes.find_by_id(product.genes_id)
            terms = gene.terms
            product_terms = product.terms
            terms_list = []
            for term in terms:
                term = {
                    'id': term.terms_id,
                    'name': term.terms_name
                }
                terms_list.append(term.copy())
            self._genes.append({
                "id": gene.id,
                "name": gene.name,
                # TODO de donde se obtiene esto?
                # "function": None,
                "terms": {
                    'multifun': terms_list,
                    "geneOntology": gene_ontology_extrac(product_terms)
                }
            })

    @property
    def regulators(self):
        return self.regulators

    @regulators.setter
    def regulators(self, regulator):
        pass


def gene_ontology_extrac(terms):
    terms_dict = {
        'biologicalProcess': [],
        'cellularComponent': [],
        'molecularFunction': []
    }
    if terms:
        for term in terms.biological_process:
            term = Term(term)
            terms_dict['biologicalProcess'].append(term.to_dict())

        for term in terms.cellular_component:
            term = Term(term)
            terms_dict['cellularComponent'].append(term.to_dict())

        for term in terms.molecular_function:
            term = Term(term)
            terms_dict['molecularFunction'].append(term.to_dict())
    return terms_dict

class Term(BiologicalBase):

    def __init__(self, term):
        super().__init__([], term.citations, [])
        self.term = term

    def to_dict(self):
        term = {
            'term_id': self.term.terms_id,
            'name': self.term.terms_name,
            #'gene_ids': self.product_ids
        }
        return term