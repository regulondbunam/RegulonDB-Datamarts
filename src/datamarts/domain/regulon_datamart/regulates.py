import multigenomic_api
from src.datamarts.domain.general.biological_base import BiologicalBase


class Regulates(BiologicalBase):
    def __init__(self, regulator):
        super().__init__(regulator.external_cross_references, regulator.citations, regulator.note)
        self.genes = regulator.products_ids
        self.regulators = regulator
        self.transcription_units = regulator.products_ids
        self.operons = regulator.products_ids
        self.sigmaFactors = regulator

    def to_dict(self):
        regulates = {
            "genes": self.genes,
            # "transcriptionFactors": self.regulators,
            # "transcriptionUnits": self.transcription_units,
            # "operons": self.operons,
            # "sigmaFactors": self.sigmaFactors
        }
        return regulates

    @property
    def genes(self):
        return self._genes

    @genes.setter
    def genes(self, products_ids):
        global transcription_units
        self._genes = []
        for product_id in products_ids:
            regulatory_interactions = multigenomic_api.regulatory_interactions.find_by_regulator_id(product_id)
            print(product_id)
            for ri in regulatory_interactions:
                print(ri.id)
                if ri.regulated_entity.type == "promoter":
                    transcription_units = multigenomic_api.transcription_units.find_by_promoter_id(ri.regulated_entity.id)
                elif ri.regulated_entity.type == "transcriptionUnit":
                    transcription_units = [].append(multigenomic_api.transcription_units.find_by_id(ri.regulated_entity.id))
                product = multigenomic_api.products.find_by_id(product_id)
                for tu in transcription_units:
                    print(tu.id)
                    for gene_id in tu.genes_ids:
                        gene = multigenomic_api.genes.find_by_id(gene_id)
                        terms = gene.terms
                        product_terms = product.terms
                        terms_list = []
                        for term in terms:
                            term = {
                                'id': term.terms_id,
                                'name': term.terms_name
                            }
                            terms_list.append(term.copy())
                        gene_object = {
                            "id": gene.id,
                            "name": gene.name,
                            "terms": {
                                "multifun": terms_list,
                                "geneOntology": gene_ontology_extrac(product_terms)
                            }
                        }
                        if gene_object not in self._genes:
                            self._genes.append(gene_object)

    @property
    def regulators(self):
        return self.regulators

    @regulators.setter
    def regulators(self, regulator):
        pass

    @property
    def transcription_units(self):
        return self._transcription_units

    @transcription_units.setter
    def transcription_units(self, products_ids):
        self._transcription_units = []
        if products_ids:
            for product_id in products_ids:
                product = multigenomic_api.products.find_by_id(product_id)
                gene = multigenomic_api.genes.find_by_id(product.genes_id)
                transcription_units = multigenomic_api.transcription_units.find_by_gene_id(gene.id)
                for tu in transcription_units:
                    self._transcription_units.append({
                        "id": tu.id,
                        "name": tu.name,
                    })


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