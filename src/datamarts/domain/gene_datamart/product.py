import multigenomic_api

from src.datamarts.domain.general.biological_base import BiologicalBase
from src.datamarts.domain.gene_datamart.motif import Motif


class Product(BiologicalBase):
    def __init__(self, product):
        super().__init__(product.external_cross_references, product.citations, product.note)
        self.product = product
        self.motifs = product.id
        self.terms = product.terms
        #self.terms = product.genes_id

    @property
    def motifs(self):
        return self._motifs

    @motifs.setter
    def motifs(self, product_id):
        self._motifs = []

        motifs = multigenomic_api.motifs.find_by_product_id(product_id)
        for motif in motifs:
            motif = Motif(motif)
            motif = motif.to_json()
            self._motifs.append(motif.copy())

    @property
    def terms(self):
        return self._terms

    @terms.setter
    def terms(self, terms):
        self._terms = {
            'biologicalProcess': [],
            'cellularComponent': [],
            'molecularFunction': []
        }
        if terms:
            for term in terms.biological_process:
                term = Term(term)
                self._terms['biologicalProcess'].append(term.to_dict())

            for term in terms.cellular_component:
                term = Term(term)
                self._terms['cellularComponent'].append(term.to_dict())

            for term in terms.molecular_function:
                term = Term(term)
                self._terms['molecularFunction'].append(term.to_dict())

    def get_regulon_id(self):
        tf = multigenomic_api.transcription_factors.find_tf_id_by_product_id(self.product.id)
        if tf:
            return tf[0].id
        return None

    def to_dict(self):
        product = {
            "anticodon": self.product.anticodon,
            "cellularLocations": self.product.locations,
            "citations": self.citations,
            "externalCrossReferences": self.external_cross_references,
            "geneOntologyTerms": self.terms,
            "id": self.product.id,
            "isoelectricPoint": self.product.isoelectric_point,
            # is_regulator: self.product.is_regulator,
            "molecularWeight": self.product.molecular_weight,
            "motifs": self.motifs,
            "name": self.product.name,
            "note": self.formatted_note,
            "regulon_id": self.get_regulon_id(),
            "sequence": self.product.sequence,
            "synonyms": self.product.synonyms,
            "type": self.product.type,
        }
        return product


class Term(BiologicalBase):

    def __init__(self, term):
        super().__init__([], term.citations, [])
        self.term = term

    def to_dict(self):
        term = {
            'id': self.term.terms_id,
            'name': self.term.terms_name,
            #'productIds': self.product_ids,
            'citations': self.citations
        }
        return term
