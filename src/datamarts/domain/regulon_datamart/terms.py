import multigenomic_api
from src.datamarts.domain.general.biological_base import BiologicalBase


class Terms():
    def __init__(self, transcription_factor):
        super().__init__()
        self.gene_ontology = transcription_factor.products_ids
        self.multifun = transcription_factor

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
    def multifun(self, transcription_factor):
        self._multifun = []
        regulated_genes = get_genes(transcription_factor.active_conformations)
        for gene in regulated_genes:
            gene_dict = multigenomic_api.genes.find_by_id(gene)
            terms = gene_dict.terms
            for term in terms:
                members_genes = get_members_from_term(term.terms_id)
                term = {
                    'id': term.terms_id,
                    'name': term.terms_name,
                    'gene_ids': intersection_genes(regulated_genes, members_genes)
                }
                if term not in self._multifun:
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
            'term_id': self.term.terms_id,
            'name': self.term.terms_name,
            'gene_ids': []
        }
        return term


def get_members_from_term(term_id):
    term_doc = multigenomic_api.terms.find_by_id(term_id)
    return term_doc.members.genes


def get_genes(active_conformations):
    all_transcription_units = []
    transcription_units = []
    genes = []
    for active_conformation in active_conformations:
        regulatory_interactions = multigenomic_api.regulatory_interactions.find_by_regulator_id(active_conformation.id)
        for ri in regulatory_interactions:
            if ri.regulated_entity.type == "promoter":
                transcription_units = multigenomic_api.transcription_units.find_by_promoter_id(
                    ri.regulated_entity.id)
            elif ri.regulated_entity.type == "transcriptionUnit":
                transcription_units = [].append(multigenomic_api.transcription_units.find_by_id(ri.regulated_entity.id))
            print(transcription_units)
            if transcription_units:
                for tu in transcription_units:
                    if tu not in all_transcription_units:
                        all_transcription_units.append(tu)
    for tu in all_transcription_units:
            for gene_id in tu.genes_ids:
                gene = multigenomic_api.genes.find_by_id(gene_id)
                if gene.id not in genes:
                    genes.append(gene.id)
    return genes


def intersection_genes(regulated_genes, gene_members):
    intersected_genes = [value for value in regulated_genes if value in gene_members]
    return intersected_genes
