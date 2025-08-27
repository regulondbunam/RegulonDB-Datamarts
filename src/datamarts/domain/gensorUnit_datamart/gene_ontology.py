import multigenomic_api
from src.datamarts.domain.general.biological_base import BiologicalBase


class GeneOntology:
    def __init__(self, regulator):
        super().__init__()
        tf = multigenomic_api.transcription_factors.find_by_id(regulator['_id'])
        regulated_genes = get_genes(tf)
        tf["regulated_genes"] = regulated_genes
        self.gene_ontology = tf

    def to_dict(self):
        return self.gene_ontology

    @property
    def gene_ontology(self):
        return self._gene_ontology

    @gene_ontology.setter
    def gene_ontology(self, regulator):
        self._gene_ontology = {
            'biologicalProcess': [],
            'cellularComponent': [],
            'molecularFunction': []
        }
        regulated_genes = regulator.regulated_genes
        for reg_gene in regulated_genes:
            products = multigenomic_api.products.find_by_gene_id(reg_gene)
            for product in products:
                terms = product.terms
                if terms:
                    if terms:
                        self._add_terms(terms.biological_process, "biologicalProcess", reg_gene)
                        self._add_terms(terms.cellular_component, "cellularComponent", reg_gene)
                        self._add_terms(terms.molecular_function, "molecularFunction", reg_gene)

    def _add_terms(self, term_list, ontology_key, regulated_gene):
        gene = multigenomic_api.genes.find_by_id(regulated_gene)

        for term in term_list:
            term_obj = Term(term, gene.name)
            term_dict = term_obj.to_dict()

            existing = next(
                (t for t in self._gene_ontology[ontology_key] if t['_id'] == term_dict['_id']),
                None
            )

            if existing is None:
                self._gene_ontology[ontology_key].append(term_dict)
            else:
                if gene.name not in existing['genes']:
                    existing['genes'].append(gene.name)

        self._gene_ontology[ontology_key].sort(
            key=lambda t: len(t['genes']),
            reverse=True
        )

class Term(BiologicalBase):

    def __init__(self, term, gene_name):
        super().__init__([], term.citations, [])
        self.gene = gene_name
        self.term = term

    def to_dict(self):
        term = {
            '_id': self.term.terms_id,
            'name': self.term.terms_name,
            'genes': [self.gene]
        }
        return term

def get_genes(tf):
    all_transcription_units = []
    transcription_units = []
    genes = []
    reg_ints = []
    for active_conformation in tf.active_conformations:
        reg_ints.extend(multigenomic_api.regulatory_interactions.find_by_regulator_id(active_conformation.id))
    for product_id in tf.products_ids:
        reg_ints.extend(multigenomic_api.regulatory_interactions.find_by_regulator_id(product_id))
    for ri in reg_ints:
        if ri.regulated_entity.type == "promoter":
            transcription_units = multigenomic_api.transcription_units.find_by_promoter_id(
                ri.regulated_entity.id)
        elif ri.regulated_entity.type == "transcriptionUnit":
            transcription_units = [].append(multigenomic_api.transcription_units.find_by_id(ri.regulated_entity.id))
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