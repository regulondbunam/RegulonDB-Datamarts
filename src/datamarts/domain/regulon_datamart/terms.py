import multigenomic_api
from src.datamarts.domain.general.biological_base import BiologicalBase


class Terms:
    def __init__(self, transcription_factor):
        super().__init__()
        regulated_genes = get_genes(transcription_factor.active_conformations)
        self.gene_ontology = [transcription_factor.products_ids, regulated_genes]
        self.multifun = regulated_genes

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
    def multifun(self, regulated_genes):
        self._multifun = []
        for gene in regulated_genes:
            gene_dict = multigenomic_api.genes.find_by_id(gene)
            terms = gene_dict.terms
            for term in terms:
                members_genes = get_members_from_term(term.terms_id)
                genes = intersection_genes(regulated_genes, members_genes)
                term = {
                    '_id': term.terms_id,
                    'name': term.terms_name,
                    'genes': gene_object(genes)
                }
                if term not in self._multifun:
                    self._multifun.append(term.copy())

    @property
    def gene_ontology(self):
        return self._gene_ontology

    @gene_ontology.setter
    def gene_ontology(self, params):
        self._gene_ontology = {
            'biologicalProcess': [],
            'cellularComponent': [],
            'molecularFunction': []
        }
        products_ids = params[0]
        regulated_genes = params[1]
        for product_id in products_ids:
            product = multigenomic_api.products.find_by_id(product_id)
            terms = product.terms
            for term in terms.biological_process:
                term = Term(term, regulated_genes)
                term_dict = term.to_dict()
                if term_dict not in self._gene_ontology['biologicalProcess']:
                    self._gene_ontology['biologicalProcess'].append(term_dict)

            for term in terms.cellular_component:
                term = Term(term, regulated_genes)
                term_dict = term.to_dict()
                if term_dict not in self._gene_ontology['cellularComponent']:
                    self._gene_ontology['cellularComponent'].append(term_dict)

            for term in terms.molecular_function:
                term = Term(term, regulated_genes)
                term_dict = term.to_dict()
                if term_dict not in self._gene_ontology['molecularFunction']:
                    self._gene_ontology['molecularFunction'].append(term_dict)


class Term(BiologicalBase):

    def __init__(self, term, regulated_genes):
        super().__init__([], term.citations, [])
        self.regulated_genes = regulated_genes
        self.term = term

    def to_dict(self):
        product_term_members = get_members_from_term(self.term.terms_id)
        genes = intersection_genes(self.regulated_genes, product_term_members)
        term = {
            '_id': self.term.terms_id,
            'name': self.term.terms_name,
            'genes': gene_object(genes)
        }
        return term


def get_members_from_term(term_id):
    term_doc = multigenomic_api.terms.find_by_id(term_id)
    if term_doc.members.genes:
        return term_doc.members.genes
    else:
        products = term_doc.members.products
        return get_genes_of_products(products)


def get_genes_of_products(products):
    genes = []
    for product in products:
        product_dict = multigenomic_api.products.find_by_id(product)
        if product_dict.genes_id not in genes:
            genes.append(product_dict.genes_id)
    return genes


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


def gene_object(genes):
    gene_list = []
    for gene_id in genes:
        gene = multigenomic_api.genes.find_by_id(gene_id)
        gene_list.append({
            "_id": gene.id,
            "name": gene.name
        })
    return gene_list
