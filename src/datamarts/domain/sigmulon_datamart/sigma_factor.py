from src.datamarts.domain.general.biological_base import BiologicalBase
import multigenomic_api


class SigmaFactor(BiologicalBase):

    def __init__(self, sigma_factor):
        super().__init__(sigma_factor.external_cross_references, sigma_factor.citations, sigma_factor.note)
        promoters = multigenomic_api.promoters.find_by_sigma_factor_id(sigma_factor.id)
        self.sigma_factor = sigma_factor
        self.gene = sigma_factor.genes_id
        self.sigmulon_regulators = promoters
        self.sigmulon_genes = promoters

    def to_dict(self):
        sigma_factor = {
            "_id": self.sigma_factor.id,
            "name": self.sigma_factor.name,
            "abbreviatedName": self.sigma_factor.abbreviated_name,
            "synonyms": self.sigma_factor.synonyms,
            "gene": self.gene,
            "sigmulonRegulators": self.sigmulon_regulators,
            "sigmulonGenes": self.sigmulon_genes
        }
        return sigma_factor

    @property
    def gene(self):
        return self._gene

    @gene.setter
    def gene(self, gene_id):
        gene = multigenomic_api.genes.find_by_id(gene_id)
        self._gene = {
            "_id": gene.id,
            "name": gene.name
        }

    @property
    def sigmulon_regulators(self):
        return self._sigmulon_regulators

    @sigmulon_regulators.setter
    def sigmulon_regulators(self, promoters):
        self._sigmulon_regulators = []
        for promoter in promoters:
            reg_ints = multigenomic_api.regulatory_interactions.find_regulatory_interactions_by_reg_entity_id(promoter.id)
            for reg_int in reg_ints:
                if reg_int.regulator:
                    trans_factor = multigenomic_api.transcription_factors.find_tf_id_by_conformation_id(reg_int.regulator.id)
                    if trans_factor is None:
                        trans_factor = multigenomic_api.transcription_factors.find_by_name(reg_int.regulator.name)
                    if trans_factor:
                        tf_object = {
                            "_id": trans_factor.id,
                            "name": trans_factor.name
                        }
                        if tf_object not in self._sigmulon_regulators:
                            self._sigmulon_regulators.append(tf_object)

    @property
    def sigmulon_genes(self):
        return self._sigmulon_genes

    @sigmulon_genes.setter
    def sigmulon_genes(self, promoters):
        self._sigmulon_genes = []
        for promoter in promoters:
            transcription_units = multigenomic_api.transcription_units.find_by_promoter_id(promoter.id)
            for tu in transcription_units:
                if tu.genes_ids:
                    for gene_id in tu.genes_ids:
                        gene = multigenomic_api.genes.find_by_id(gene_id)
                        products = multigenomic_api.products.find_by_gene_id(gene_id)
                        terms = {
                            'biologicalProcess': [],
                            'cellularComponent': [],
                            'molecularFunction': []
                        }
                        for product in products:
                            if product.terms:
                                for term in product.terms.biological_process:
                                    term = Term(term)
                                    if term.to_dict() not in terms['biologicalProcess']:
                                        terms['biologicalProcess'].append(term.to_dict())

                                for term in product.terms.cellular_component:
                                    term = Term(term)
                                    if term.to_dict() not in terms['cellularComponent']:
                                        terms['cellularComponent'].append(term.to_dict())

                                for term in product.terms.molecular_function:
                                    term = Term(term)
                                    if term.to_dict() not in terms['molecularFunction']:
                                        terms['molecularFunction'].append(term.to_dict())
                        gene_object = {
                            "_id": gene.id,
                            "name": gene.name,
                            "geneOntologyTerms": terms
                        }
                        if gene_object not in self._sigmulon_genes:
                            self.sigmulon_genes.append(gene_object)


class Term(BiologicalBase):

    def __init__(self, term):
        super().__init__([], [], [])
        self.term = term

    def to_dict(self):
        term = {
            '_id': self.term.terms_id,
            'name': self.term.terms_name
        }
        return term
