import multigenomic_api


class Statistics:

    def __init__(self, sigmulon_id):
        promoters = multigenomic_api.promoters.find_by_sigma_factor_id(sigmulon_id)
        transcription_units = get_tus_from_promoters(promoters)
        self.genes = transcription_units
        self.transcription_factors = transcription_units
        self.promoters = promoters
        self.transcription_units = transcription_units
        self.cotranscription_factors = promoters
        self.sigma_factors = transcription_units

    def to_dict(self):
        statistics = {
            "genes": self.genes,
            "transcriptionFactors": self.transcription_factors,
            "promoters": self.promoters,
            "transcriptionUnits": self.transcription_units,
            "cotranscriptionFactors": self.cotranscription_factors,
            "sigmaFactors": self.sigma_factors
        }
        return statistics

    @property
    def genes(self):
        return len(self._genes)

    @genes.setter
    def genes(self, trans_units):
        self._genes = []
        for tu in trans_units:
            if tu.genes_ids:
                for gene_id in tu.genes_ids:
                    if gene_id not in self._genes:
                        self._genes.append(gene_id)

    @property
    def promoters(self):
        return len(self._promoters)

    @promoters.setter
    def promoters(self, promoters):
        self._promoters = promoters

    @property
    def transcription_factors(self):
        return len(self._transcription_factors)

    @transcription_factors.setter
    def transcription_factors(self, trans_units):
        self._transcription_factors = []
        for tu in trans_units:
            if tu.genes_ids:
                for gene_id in tu.genes_ids:
                    products = multigenomic_api.products.find_by_gene_id(gene_id)
                    for product in products:
                        trans_factors = multigenomic_api.transcription_factors.find_tf_id_by_product_id(product.id)
                        for tf in trans_factors:
                            if tf.id not in self._transcription_factors:
                                self._transcription_factors.append(tf.id)

    @property
    def cotranscription_factors(self):
        return len(self._cotranscription_factors)

    @cotranscription_factors.setter
    def cotranscription_factors(self, promoters):
        self._cotranscription_factors = []
        for promoter in promoters:
            reg_interactions = multigenomic_api.regulatory_interactions.find_regulatory_interactions_by_reg_entity_id(promoter.id)
            for ri in reg_interactions:
                if ri.regulator:
                    tf = multigenomic_api.transcription_factors.find_tf_id_by_conformation_id(ri.regulator.id)
                    if tf is None:
                        tf = multigenomic_api.transcription_factors.find_by_name(ri.regulator.name)
                    if tf:
                        if tf.id not in self._cotranscription_factors:
                            self._cotranscription_factors.append(tf.id)

    @property
    def transcription_units(self):
        return len(self._transcription_units)

    @transcription_units.setter
    def transcription_units(self, trans_units):
        self._transcription_units = trans_units

    @property
    def sigma_factors(self):
        return len(self._sigma_factors)

    @sigma_factors.setter
    def sigma_factors(self, trans_units):
        self._sigma_factors = []
        for tu in trans_units:
            if not tu.genes_ids:
                continue
            for gene_id in tu.genes_ids:
                sigma = multigenomic_api.sigma_factors.find_by_gene_id(gene_id)
                if sigma:
                    if sigma.id not in self._sigma_factors:
                        self._sigma_factors.append(sigma.id)


def get_tus_from_promoters(promoters):
    tus = []
    for promoter in promoters:
        trans_units = multigenomic_api.transcription_units.find_by_promoter_id(promoter.id)
        for tu in trans_units:
            if tu not in tus:
                tus.append(tu)
    return tus
