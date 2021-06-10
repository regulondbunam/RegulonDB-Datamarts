import multigenomic_api
from src.datamarts.domain.regulon_datamart.summary import insert_ri_bs_counts, get_counts_of_regulated_object

class Summary:
    def __init__(self, srna, reg_ints):
        transcription_units = get_all_transcription_units(srna.id)
        self.genes = transcription_units
        self.transcription_factors = transcription_units
        self.transcription_units = transcription_units
        self.sigma_factors = transcription_units
        self.regulatory_interactions = reg_ints

    def to_dict(self):
        summary = {
            "genes": self.genes,
            "transcriptionFactors": self.transcription_factors,
            "transcriptionUnits": self.transcription_units,
            "sigmaFactors": self.sigma_factors,
            "regulatoryInteractions": self.regulatory_interactions
        }
        for obj, regulated_list in summary.items():
            count = get_counts_of_regulated_object(regulated_list)
            summary.update({obj: count})
        return insert_ri_bs_counts(self.regulatory_interactions, summary)

    @property
    def genes(self):
        return self._genes

    @genes.setter
    def genes(self, trans_units):
        self._genes = []
        for tu_object in trans_units:
            tu = tu_object["transcription_unit"]
            for gene_id in tu.genes_ids:
                gene = multigenomic_api.genes.find_by_id(gene_id)
                gene_object = {
                    "id": gene.id,
                    "name": gene.name,
                    "function": tu_object["function"]
                }
                self._genes = insert_element_on(gene_object, self._genes)

    @property
    def transcription_factors(self):
        return self._transcription_factors

    @transcription_factors.setter
    def transcription_factors(self, trans_units):
        self._transcription_factors = []
        for tu_object in trans_units:
            tu = tu_object["transcription_unit"]
            for gene_id in tu.genes_ids:
                products = multigenomic_api.products.find_by_gene_id(gene_id)
                for product in products:
                    trans_factors = multigenomic_api.transcription_factors.find_tf_id_by_active_conformation_id(
                        product.id)
                    for tf in trans_factors:
                        tf_object = {
                            "id": tf.id,
                            "name": tf.name,
                            "function": tu_object["function"]
                        }
                        self._transcription_factors = insert_element_on(tf_object, self._transcription_factors)

    @property
    def transcription_units(self):
        return self._transcription_units

    @transcription_units.setter
    def transcription_units(self, trans_units):
        self._transcription_units = []
        for tu_object in trans_units:
            tu = tu_object["transcription_unit"]
            tu_object = {
                "id": tu.id,
                "name": tu.name,
                "function": tu_object["function"]
            }
            self._transcription_units = insert_element_on(tu_object, self._transcription_units)

    @property
    def sigma_factors(self):
        return self._sigma_factors

    @sigma_factors.setter
    def sigma_factors(self, trans_units):
        self._sigma_factors = []
        for tu_object in trans_units:
            tu = tu_object["transcription_unit"]
            if tu.promoters_id:
                promoter = multigenomic_api.promoters.find_by_id(tu.promoters_id)
                if promoter.binds_sigma_factor:
                    sigma_factor = multigenomic_api.sigma_factors.find_by_id(
                        promoter.binds_sigma_factor.sigma_factors_id)
                    sigma_factor_object = {
                        "id": sigma_factor.id,
                        "name": sigma_factor.name,
                        "function": tu_object["function"]
                    }
                    self._sigma_factors = insert_element_on(sigma_factor_object, self._sigma_factors)


def get_all_transcription_units(srna_id):
    all_transcription_units = []
    regulatory_interactions = multigenomic_api.regulatory_interactions.find_by_regulator_id(srna_id)
    for ri in regulatory_interactions:
        if ri.regulated_entity.type == "transcriptionUnit":
            transcription_unit = multigenomic_api.transcription_units.find_by_id(ri.regulated_entity.id)
            tu_object = {
                "transcription_unit": transcription_unit,
                "function": ri.function
            }
            if tu_object not in all_transcription_units:
                all_transcription_units.append(tu_object)
    return all_transcription_units


def insert_element_on(reg_object, regulated_list):
    for regulated in regulated_list:
        if reg_object["id"] == regulated["id"]:
            if reg_object['function'] != regulated['function']:
                reg_object["function"] = "dual"
                regulated["function"] = "dual"
    if reg_object not in regulated_list:
        regulated_list.append(reg_object)
    return regulated_list
