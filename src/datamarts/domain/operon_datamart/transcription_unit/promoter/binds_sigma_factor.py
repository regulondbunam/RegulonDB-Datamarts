from src.datamarts.domain.general.biological_base import BiologicalBase
import multigenomic_api


class BindsSigmaFactor(BiologicalBase):
    def __init__(self, binds_sigma_factor):
        super().__init__([], binds_sigma_factor.citations, "")
        self.sigma_factor = multigenomic_api.sigma_factors.find_by_id(binds_sigma_factor.sigma_factors_id)

    def to_dict(self):
        binds_sigma_factor_dict = {
            "_id": self.sigma_factor.id,
            "citations": self.citations,
            "name": self.sigma_factor.name,
            "abbreviatedName": self.sigma_factor.abbreviated_name
        }
        return binds_sigma_factor_dict
