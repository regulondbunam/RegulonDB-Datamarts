from src.datamarts.domain.general.biological_base import BiologicalBase


class BindsSigmaFactor(BiologicalBase):
    def __init__(self, binds_sigma_factor):
        super().__init__(binds_sigma_factor.external_cross_references, binds_sigma_factor.citations, binds_sigma_factor.note)
        self.binds_sigma_factor = binds_sigma_factor

    def to_dict(self):
        binds_sigma_factor_dict = {
            "_id": self.binds_sigma_factor.id,
            "citations": self.citations,
            "name": self.binds_sigma_factor.name,
            "abbreviatedName": self.binds_sigma_factor.abbreviated_name
        }
        return binds_sigma_factor_dict
