from src.datamarts.domain.general.biological_base import BiologicalBase


class BindsSigmaFactor(BiologicalBase):
    def __init__(self, binds_sigma_factor):
        super().__init__(binds_sigma_factor.external_cross_references, binds_sigma_factor.citations, binds_sigma_factor.note)
        self.binds_sigma_factor = binds_sigma_factor

    def to_dict(self):
        binds_sigma_factor_dict = {
            "sigmaFactor_id": self.binds_sigma_factor.id,
            "citations": self.citations,
            "sigmaFactor_name": self.binds_sigma_factor.name
        }
        return binds_sigma_factor_dict