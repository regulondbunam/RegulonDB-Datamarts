from src.datamarts.domain.general.biological_base import BiologicalBase
from src.datamarts.domain.general.additiveEvidences import AdditiveEvidences


class RegulatoryInteractions(BiologicalBase):
    def __init__(self, reg_int, regulatory_sites):
        super().__init__(reg_int.external_cross_references, reg_int.citations, reg_int.note)
        self.regulatory_interactions = reg_int
        self.regulatory_sites = regulatory_sites

    def to_dict(self):
        additive_evs = AdditiveEvidences(self.citations + self.regulatory_sites.get("citations", []))
        reg_int_dict = {
            "_id": self.regulatory_interactions.id,
            "relativeCenterPosition": self.regulatory_interactions.dist_site_promoter,
            "citations": self.citations,
            "function": self.regulatory_interactions.function,
            "note": self.formatted_note,
            "mechanism": self.regulatory_interactions.mechanism,
            "regulatorySite": self.regulatory_sites,
            "additiveEvidences": additive_evs.to_dict(),
            "confidenceLevel": self.regulatory_interactions.confidence_level
        }
        return reg_int_dict
