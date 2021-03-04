from src.datamarts.domain.general.biological_base import BiologicalBase

class RegulatoryInteractions(BiologicalBase):
    def __init__(self, reg_int, regulatory_sites):
        super().__init__(reg_int.external_cross_references, reg_int.citations, reg_int.note)
        self.regulatory_interactions = reg_int
        self.regulatory_sites = regulatory_sites

    def to_dict(self):
        reg_int_dict = {
            "_id": self.regulatory_interactions.id,
            "centerPosition": self.regulatory_interactions.center_position,
            # TODO: Check this later
            "citations": self.citations,
            "function": self.regulatory_interactions.function,
            "note": self.regulatory_interactions.note,
            "mechanism": self.regulatory_interactions.mechanism,
            "regulatorySite": self.regulatory_sites
        }
        return reg_int_dict