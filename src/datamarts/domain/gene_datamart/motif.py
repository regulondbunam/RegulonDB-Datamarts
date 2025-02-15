from src.datamarts.domain.general.biological_base import BiologicalBase


class Motif(BiologicalBase):

    def __init__(self, motif):
        super().__init__(motif.external_cross_references, motif.citations, motif.note)
        self.motif = motif

    def to_json(self):
        motif = {
            "description": self.motif.description,
            "_id": self.motif.id,
            "leftEndPosition": self.motif.left_end_position,
            "note": self.formatted_note,
            "rightEndPosition": self.motif.right_end_position,
            "sequence": self.motif.sequence,
            "type": self.motif.class_,
            "dataSource": self.motif.data_source,
            "citations": self.citations
        }
        return motif
