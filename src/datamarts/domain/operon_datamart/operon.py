from src.datamarts.domain.general.biological_base import BiologicalBase


class Operon(BiologicalBase):

    def __init__(self, operon):
        super().__init__(operon.external_cross_references, operon.citations, operon.note)
        super().__init__([], operon.citations, operon.note)
        self.operon = operon
        self.regulation_positions = operon.regulation_positions

    @property
    def regulation_positions(self):
        return self._regulation_positions

    @regulation_positions.setter
    def regulation_positions(self, regulation_positions):
        self._regulation_positions = {
            "leftEndPosition": regulation_positions.left_end_position or None,
            "rightEndPosition": regulation_positions.right_end_position or None
        }

    def to_dict(self):
        operon = {
            "id": self.operon.id,
            # TODO: Operon doesn't have citations
            "citations": self.operon.citations,
            "name": self.operon.name,
            "regulationPositions": self.regulation_positions,
            "strand": self.operon.strand
            # TODO: add statistics
        }
        return operon