import multigenomic_api

from src.datamarts.domain.general.biological_base import BiologicalBase


class Operon(BiologicalBase):

    def __init__(self, operon):
        super().__init__(operon.external_cross_references, operon.citations, operon.note)
        super().__init__([], operon.citations, operon.note)
        trans_units = multigenomic_api.transcription_units.find_by_operon_id(operon.id)
        self.operon = operon
        self.transcription_units = trans_units
        self.promoters = trans_units
        self.genes = trans_units
        self.regulation_positions = operon.regulation_positions

    @property
    def transcription_units(self):
        return self._transcription_units

    @transcription_units.setter
    def transcription_units(self, trans_units):
        self._transcription_units = trans_units

    @property
    def promoters(self):
        return self._promoters

    @promoters.setter
    def promoters(self, trans_units):
        self._promoters = []
        for tu in trans_units:
            if tu.promoters_id:
                self._promoters.append(tu.promoters_id)
        self._promoters = set(list(self._promoters))

    @property
    def genes(self):
        return self._genes

    @genes.setter
    def genes(self, trans_units):
        self._genes = []
        for tu in trans_units:
            if tu.genes_ids:
                self._genes.extend(tu.genes_ids)

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
            "citations": self.operon.citations,
            "name": self.operon.name,
            "regulationPositions": self.regulation_positions,
            "strand": self.operon.strand,
            "statistics": {
                "transcriptionUnits": len(self.transcription_units),
                "promoters": len(self.promoters),
                "genes": len(self.genes)
            }
        }
        return operon