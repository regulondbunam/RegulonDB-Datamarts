from src.datamarts.domain.general.biological_base import BiologicalBase


class Terminators(BiologicalBase):
    def __init__(self, terminator):
        super().__init__(terminator.external_cross_references, terminator.citations, terminator.note)
        self.terminator = terminator
        self.transcription_termination_site = terminator

    @property
    def transcription_termination_site(self):
        return self._transcription_termination_site

    @transcription_termination_site.setter
    def transcription_termination_site(self, terminator):
        if terminator.transcriptionTerminationSite:
            self._transcription_termination_site = {
                "leftEndPosition": terminator.transcriptionTerminationSite.left_end_position,
                "rightEndPosition": terminator.transcriptionTerminationSite.right_end_position,
                "type:": terminator.transcriptionTerminationSite.type
            }

    def to_dict(self):
        terminator_dict = {
            "_id": self.terminator.id,
            "citations": self.citations,
            "sequence": self.terminator.sequence,
            "transcriptionTerminationSite": self.transcription_termination_site,
            "class": self.terminator.class_
        }
        return terminator_dict
