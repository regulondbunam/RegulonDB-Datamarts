import multigenomic_api

from src.datamarts.domain.general.biological_base import BiologicalBase
from src.datamarts.domain.operon_datamart.regulator_binding_sites import RegulatoryBindingSites
from src.datamarts.domain.operon_datamart.transcription_unit.promoter.binds_sigma_factor import BindsSigmaFactor
from src.datamarts.domain.general.additiveEvidences import AdditiveEvidences


class Promoters(BiologicalBase):
    def __init__(self, promoter):
        super().__init__(promoter.external_cross_references, promoter.citations, promoter.note)
        self.promoter = promoter
        self.regulator_binding_sites = promoter.id
        self.binds_sigma_factor = promoter
        self.boxes = promoter
        self.transcription_start_site = promoter

    @property
    def binds_sigma_factor(self):
        return self._binds_sigma_factor

    @binds_sigma_factor.setter
    def binds_sigma_factor(self, promoter):
        self._binds_sigma_factor = {}
        if promoter.binds_sigma_factor:
            binds_sigma_factor = BindsSigmaFactor(promoter.binds_sigma_factor)
            self._binds_sigma_factor = binds_sigma_factor.to_dict()

    @property
    def boxes(self):
        return self._boxes

    @boxes.setter
    def boxes(self, promoter):
        self._boxes = []
        boxes = promoter.boxes
        for box in boxes:
            promoter_boxes = ({
                "leftEndPosition": box.left_end_position,
                "rightEndPosition": box.right_end_position,
                "sequence": box.sequence,
                "type": box.type
            })
            self._boxes.append(promoter_boxes.copy())

    @property
    def transcription_start_site(self):
        return self._transcription_start_site

    @transcription_start_site.setter
    def transcription_start_site(self, promoter):
        self._transcription_start_site = {}
        if promoter.transcription_start_site:
            trans_ss = {
                "leftEndPosition": promoter.transcription_start_site.left_end_position,
                "rightEndPosition": promoter.transcription_start_site.right_end_position,
                "range": promoter.transcription_start_site.range,
                # TODO: Type doesn't appear
                "type": promoter.transcription_start_site.type
            }
            self._transcription_start_site = trans_ss

    @property
    def regulator_binding_sites(self):
        return self._regulator_binding_sites

    @regulator_binding_sites.setter
    def regulator_binding_sites(self, promoter_id):
        self._regulator_binding_sites = []
        tf_binding_sites_dict = RegulatoryBindingSites(promoter_id)
        self._regulator_binding_sites = tf_binding_sites_dict.to_dict()

    def to_dict(self):
        citations = self.citations
        additive_evs = AdditiveEvidences(citations)
        promoter_dict = {
            "_id": self.promoter.id,
            "citations": citations,
            "bindsSigmaFactor": self.binds_sigma_factor,
            "name": self.promoter.name,
            "note": self.formatted_note,
            "boxes": self.boxes,
            "score": self.promoter.score,
            "sequence": self.promoter.sequence,
            "synonyms": self.promoter.synonyms,
            "regulatorBindingSites": self.regulator_binding_sites,
            "transcriptionStartSite": self.transcription_start_site,
            "additiveEvidences": additive_evs.to_dict(),
            "confidenceLevel": additive_evs.get_confidence_level()
        }
        return promoter_dict
