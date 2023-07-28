import multigenomic_api
import re
from src.datamarts.domain.general.biological_base import BiologicalBase
from src.datamarts.domain.operon_datamart.reg_binding_sites.regulatory_interactions import RegulatoryInteractions


class RegulatoryBindingSites(BiologicalBase):
    def __init__(self, reg_entity):
        super().__init__([], [], None)
        self.tf_binding_sites = reg_entity

    @property
    def tf_binding_sites(self):
        return self._tf_binding_sites

    @tf_binding_sites.setter
    def tf_binding_sites(self, reg_entity):
        self._tf_binding_sites = []
        tf_ri_dict = {}
        '''
        Update this when sRNA is pushed on regulondbmultigenomic with mechanism in RI's
        '''
        regulatory_ints = \
            multigenomic_api.regulatory_interactions.find_regulatory_interactions_by_reg_entity_id(reg_entity)
        tf_binding_sites_dict = self.fill_tf_binding_sites_dict(regulatory_ints)
        if tf_binding_sites_dict:
            self._tf_binding_sites = tf_binding_sites_dict

    def to_dict(self):
        return self._tf_binding_sites

    @staticmethod
    def fill_tf_binding_sites_dict(ris):
        transcription_factor_binding_sites = []
        mechanism = ""
        for ri in ris:
            repressor_ris = []
            activator_ris = []
            reg_sites_dict = {}
            mechanism = ri.mechanism
            if ri.regulatory_sites_id:
                reg_site = multigenomic_api.regulatory_sites.find_by_id(ri.regulatory_sites_id)
                reg_sites_dict = RegulatorySites(reg_site).to_dict()
            reg_int = RegulatoryInteractions(ri, reg_sites_dict).to_dict()
            if ri.function == "repressor":
                repressor_ris.append(reg_int)
            elif ri.function == "activator":
                activator_ris.append(reg_int)
            if len(repressor_ris) != 0:
                transcription_factor_binding_sites.append({
                    "regulator": {
                        "_id": ri.regulator.id,
                        "name": ri.regulator.name,
                        "function": "repressor"
                    },
                    "regulatoryInteractions": repressor_ris,
                    "function": "repressor",
                    "mechanism": mechanism
                })
            if len(activator_ris) != 0:
                transcription_factor_binding_sites.append({
                    "regulator": {
                        "_id": ri.regulator.id,
                        "name": ri.regulator.name,
                        "function": "activator"
                    },
                    "regulatoryInteractions": activator_ris,
                    "function": "activator",
                    "mechanism": mechanism
                })
        return transcription_factor_binding_sites


class RegulatorySites(BiologicalBase):
    def __init__(self, reg_site):
        super().__init__([], reg_site.citations, reg_site.note)
        self.reg_site = reg_site

    def to_dict(self):
        reg_sites_dict = {
            "_id": self.reg_site.id,
            "centerEndPosition": self.reg_site.absolute_position,
            "citations": self.citations,
            "leftEndPosition": self.reg_site.left_end_position,
            "length": self.reg_site.length,
            "note": self.formatted_note,
            "rightEndPosition": self.reg_site.right_end_position,
            "sequence": self.reg_site.sequence
        }
        return reg_sites_dict
