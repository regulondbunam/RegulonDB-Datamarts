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
        regulatory_int = \
            multigenomic_api.regulatory_interactions.find_regulatory_interactions_by_reg_entity_id(reg_entity)
        for ri in regulatory_int:
            if ri.regulator:
                if ri.mechanism == "Translation":
                    tf_ri_dict.setdefault(ri.regulator.id, []).append(ri)
                else:
                    trans_factors = \
                        multigenomic_api.transcription_factors.find_tf_id_by_conformation_id(ri.regulator.id)
                    for trans_factor in trans_factors:
                        tf_ri_dict.setdefault(trans_factor.id, []).append(ri)
        tf_binding_sites_dict = self.fill_tf_binding_sites_dict(tf_ri_dict)
        if tf_binding_sites_dict:
            self._tf_binding_sites = tf_binding_sites_dict

    def to_dict(self):
        return self._tf_binding_sites

    @staticmethod
    def fill_tf_binding_sites_dict(first_dict):
        transcription_factor_binding_sites = []
        mechanism = ""
        for regulator, ris in first_dict.items():
            repressor_ris = []
            activator_ris = []
            reg_sites_dict = {}
            if re.match(r"^RDB[A-Z\d_]{5}PDC[\dA-Z]{5}$", regulator):
                trans_factor = multigenomic_api.products.find_by_id(regulator)
            else:
                trans_factor = multigenomic_api.transcription_factors.find_by_id(regulator)
            for ri in ris:
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
                        "_id": trans_factor.id,
                        "name": trans_factor.name,
                        "function": "repressor"
                    },
                    "regulatoryInteractions": repressor_ris,
                    "function": "repressor",
                    "mechanism": mechanism
                })
            if len(activator_ris) != 0:
                transcription_factor_binding_sites.append({
                    "regulator": {
                        "_id": trans_factor.id,
                        "name": trans_factor.name,
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
            "absolutePosition": self.reg_site.absolute_position,
            "citations": self.citations,
            "leftEndPosition": self.reg_site.left_end_position,
            "length": self.reg_site.length,
            "note": self.formatted_note,
            "rightEndPosition": self.reg_site.right_end_position,
            "sequence": self.reg_site.sequence
        }
        return reg_sites_dict
