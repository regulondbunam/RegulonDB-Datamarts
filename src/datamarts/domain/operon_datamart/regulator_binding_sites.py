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
        ri_groups = {}
        for ri in ris:
            # Getting the regulator data
            regulator = get_abbr_regulator_name(ri.regulator)
            # Creating a unique key for the id and name combination
            key = (regulator["name"], regulator["_id"], ri["function"])
            # If key doesn't exist in dict, its initilazed as dict with empty ri list
            if key not in ri_groups:
                ri_groups[key] = {
                    "regulator": regulator,
                    "function": ri["function"],
                    "mechanism": ri.mechanism,
                    "regulatoryInteractions": []
                }
            # Adding the RI item to the list in the corresponding key
            reg_sites_dict = {}
            if ri.regulatory_sites_id:
                reg_site = multigenomic_api.regulatory_sites.find_by_id(ri.regulatory_sites_id)
                reg_sites_dict = RegulatorySites(reg_site).to_dict()
            reg_int = RegulatoryInteractions(ri, reg_sites_dict).to_dict()
            ri_groups[key]["regulatoryInteractions"].append(reg_int)
        return list(ri_groups.values())


def get_abbr_regulator_name(regulator):
    reg_id = regulator.id
    reg_name = regulator.name

    if regulator.type == "product":
        reg = multigenomic_api.products.find_by_id(regulator.id)
        if reg.abbreviated_name:
            reg_name = reg.abbreviated_name

    tf = multigenomic_api.transcription_factors.find_tf_id_by_conformation_id(regulator.id)
    if tf is None:
        tf = multigenomic_api.transcription_factors.find_by_name(regulator.name)
        if tf:
            reg_id = tf.id
            reg_name = tf.abbreviated_name
    else:
        reg_id = tf.id
        reg_name = tf.abbreviated_name
    return {
        "_id": reg_id,
        "name": reg_name
    }


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
