from src.datamarts.domain.general.biological_base import BiologicalBase
import multigenomic_api

from src.datamarts.domain.regulon_datamart.transcription_factor import TranscriptionFactor
from src.datamarts.domain.regulon_datamart.regulates import Regulates
from src.datamarts.domain.regulon_datamart.terms import Terms
from src.datamarts.domain.regulon_datamart.regulatory_interactions import RegulatoryInteractions
from src.datamarts.domain.regulon_datamart.summary import Summary

from src.datamarts.domain.general.remove_items import remove_empty_items


class RegulonDatamarts:

    @property
    def objects(self):
        regulator_objects = multigenomic_api.transcription_factors.get_all()
        for regulator_object in regulator_objects:
            print(regulator_object.id)
            regulon_datamart = RegulonDatamarts.RegulonDatamart(regulator_object)
            yield regulon_datamart
        del regulator_objects

    class RegulonDatamart:

        def __init__(self, transcription_factor):
            self.id = transcription_factor.id
            self.transcription_factor = transcription_factor
            self.terms = transcription_factor
            self.regulates = transcription_factor
            self.regulatory_interactions = transcription_factor.active_conformations
            self.alignmentMatrix = transcription_factor.id
            self.evolutionaryConservation = transcription_factor.id
            self.organism = transcription_factor.organisms_id
            self.summary = [self.regulates, self._regulatory_interactions]

        @property
        def transcription_factor(self):
            return self._transcription_factor

        @transcription_factor.setter
        def transcription_factor(self, transcription_factor):
            transcription_factor = TranscriptionFactor(transcription_factor)
            self._transcription_factor = transcription_factor.to_dict()

        @property
        def terms(self):
            return self._terms

        @terms.setter
        def terms(self, transcription_factor):
            terms = Terms(transcription_factor)
            self._terms = terms.to_dict()

        @property
        def regulates(self):
            return self._regulates

        @regulates.setter
        def regulates(self, regulator):
            regulates = Regulates(regulator)
            self._regulates = regulates.to_dict()

        @property
        def regulatory_interactions(self):
            return self._regulatory_interactions

        @regulatory_interactions.setter
        def regulatory_interactions(self, tf_active_conformations):
            self._regulatory_interactions = []
            for active_conformation in tf_active_conformations:
                regulatory_interactions = multigenomic_api.regulatory_interactions.find_by_regulator_id(
                    active_conformation.id)
                self._regulatory_interactions = get_ri_objects(regulatory_interactions, self._regulatory_interactions)
                if active_conformation.type == "regulatoryComplex":
                    reg_complex = multigenomic_api.regulatory_complexes.find_by_id(active_conformation.id)
                    for reg_cont_id in reg_complex.regulatory_continuants_ids:
                        regulatory_interactions = multigenomic_api.regulatory_interactions.find_by_regulator_id(reg_cont_id)
                        self._regulatory_interactions = get_ri_objects(regulatory_interactions, self._regulatory_interactions)

        @property
        def organism(self):
            return self._organism

        @organism.setter
        def organism(self, organism_id):
            organism = multigenomic_api.organisms.find_by_id(organism_id)
            if organism:
                self._organism = {
                    "_id": organism.id,
                    "name": organism.name
                }
            else:
                self._organism = {
                    "_id": organism_id
                }

        @property
        def summary(self):
            return self._summary

        @summary.setter
        def summary(self, regulated_objects):
            self._summary = Summary(regulated_objects[0], regulated_objects[1]).get_summary()

        def to_dict(self):
            regulon_datamart = {
                "_id": self.id,
                "transcriptionFactor": self.transcription_factor,
                "terms": self.terms,
                "regulates": self.regulates,
                "regulatoryInteractions": self.regulatory_interactions,
                # "alignmentMatrix": [],
                # "evolutionaryConservation":[],
                "organism": self.organism,
                "summary": self.summary,
                "allCitations": BiologicalBase.get_all_citations()
            }
            return regulon_datamart


def all_regulon_datamarts():
    regulons = RegulonDatamarts()
    json_regulons = []
    for regulon in regulons.objects:
        regulon_dict = remove_empty_items(regulon.to_dict().copy())
        json_regulons.append(remove_empty_items(regulon_dict))
    return json_regulons


def get_ri_objects(reg_ints, ri_list):
    for ri in reg_ints:
        reg_int = RegulatoryInteractions(ri).to_dict()
        if reg_int not in ri_list:
            ri_list.append(reg_int)
    return ri_list
