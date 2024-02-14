from mongoengine.errors import DoesNotExist

from src.datamarts.domain.general.biological_base import BiologicalBase
import multigenomic_api

from src.datamarts.domain.regulon_datamart.regulator import Regulator
from src.datamarts.domain.regulon_datamart.regulates import Regulates
from src.datamarts.domain.regulon_datamart.terms import Terms
from src.datamarts.domain.regulon_datamart.regulatory_interactions import RegulatoryInteractions
from src.datamarts.domain.regulon_datamart.summary import Summary

from src.datamarts.domain.general.remove_items import remove_empty_items


class RegulonDatamarts:

    @property
    def objects(self):
        regulator_objects = multigenomic_api.regulators.get_all()
        for regulator in regulator_objects:
            print(regulator.id)
            if regulator.regulation_type:
                if regulator.regulation_type == "Transcription-Factor-Binding" or regulator.id == "RDBECOLITFC00039":
                    regulator = multigenomic_api.transcription_factors.find_by_id(regulator.id)
                    if regulator.id == "RDBECOLITFC00039":
                        regulator["regulation_type"] = "Allosteric-Regulation-of-RNAP"
                    else:
                        regulator["regulation_type"] = "Transcription-Factor-Binding"
                regulon_datamart = RegulonDatamarts.RegulonDatamart(regulator)
                yield regulon_datamart
        del regulator_objects

    class RegulonDatamart:

        def __init__(self, regulator):
            self.id = regulator.id
            self.regulator = regulator
            self.terms = regulator
            self.regulates = regulator
            self.regulatory_interactions = regulator
            # self.alignmentMatrix = regulator.id
            # self.evolutionaryConservation = regulator.id
            self.organism = regulator.organisms_id
            self.summary = [self.regulates, self._regulatory_interactions]

        @property
        def regulator(self):
            return self._regulator

        @regulator.setter
        def regulator(self, regulator):
            regulator = Regulator(regulator)
            self._regulator = regulator.to_dict()

        @property
        def terms(self):
            return self._terms

        @terms.setter
        def terms(self, regulator):
            self._terms = []
            terms = []
            if regulator.regulation_type != "Allosteric-Regulation-of-RNAP" or regulator.id == "RDBECOLITFC00039":
                terms = Terms(regulator).to_dict()
            self._terms = terms

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
        def regulatory_interactions(self, regulator):
            self._regulatory_interactions = []
            reg_complex = None
            product = None
            if regulator.regulation_type == "Transcription-Factor-Binding" or regulator.id == "RDBECOLITFC00039":
                try:
                    reg_complex = multigenomic_api.regulatory_complexes.find_by_name(regulator.name)
                except DoesNotExist:
                    pass
                if reg_complex is not None:
                    regulatory_interactions = multigenomic_api.regulatory_interactions.find_by_regulator_id(reg_complex.id)
                    self._regulatory_interactions = get_ri_objects(regulatory_interactions, self._regulatory_interactions)
                try:
                    product = multigenomic_api.products.find_by_name(regulator.name)
                except DoesNotExist:
                    pass
                if product is not None:
                    regulatory_interactions = multigenomic_api.regulatory_interactions.find_by_regulator_id(product.id)
                    self._regulatory_interactions = get_ri_objects(regulatory_interactions, self._regulatory_interactions)
                for active_conformation in regulator.active_conformations:
                    regulatory_interactions = multigenomic_api.regulatory_interactions.find_by_regulator_id(active_conformation.id)
                    self._regulatory_interactions = get_ri_objects(regulatory_interactions, self._regulatory_interactions)
            regulatory_interactions = multigenomic_api.regulatory_interactions.find_by_regulator_id(regulator.id)
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
                "regulator": self.regulator,
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

