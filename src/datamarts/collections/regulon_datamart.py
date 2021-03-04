import multigenomic_api
from src.datamarts.domain.regulon_datamart.regulator import Regulator
from src.datamarts.domain.regulon_datamart.regulates import Regulates
from src.datamarts.domain.regulon_datamart.terms import Terms
from src.datamarts.domain.regulon_datamart.regulatory_interactions import RegulatoryInteractions

class RegulonDatamarts:

    @property
    def objects(self):
        regulator_objects = multigenomic_api.transcription_factors.get_all()
        for regulator_object in regulator_objects[0:1]:
            regulon_datamart = RegulonDatamarts.RegulonDatamart(regulator_object)
            yield regulon_datamart

    class RegulonDatamart:

        def __init__(self, regulator):
            self.id = regulator.id
            self.regulator = regulator
            self.terms = regulator.products_ids
            self.regulates = regulator
            self.regulatory_interactions = regulator.products_ids
            self.alignmentMatrix = regulator.id
            self.evolutionaryConservation = regulator.id
            self.organism = regulator.organisms_id
            self.summary = regulator.id

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
        def terms(self, products_ids):
            terms = Terms(products_ids)
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
        def regulatory_interactions(self, products_ids):
            regulatory_interactions = RegulatoryInteractions(products_ids)
            self._regulatory_interactions = regulatory_interactions.to_dict()

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
                # "summary": []
            }
            return regulon_datamart

def all_regulon_datamarts():
    regulons = RegulonDatamarts()
    json_regulons = []
    for regulon in regulons.objects:
        json_regulons.append(regulon.to_dict().copy())
    return json_regulons

