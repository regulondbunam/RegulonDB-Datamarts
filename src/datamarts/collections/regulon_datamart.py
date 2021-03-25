import multigenomic_api
from src.datamarts.domain.regulon_datamart.transcription_factor import TranscriptionFactor
from src.datamarts.domain.regulon_datamart.regulates import Regulates
from src.datamarts.domain.regulon_datamart.terms import Terms
from src.datamarts.domain.regulon_datamart.regulatory_interactions import RegulatoryInteractions

class RegulonDatamarts:

    @property
    def objects(self):
        regulator_objects = multigenomic_api.transcription_factors.get_all()
        for regulator_object in regulator_objects[2:3]:
            regulon_datamart = RegulonDatamarts.RegulonDatamart(regulator_object)
            yield regulon_datamart

    class RegulonDatamart:

        def __init__(self, transcription_factor):
            self.id = transcription_factor.id
            self.transcription_factor = transcription_factor
            self.terms = transcription_factor.products_ids
            # self.regulates = transcription_factor
            self.regulatory_interactions = transcription_factor
            self.alignmentMatrix = transcription_factor.id
            self.evolutionaryConservation = transcription_factor.id
            self.organism = transcription_factor.organisms_id
            self.summary = transcription_factor.id

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
        def regulatory_interactions(self, transcription_factor):
            self._regulatory_interactions = []
            for product_id in transcription_factor.products_ids:
                regulatory_interactions = multigenomic_api.regulatory_interactions.find_by_regulator_id(product_id)
                for ri in regulatory_interactions:
                    self._regulatory_interactions.append(RegulatoryInteractions(ri).to_dict())

        def to_dict(self):
            regulon_datamart = {
                "_id": self.id,
                "transcriptionFactor": self.transcription_factor,
                "terms": self.terms,
                # "regulates": self.regulates,
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

