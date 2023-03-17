import multigenomic_api
from src.datamarts.domain.operon_datamart.operon import Operon
from src.datamarts.domain.operon_datamart.transcription_units import TranscriptionUnit
from src.datamarts.domain.general.biological_base import BiologicalBase

from src.datamarts.domain.general.remove_items import remove_empty_items


class OperonDatamarts:

    @property
    def objects(self):
        operon_objects = multigenomic_api.operons.get_all()
        for operon_object in operon_objects:
            print(operon_object.id)
            operon_datamart = OperonDatamarts.OperonDatamart(operon_object)
            yield operon_datamart
        del operon_objects

    class OperonDatamart:

        def __init__(self, operon):
            self.id = operon.id
            self.operon = operon
            self.transcription_units = multigenomic_api.transcription_units.find_by_operon_id(operon.id)
            self.organism = operon.organisms_id

        @property
        def operon(self):
            return self._operon

        @operon.setter
        def operon(self, operon):
            operon = Operon(operon)
            self._operon = operon.to_dict()

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
        def transcription_units(self):
            return self._transcription_units

        @transcription_units.setter
        def transcription_units(self, transcription_units):
            self._transcription_units = []
            for transcription_unit in transcription_units:
                transcription_unit = TranscriptionUnit(transcription_unit)
                self._transcription_units.append(transcription_unit.to_dict().copy())

        def to_dict(self):
            operon_datamart = {
                "_id": self.id,
                "operon": self.operon,
                "transcriptionUnits": self.transcription_units,
                "organism": self.organism,
                "allCitations": BiologicalBase.get_all_citations()
                # TODO: Añadir relatedIds
            }
            return operon_datamart


def all_operon_datamarts():
    operons = OperonDatamarts()
    json_operon = []
    for operon in operons.objects:
        operon_dict = remove_empty_items(operon.to_dict().copy())
        json_operon.append(remove_empty_items(operon_dict))
    return json_operon
