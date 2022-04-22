import multigenomic_api


class Regulation:

    def __init__(self, operon):
        self.operon = operon
        self.regulators = self.operon
        self.promoters = self.operon
        self.regulatory_interactions = self.operon

    @property
    def operon(self):
        return self._operon

    @operon.setter
    def operon(self, operon):
        operon = Operon(operon)
        self._operon = operon

    @property
    def regulators(self):
        return self._regulators

    @regulators.setter
    def regulators(self, operon):
        self._regulators = operon.regulators

    @property
    def promoters(self):
        return self._promoters

    @promoters.setter
    def promoters(self, operon):
        self._promoters = operon.promoters

    @property
    def regulatory_interactions(self):
        return self._regulatory_interactions

    @regulatory_interactions.setter
    def regulatory_interactions(self, operon):
        self._regulatory_interactions = operon.regulatory_interactions

    def to_dict(self):
        regulation = {
            "operon": self.operon.to_dict(),
            "regulators": self.regulators,
            "statistics": {
                "regulators": len(self.regulators),
                "promoters": len(self.promoters),
                "regulatoryInteractions": len(self.regulatory_interactions)
            }
        }
        return regulation


class Operon:

    def __init__(self, operon):
        self.promoters = []
        self.regulators = []
        self.operon = operon
        self.transcription_units = multigenomic_api.transcription_units.find_by_operon_id(operon.id)
        self.arrangement = self.transcription_units
        self.regulatory_interactions = self.transcription_units

    @property
    def arrangement(self):
        return self._arrangement

    @arrangement.setter
    def arrangement(self, transcription_units):
        self._arrangement = []
        for transcription_unit in transcription_units:
            arrangement = {
                'transcriptionUnit': {
                    "id": transcription_unit.id,
                    "name": transcription_unit.name
                },
                'promoters': self.get_promoters(transcription_unit.promoters_id),
                'regulators': self.get_regulators(transcription_unit.id, transcription_unit.promoters_id)
            }
            self._arrangement.append(arrangement.copy())

    @property
    def regulatory_interactions(self):
        return self._regulatory_interactions

    @regulatory_interactions.setter
    def regulatory_interactions(self, transcription_units):
        reg_entities = []
        for transcription_unit in transcription_units:
            promoter_id = transcription_unit.promoters_id
            reg_entities = [promoter_id] + [transcription_unit.id]
        reg_entities = set(reg_entities)
        self._regulatory_interactions = \
            multigenomic_api.regulatory_interactions.find_by_regulated_entity_ids(reg_entities)

    def get_promoters(self, promoter_id):
        promoters = []
        if promoter_id:
            promoter = multigenomic_api.promoters.find_by_id(promoter_id)
            new_promoter = {
                'id': promoter.id,
                'name': promoter.name
            }
            promoters.append(new_promoter.copy())

            if new_promoter not in self.promoters:
                self.promoters.append(new_promoter.copy())

        return promoters

    def get_regulators(self, transcription_unit_id, promoter_id):
        arrangement_regulators = []
        reg_entities = [promoter_id] + [transcription_unit_id]
        regulators = multigenomic_api.regulatory_interactions.find_regulators_by_regulated_entity_ids(reg_entities)
        for regulator in regulators:
            new_regulator = {
                "id": regulator.id,
                "name": regulator.name,
                "type": regulator.type,
                "function": regulator.function
            }
            arrangement_regulators.append(new_regulator.copy())

            if new_regulator not in self.regulators:
                self.regulators.append(new_regulator.copy())

        return arrangement_regulators

    def to_dict(self):
        operon = {
            "id": self.operon.id,
            "name": self.operon.name,
            "arrangement": self.arrangement
        }
        return operon
