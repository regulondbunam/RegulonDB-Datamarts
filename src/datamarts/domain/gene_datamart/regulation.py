import multigenomic_api


class Regulation:

    def __init__(self, operon, gene_id):
        self.operon = [operon, gene_id]
        self.regulators = self.operon
        self.promoters = self.operon
        self.regulatory_interactions = self.operon

    @property
    def operon(self):
        return self._operon

    @operon.setter
    def operon(self, operon_gene):
        operon = Operon(operon_gene)
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

    def __init__(self, operon_gene):
        operon = operon_gene[0]
        gene_id = operon_gene[1]
        self.promoters = []
        self.regulators = []
        self.operon = operon
        self.transcription_units = multigenomic_api.transcription_units.find_by_operon_id(operon.id)
        self.arrangement = [self.transcription_units, gene_id]
        self.regulatory_interactions = [self.transcription_units, gene_id]

    @property
    def arrangement(self):
        return self._arrangement

    @arrangement.setter
    def arrangement(self, regulated_entities):
        self._arrangement = []
        transcription_units = regulated_entities[0]
        gene_id = regulated_entities[1]
        for transcription_unit in transcription_units:
            arrangement = {
                'transcriptionUnit': {
                    "_id": transcription_unit.id,
                    "name": transcription_unit.name
                },
                'promoters': self.get_promoters(transcription_unit.promoters_id),
                'regulators': self.get_regulators(transcription_unit.id, transcription_unit.promoters_id, gene_id)
            }
            self._arrangement.append(arrangement.copy())

    @property
    def regulatory_interactions(self):
        return self._regulatory_interactions

    @regulatory_interactions.setter
    def regulatory_interactions(self, regulated_entities):
        transcription_units = regulated_entities[0]
        gene_id = regulated_entities[1]
        reg_entities = []
        for transcription_unit in transcription_units:
            promoter_id = transcription_unit.promoters_id
            reg_entities.append(promoter_id)
            reg_entities.append(transcription_unit.id)
        reg_entities.append(gene_id)
        reg_entities = set(reg_entities)
        ris = multigenomic_api.regulatory_interactions.find_by_regulated_entity_ids(reg_entities)
        self._regulatory_interactions = ris

    def get_promoters(self, promoter_id):
        promoters = []
        if promoter_id:
            promoter = multigenomic_api.promoters.find_by_id(promoter_id)
            new_promoter = {
                '_id': promoter.id,
                'name': promoter.name
            }
            promoters.append(new_promoter.copy())

            if new_promoter not in self.promoters:
                self.promoters.append(new_promoter.copy())

        return promoters

    def get_regulators(self, transcription_unit_id, promoter_id, gene_id):
        arrangement_regulators = []
        # TODO: Añadir genes y nombre del TF
        reg_entities = [promoter_id] + [transcription_unit_id] + [gene_id]
        regulators = multigenomic_api.regulatory_interactions.find_regulators_by_regulated_entity_ids(reg_entities)
        for regulator in regulators:
            trans_factor = multigenomic_api.transcription_factors.find_tf_id_by_conformation_id(regulator.id)
            for tf in trans_factor:
                tf = {
                    "_id": tf.id,
                    "name": tf.abbreviated_name,
                    "function": regulator.function
                }
                arrangement_regulators = identify_dual(arrangement_regulators, tf)
                self.regulators = identify_dual(self.regulators, tf)
        return arrangement_regulators

    def to_dict(self):
        operon = {
            "_id": self.operon.id,
            "name": self.operon.name,
            "arrangement": self.arrangement
        }
        return operon


def identify_dual(regulators, item):
    for reg in regulators:
        if reg["_id"] == item["_id"]:
            if reg['function'] != item['function']:
                reg["function"] = "dual"
                item["function"] = "dual"
    if item not in regulators:
        regulators.append(item.copy())
    return regulators
