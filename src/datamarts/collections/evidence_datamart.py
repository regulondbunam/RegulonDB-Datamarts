import multigenomic_api

from src.datamarts.domain.general.remove_items import remove_empty_items

class EvidenceDatamarts:

    @property
    def objects(self):
        evidence_objects = multigenomic_api.evidences.get_all()
        for ev_item in evidence_objects:
            print(ev_item.id)
            gene_datamart = EvidenceDatamarts.EvidenceDatamart(ev_item)
            yield gene_datamart
        del evidence_objects

    class EvidenceDatamart:

        def __init__(self, evidence):
            self.id = evidence.id
            self.evidence = evidence
            self.code = evidence.code
            self.pertains_to = evidence.id
            self.children = evidence

        @property
        def code(self):
            return self._code

        @code.setter
        def code(self, code):
            if code.startswith("|") and code.endswith("|"):
                code = code[1:-1]
            self._code = code

        @property
        def children(self):
            return self._children

        @children.setter
        def children(self, evidence):
            self._children = []
            children = multigenomic_api.evidences.find_by_head(evidence.id)
            for node in children:
                if node.id != evidence.id:
                    self._children.append(node.id)

        @property
        def pertains_to(self):
            return self._pertains_to

        @pertains_to.setter
        def pertains_to(self, ev_id):
            self._pertains_to = []
            promoters = multigenomic_api.promoters.find_by_evidence_id(self.evidence.id)
            if len(promoters) > 0:
                self._pertains_to.append("Promoter")
            trans_units = multigenomic_api.transcription_units.find_by_evidence_id(self.evidence.id)
            if len(trans_units) > 0:
                self._pertains_to.append("Transcription Units")
            reg_sites = multigenomic_api.regulatory_sites.find_by_evidence_id(self.evidence.id)
            if len(reg_sites) > 0:
                self._pertains_to.append("Regulatory Sites")
            reg_ints = multigenomic_api.regulatory_interactions.find_by_evidence_id(self.evidence.id)
            if len(reg_ints) > 0:
                self._pertains_to.append("Regulatory Interactions")
            tfs = multigenomic_api.transcription_factors.find_by_evidence_id(self.evidence.id)
            if len(tfs) > 0:
                self._pertains_to.append("Transcription Factors")
            terminators = multigenomic_api.terminators.find_by_evidence_id(self.evidence.id)
            if len(terminators) > 0:
                self._pertains_to.append("Terminators")
            regulators = multigenomic_api.regulators.find_by_evidence_id(self.evidence.id)
            if len(regulators) > 0:
                self._pertains_to.append("Regulators")
            continuants = multigenomic_api.regulatory_continuants.find_by_evidence_id(self.evidence.id)
            if len(continuants) > 0:
                self._pertains_to.append("Regulatory Continuants")

        def to_dict(self):
            return {
                "_id": self.id,
                "name": self.evidence.name,
                "code": self.code,
                "superClassOf": self.children,
                "category": self.evidence.category,
                "codeRule": self.evidence.cv_code_rule,
                "note": self.evidence.note,
                "pertainsTo": self.pertains_to,
                "type": self.evidence.type,
            }


def all_evidences_datamarts():
    evidences = EvidenceDatamarts()
    json_evs = []
    for evidence in evidences.objects:
        evidence_dict = remove_empty_items(evidence.to_dict().copy())
        json_evs.append(remove_empty_items(evidence_dict))
    return json_evs