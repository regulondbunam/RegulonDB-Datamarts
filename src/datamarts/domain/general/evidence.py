class Evidence:

    def __init__(self, evidence, _type=None):
        self.evidence = evidence
        self.type = _type

    @property
    def type(self):
        return self._type

    @type.setter
    def type(self, _type):
        self._type = _type

    def to_dict(self):
        evidence = {
            "type": self.evidence.type,
            "code": self.evidence.code,
            "id": self.evidence.id,
            "name": self.evidence.name,
            "additiveEvidenceCodeRule": self.evidence.cv_code_rule
        }
        return evidence
