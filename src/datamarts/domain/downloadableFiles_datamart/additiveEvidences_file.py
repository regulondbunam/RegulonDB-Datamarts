import multigenomic_api
from datetime import datetime
import json

def load_rules():
    with open("./lib/evidences_rules.json", 'r') as add_json_rules:
        rules = json.load(add_json_rules)
        return rules


class ObjectEvidences:

    @property
    def objects(self):
        all_rules = load_rules()
        for rule in all_rules:
            # print(ev_obj.id)
            add_ev_row = ObjectEvidences.ObjectEvidenceFile(rule)
            yield add_ev_row

    class ObjectEvidenceFile:
        def __init__(self, rule):
            self.rule = rule
            self.rules = rule["rules_values"]
            self.object_type = rule["pertainsTo"]
            self.rules_matched = rule["rules_values"]

        @property
        def rules(self):
            return self._rules

        @rules.setter
        def rules(self, rules):
            self._rules = "Additive Evidences("
            for code_rule in rules:
                self._rules += f"{code_rule}/"
            self._rules = self._rules[:-1]
            self._rules += ")"

        @property
        def object_type(self):
            return self._object_type

        @object_type.setter
        def object_type(self, object_type):
            self._object_type = object_type[0]

        @property
        def rules_matched(self):
            return self._rules_matched

        @rules_matched.setter
        def rules_matched(self, rules):
            self._rules_matched = ""
            for code_rule in rules:
                self._rules_matched += f"{code_rule}: ("
                evidences = multigenomic_api.evidences.find_one_by_code_rule(code_rule)
                for ev in evidences:
                    self._rules_matched += f"{ev.code} |"
                self._rules_matched = self._rules_matched[:-1]
                self._rules_matched += ") </br>"

        def to_row(self):
            return f"{self.rules}" \
                   f"\t{self.rules_matched}" \
                   f"\t{self.object_type}" \
                   f"\t{self.rule['type'] or ''}"


def all_evidences_rows():
    evidences = ObjectEvidences()
    evidences_content = ["1)rule_title\t2)matched_rules\t3)object_type\t4)confidence_level"]
    for ev in evidences.objects:
        evidences_content.append(ev.to_row())
    creation_date = datetime.now()
    evidence_doc = {
        "_id": "RDBECOLIDLF00014",
        "fileName": "AdditiveEvidenceSet",
        "title": "Complete Additive Evidences Set",
        "fileFormat": "rif-version 1",
        "license": "# RegulonDB is free for academic/noncommercial use\n# User is not entitled to change or erase data sets of the RegulonDB\n# database or to eliminate copyright notices from RegulonDB. Furthermore,\n# User is not entitled to expand RegulonDB or to integrate RegulonDB partly\n# or as a whole into other databank systems, without prior written consent\n# from CCG-UNAM.\n# Please check the license at https://regulondb.ccg.unam.mx/manual/aboutUs/terms-conditions",
        "citation": "# Heladia Salgado, Socorro Gama-Castro, et al., RegulonDB v12.0: a comprehensive resource of transcriptional regulation in E. coli K-12,\n# Nucleic Acids Research, 2023;, gkad1072, https://doi.org/10.1093/nar/gkad1072",
        "contact": {
            "person": "RegulonDB Team",
            "webPage": None,
            "email": "regulondb@ccg.unam.mx"
        },
        "version": "1.0",
        "creationDate": f"{creation_date.strftime('%m-%d-%Y')}",
        "columnsDetails": "# Columns:\n# (1) Title of the evidence\n# (2) Matched evidences with their respective rule\n# (3) Object(s) were this evidence appears (Promoter,Transcription Units,Regulatory Interactions,Transcription Factors)\n# (4) Confidence Level of the evidence",
        "content": " \n".join(evidences_content),
        "rdbVersion": "12.0",
        "description": "Additive Evidence Catalog.",
        "group": "EVIDENCE"
    }
    return evidence_doc