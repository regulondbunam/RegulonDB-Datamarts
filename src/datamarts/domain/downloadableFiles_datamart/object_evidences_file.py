import multigenomic_api
from datetime import datetime


class ObjectEvidences:

    @property
    def objects(self):
        evidences_objects = multigenomic_api.evidences.get_all()
        for ev_obj in evidences_objects:
            # print(ev_obj.id)
            ev_row = ObjectEvidences.ObjectEvidenceFile(ev_obj)
            yield ev_row
        del evidences_objects

    class ObjectEvidenceFile:
        def __init__(self, evidence):
            self.evidence = evidence
            self.code = evidence.code

        @property
        def code(self):
            return self._code

        @code.setter
        def code(self, code):
            if code.startswith("|") and code.endswith("|"):
                code = code[1:-1]
            self._code = code

        def get_rows(self):
            rows = []
            obj_type = []
            promoters = multigenomic_api.promoters.find_by_evidence_id(self.evidence.id)
            if len(promoters) > 0:
                obj_type.append("Promoter")
            trans_units = multigenomic_api.transcription_units.find_by_evidence_id(self.evidence.id)
            if len(trans_units) > 0:
                obj_type.append("Transcription Units")
            reg_ints = multigenomic_api.regulatory_interactions.find_by_evidence_id(self.evidence.id)
            reg_sites = multigenomic_api.regulatory_sites.find_by_evidence_id(self.evidence.id)
            if len(reg_ints) > 0 or len(reg_sites) > 0:
                obj_type.append("Regulatory Interactions")
            tfs = multigenomic_api.transcription_factors.find_by_evidence_id(self.evidence.id)
            if len(tfs) > 0:
                obj_type.append("Transcription Factors")
            obj_type = ";".join(obj_type)
            if obj_type != "":
                rows.append(self.to_row(obj_type))
            return rows

        def to_row(self, object_type):
            return f"{self.code}" \
                   f"\t{self.evidence.name}" \
                   f"\t{self.evidence.cv_code_rule or ''}" \
                   f"\t{self.evidence.type or ''}" \
                   f"\t{object_type}" \
                   f"\t{self.evidence.category or ''}"


def all_evidences_rows():
    evidences = ObjectEvidences()
    evidences_content = [
        "1)evidence_code\t2)evidence_name\t3)evidence_group\t4)confidence_level\t5)object_type\t6)evidence_category"]
    for ev in evidences.objects:
        evidences_content.extend(ev.get_rows())
    creation_date = datetime.now()
    evidence_doc = {
        "_id": "RDBECOLIDLF00013",
        "fileName": "EvidenceSet",
        "title": "Complete Object Evidence Set",
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
        "columnsDetails": "# Columns:\n# (1) Code associated to evidence\n# (2) Name of the Evidence\n# (3) Confidence Level of the evidence\n# (4) Group of the evidence (Code rule)\n# (5) Object(s) were this evidence appears (Promoter,Transcription Units,Regulatory Interactions,Transcription Factors)\n# (6) Category of the Evidence",
        "content": " \n".join(evidences_content),
        "rdbVersion": "12.0"
    }
    return evidence_doc
