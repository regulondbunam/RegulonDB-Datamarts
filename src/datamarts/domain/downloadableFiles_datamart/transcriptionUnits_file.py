import multigenomic_api
from datetime import datetime
from . import common_functions


class TranscriptionUnit:

    @property
    def objects(self):
        tu_objects = multigenomic_api.transcription_units.get_all()
        for tu_object in tu_objects:
            # print(tu_object.id)
            ri_row = TranscriptionUnit.TUDatamart(tu_object)
            yield ri_row
        del tu_objects

    class TUDatamart:
        def __init__(self, tu):
            self.tu = tu
            self.operon = tu.operons_id
            self.genes = tu.genes_ids
            self.promoter = tu.promoters_id
            self.terminator_ids = tu.terminators_ids
            self.tu_evidences = tu.citations
            self.additive_evidences = tu.additive_evidences_ids
            self.pmids = tu

        @property
        def operon(self):
            return self._operon

        @operon.setter
        def operon(self, operon_id):
            self._operon = {
                "_id": "",
                "name": ""
            }
            if operon_id:
                self._operon = multigenomic_api.operons.find_by_id(operon_id)

        @property
        def genes(self):
            return self._genes

        @genes.setter
        def genes(self, genes_ids):
            self._genes = ""
            if len(genes_ids) > 0:
                for gene_id in genes_ids:
                    gene = multigenomic_api.genes.find_by_id(gene_id)
                    self._genes += f"{gene.name};"

        @property
        def promoter(self):
            return self._promoter

        @promoter.setter
        def promoter(self, promoter_id):
            self._promoter = ""
            if promoter_id:
                promoter = multigenomic_api.promoters.find_by_id(promoter_id)
                self._promoter = promoter.name

        @property
        def terminator_ids(self):
            return self._terminator_ids

        @terminator_ids.setter
        def terminator_ids(self, terminators_ids):
            self._terminator_ids = ""
            if len(terminators_ids) > 0:
                self._terminator_ids = ",".join(terminators_ids)

        @property
        def tu_evidences(self):
            return self._tu_evidences

        @tu_evidences.setter
        def tu_evidences(self, citations):
            self._tu_evidences = []
            for citation in citations:
                if citation.evidences_id:
                    citation_dict = multigenomic_api.evidences.find_by_id(citation.evidences_id)
                    citation_item = f"[{citation_dict.code}:{citation_dict.type}]"
                    if citation_item not in self._tu_evidences:
                        self._tu_evidences.append(citation_item)
            self._tu_evidences = "".join(self._tu_evidences)

        @property
        def additive_evidences(self):
            return self._additive_evidences

        @additive_evidences.setter
        def additive_evidences(self, additive_evs_ids):
            self._additive_evidences = ""
            if len(additive_evs_ids) > 0:
                self._additive_evidences = ""
                for additive_evs_id in additive_evs_ids:
                    additive_evidence_dict = multigenomic_api.additive_evidences.find_by_id(additive_evs_id)
                    self._additive_evidences += f"[{additive_evidence_dict.code}:{additive_evidence_dict.confidence_level}]"

        @property
        def pmids(self):
            return self._pmids

        @pmids.setter
        def pmids(self, promoter):
            self._pmids = common_functions.get_pmids(promoter)

        def to_row(self):
            return f"{self.tu.id}" \
                   f"\t{self.tu.name}" \
                   f"\t{self.operon.id}" \
                   f"\t{self.operon.name}" \
                   f"\t{self.genes}" \
                   f"\t{self.promoter}" \
                   f"\t{self.terminator_ids}" \
                   f"\t{self.tu_evidences}" \
                   f"\t{self.additive_evidences}" \
                   f"\t{self.tu.confidence_level}" \
                   f"\t{self.pmids}"


def all_tus_rows(rdb_version, citation):
    trans_units = TranscriptionUnit()
    tus_content = ["1)tuId\t2)tuName\t3)operonId\t4)operonName\t5)tuGenes\t6)pmName\t7)terminatorIds\t8)tuEvidence\t9)addEvidence\t10)confidenceLevel\t11)pmids"]
    for tu in trans_units.objects:
        tus_content.append(tu.to_row())
    creation_date = datetime.now()
    tus_doc = {
        "_id": "RDBECOLIDLF00010",
        "fileName": "TUSet",
        "title": "Complete Transcription Unit Set",
        "fileFormat": "rif-version 1",
        "license": "# RegulonDB is free for academic/noncommercial use\n# User is not entitled to change or erase data sets of the RegulonDB\n# database or to eliminate copyright notices from RegulonDB. Furthermore,\n# User is not entitled to expand RegulonDB or to integrate RegulonDB partly\n# or as a whole into other databank systems, without prior written consent\n# from CCG-UNAM.\n# Please check the license at https://regulondb.ccg.unam.mx/manual/aboutUs/terms-conditions",
        "citation": citation,
        "contact": {
            "person": "RegulonDB Team",
            "webPage": None,
            "email": "regulondb@ccg.unam.mx"
        },
        "version": "1.0",
        "creationDate": f"{creation_date.strftime('%m-%d-%Y')}",
        "columnsDetails": "# Columns:\n"
                          "# (1) tuId. Transcription Unit identifier assigned by RegulonDB\n"
                          "# (2) tuName. Transcription unit name \n"
                          "# (3) operonId. Operon Id containing the transcription unit\n"
                          "# (4) operonName. Operon name containing the transcription unit\n"
                          "# (5) tuGenes. Name of the gene(s) contained in the transcription unit\n"
                          "# (6) pmName. Promoter Name\n"
                          "# (7) terminatorIds. Ids of the terminators associated to transcription unit\n"
                          "# (8) tuEvidence. Evidence that supports the existence of the transcription unit\n"
                          "# (9) addEvidence. Additive Evidence [CV(EvidenceCode1/EvidenceCodeN)|Confidence Level]\n"
                          "# (10) confidenceLevel. TU confidence level (Values: Confirmed, Strong, Weak)\n"
                          "# (11) pmids associated to object",
        "content": " \n".join(tus_content),
        "rdbVersion": rdb_version,
        "description": "Transcription units with information of operon, promoter and terminator.",
        "group": "OPERON STRUCTURE"
    }
    return tus_doc
