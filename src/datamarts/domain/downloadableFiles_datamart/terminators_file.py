import multigenomic_api
from datetime import datetime


class Terminator:

    @property
    def objects(self):
        terminator_objects = multigenomic_api.terminators.get_all()
        for terminator_obj in terminator_objects:
            # print(terminator_obj.id)
            ri_row = Terminator.TerminatorDatamart(terminator_obj)
            yield ri_row
        del terminator_objects

    class TerminatorDatamart:
        def __init__(self, terminator):
            self.terminator = terminator
            self.operon = terminator.id
            self.strand = terminator.id
            self.rel_tus = terminator.id
            self.pmids = terminator.citations
            self.terminator_evidences = terminator.citations

        @property
        def operon(self):
            return self._operon

        @operon.setter
        def operon(self, terminator_id):
            self._operon = []
            tus = multigenomic_api.transcription_units.find_by_terminator_id(terminator_id)
            for tu in tus:
                operon = multigenomic_api.operons.find_by_id(tu.operons_id)
                if f"{operon.name};" not in self._operon:
                    self._operon.append(f"{operon.name};")
            self._operon = "".join(self._operon)[:-1]

        @property
        def strand(self):
            return self._strand

        @strand.setter
        def strand(self, terminator_id):
            self._strand = ""
            tus = multigenomic_api.transcription_units.find_by_terminator_id(terminator_id)
            operon = multigenomic_api.operons.find_by_id(tus[0].operons_id)
            self._strand = operon.strand

        @property
        def rel_tus(self):
            return self._rel_tus

        @rel_tus.setter
        def rel_tus(self, terminator_id):
            self._rel_tus = []
            tus = multigenomic_api.transcription_units.find_by_terminator_id(terminator_id)
            for tu in tus:
                if tu.promoters_id:
                    promoter = multigenomic_api.promoters.find_by_id(tu.promoters_id)
                    self._rel_tus.append(f"{tu.name}:{promoter.name};")
                else:
                    self._rel_tus.append(f"{tu.name};")
            self._rel_tus = "".join(self._rel_tus)[:-1]

        @property
        def terminator_evidences(self):
            return self._terminator_evidences

        @terminator_evidences.setter
        def terminator_evidences(self, citations):
            self._terminator_evidences = []
            for citation in citations:
                if citation.evidences_id:
                    citation_dict = multigenomic_api.evidences.find_by_id(citation.evidences_id)
                    citation_item = f"[{citation_dict.code}:{citation_dict.type}]"
                    if citation_item not in self._terminator_evidences:
                        self._terminator_evidences.append(citation_item)
            self._terminator_evidences = "".join(self._terminator_evidences)

        @property
        def pmids(self):
            return self._pmids

        @pmids.setter
        def pmids(self, citations):
            self._pmids = []
            if len(citations) > 0:
                for citation in citations:
                    if citation.publications_id:
                        publication = multigenomic_api.publications.find_by_id(citation.publications_id)
                        if f"{publication.pmid}, " not in self._pmids:
                            self._pmids.append(f"{publication.pmid}, ")
            self._pmids = "".join(self._pmids)[:-2]

        def to_row(self):
            return f"{self.terminator.id}" \
                   f"\t{self.terminator.transcriptionTerminationSite.left_end_position}" \
                   f"\t{self.terminator.transcriptionTerminationSite.right_end_position}" \
                   f"\t{self.strand}" \
                   f"\t{self.terminator.sequence}" \
                   f"\t{self.rel_tus}" \
                   f"\t{self.terminator.class_}" \
                   f"\t{self.operon}" \
                   f"\t{self.pmids}" \
                   f"\t{self.terminator_evidences}" \
                   f"\t{self.terminator.confidence_level}" \



def all_terminators_rows():
    terminators = Terminator()
    terminators_content = ["1)terminatorId\t2)leftEndPos\t3)rightEndPos\t4)strand\t5)sequence\t6)relatedTus\t7)type\t8)operonName\t9)pmids\t10)evidences\t11)confidenceLevel"]
    for terminator in terminators.objects:
        terminators_content.append(terminator.to_row())
    creation_date = datetime.now()
    terminators_doc = {
        "_id": "RDBECOLIDLF00012",
        "fileName": "TerminatorSet",
        "title": "Complete Terminators Set",
        "fileFormat": "rif-version 1",
        "license": "# RegulonDB is free for academic/noncommercial use\n# User is not entitled to change or erase data sets of the RegulonDB\n# database or to eliminate copyright notices from RegulonDB. Furthermore,\n# User is not entitled to expand RegulonDB or to integrate RegulonDB partly\n# or as a whole into other databank systems, without prior written consent\n# from CCG-UNAM.\n# Please check the license at https://regulondb.ccg.unam.mx/manual/aboutUs/terms-conditions",
        "citation": "Heladia Salgado, Socorro Gama-Castro, et al., RegulonDB v12.0: a comprehensive resource of transcriptional regulation in E. coli K-12, Nucleic Acids Research, 2023;, gkad1072, https://doi.org/10.1093/nar/gkad1072",
        "contact": {
            "person": "RegulonDB Team",
            "webPage": None,
            "email": "regulondb@ccg.unam.mx"
        },
        "version": "1.0",
        "creationDate": f"{creation_date.strftime('%m-%d-%Y')}",
        "columnsDetails": "# Columns:\n# (1) Terminator identifier assigned by RegulonDB\n# (2) Terminator left end position in the genome\n# (3) Terminator right end position in the genome\n# (4) DNA strand where the terminator is located\n# (5) Terminator sequence\n# (6) Transcription unit(s) related to the terminator   *see note1\n# (7) Terminator type\n# (8) Operon name\n# (9) References (PubMed ID)\n# (10) Evidence that supports the existence of the terminator \n# (11) Evidence confidence level (Confirmed, Strong, Weak)\n# note1: If there are more than one TU with the same name, it is because the promoter associated to the TU is different",
        "content": " \n".join(terminators_content),
        "rdbVersion": "12.0",
        "description": "Terminator where transcription ends and RNAP unbinds from DNA",
        "group": "OPERON STRUCTURE"
    }
    return terminators_doc
