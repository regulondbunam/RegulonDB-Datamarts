import multigenomic_api
from datetime import datetime


class RegulatoryInteractions:

    @property
    def objects(self):
        ri_objects = multigenomic_api.regulatory_interactions.get_all()
        for ri_object in ri_objects:
            print(ri_object.id)
            if ri_object.regulator.id != "RDBECOLIPDC00432" and ri_object.regulator.id != "RDBECOLIRCC00094":
                ri_row = RegulatoryInteractions.RIDatamart(ri_object)
                yield ri_row
        del ri_objects

    class RIDatamart:
        def __init__(self, ri):
            self.ri = ri
            self.type = ri
            self.transcription_factor = ri.regulator
            self.ri_pmids = ri.citations

        @property
        def type(self):
            return self._type

        @type.setter
        def type(self, ri):
            self._type = ""
            tf = multigenomic_api.transcription_factors.find_tf_id_by_conformation_id(ri.regulator.id)
            if tf is None:
                tf = multigenomic_api.transcription_factors.find_by_name(ri.regulator.name)
            if tf:
                regulator_type = "tf"
            elif ri.regulator.type == "product":
                regulator_type = "srna"
            elif ri.regulator.type == "regulatoryContinuant":
                regulator_type = "compound"
            else:
                regulator_type = "regulator"

            if ri.regulated_entity.type == "gene":
                self._type = f"{regulator_type}-gene"
            elif ri.regulated_entity.type == "promoter":
                self._type = f"{regulator_type}-promoter"
            if ri.regulated_entity.type == "transcriptionUnit":
                self._type = f"{regulator_type}-tu"

        @property
        def transcription_factor(self):
            return self._transcription_factor

        @transcription_factor.setter
        def transcription_factor(self, regulator):
            self._transcription_factor = {
                "id": "",
                "abbreviated_name": ""
            }
            tf = multigenomic_api.transcription_factors.find_tf_id_by_conformation_id(regulator.id)
            if tf is None:
                tf = multigenomic_api.transcription_factors.find_by_name(regulator.name)
                if tf is not None:
                    self._transcription_factor = tf
                else:
                    if regulator.type == "product":
                        product = multigenomic_api.products.find_by_id(regulator.id)
                        abb_name = product.abbreviated_name
                    else:
                        abb_name = regulator.name
                    self._transcription_factor = {
                        "id": regulator.id,
                        "abbreviated_name": abb_name
                    }
            else:
                self._transcription_factor = {
                    "id": tf.id,
                    "abbreviated_name": tf.abbreviated_name or tf.name
                }

        @property
        def ri_pmids(self):
            return self._ri_pmids

        @ri_pmids.setter
        def ri_pmids(self, citations):
            self._ri_pmids = []
            for citation in citations:
                if citation.publications_id:
                    publication = multigenomic_api.publications.find_by_id(citation.publications_id)
                    if publication.pmid not in self._ri_pmids:
                        self._ri_pmids.append(publication.pmid)

        def to_row(self):
            response = []
            for pmid in self.ri_pmids:
                row = f"{self.transcription_factor['abbreviated_name']}" \
                      f"\t{self.ri.regulated_entity.type}" \
                      f"\t{self.ri.regulated_entity.name}" \
                      f"\t{self.ri.function or ''}" \
                      f"\t{pmid}"
                response.append(row)
            return response


def all_ris_rows():
    ris = RegulatoryInteractions()
    ris_content = []
    ris_content.append("1)TF_Name\t2)Target_Type\t3)Target_Name\t4)Effect\t5)REFERENCE")
    for ri in ris.objects:
        if ri.type == "tf-gene" or ri.type == "tf-promoter" or ri.type == "tf-tu":
            ris_content.extend(ri.to_row())
    final_list = []
    for ri in ris_content:
        if ri not in final_list:
            final_list.append(ri)
    creation_date = datetime.now()
    ri_doc = {
        "_id": "RDBECOLIDLF00099",
        "fileName": "TF-RISet-Soco_internal",
        "title": "Complete TF RIs Set for Soco",
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
        "columnsDetails": "# Columns:\n# (1) riId. Regulatory interaction (RI) identifier assigned by RegulonDB\n# (2) riType. Regulatory interaction type [tf-promoter, tf-tu, tf-gene,srna-promoter, srna-tu, srna-gene, compound-promoter, compound-tu, compound-gene]\n# (3) tfId. Transcription Factor (TF) identifier assigned by RegulonDB\n# (4) regulatorName. Regulator name\n# (5) cnfName. regulator active conformation name\n# (6) tfrsId. regulator regulatory site (TFRS) identifier assigned by RegulonDB\n# (7) tfrsLeft. TFRS left end position in the genome\n# (8) tfrsRight. TFRS right end position in the genome\n# (9) strand. DNA strand where the TFRS is located\n# (10) tfrsSeq. TFRS sequence (upper case)\n# (11) riFunction. Gene expression effect caused by the TF bound to the TFRS\n# (12) promoterId. Promoter Identifier assigned by RegulonDB\n# (13) promoterName. Promoter name\n# (14) tss. Transcription start site (+1) position in the genome\n# (15) sigmaF. Sigma Factor that recognize the promoter\n# (16) DistToPm. Relative distance from the center position of TFRS to the Transcription Start Site\n# (17) firstGene. first transcribed gene name\n# (18) tfrsDistTo1Gene. Relative distance from center position of TFRS to the start of first gene\n# (19) targetTuOrGene. Transcription unit or gene (id:name) regulated by the TF\n# (20) confidenceLevel. RI confidence level (Values: Confirmed, Strong, Weak)\n# (21) tfrsEvidence. Evidence that supports the existence of the TFRS [EvidenceCode|EvidenceType(C:confirmed S:strong W:weak)]]\n# (22) riEvidence. Evidence that supports the RI function [EvidenceCode|EvidenceType(C:confirmed S:strong W:weak)]\n# (23) addEvidence. Additive Evidence [CV(EvidenceCode1/EvidenceCodeN)|Confidence Level]\n# (24) riEvTech. Evidence related to the type of technology used to determine the RI\n# (25) riEvCategory. Evidence  were categorized in classical, ht or non-experimental",
        "content": "\n".join(final_list),
        "rdbVersion": "12.0",
        "description": "Regulatory interactions between a Transcription Factor and its regulated entity - promoter, gene, or Transcriptional Units (subset of RISet).",
        "group": "REGULATORY INTERACTION"
    }
    return ri_doc
