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
            self.genes_names = tu.genes_ids
            self.promoter = tu.promoters_id
            self.tfIds = tu
            self.tf = tu
            self.terminator_ids = tu.terminators_ids
            self.sigma = self.promoter

        @property
        def operon(self):
            return self._operon

        @operon.setter
        def operon(self, operon_id):
            self._operon = ""
            if operon_id:
                operon = multigenomic_api.operons.find_by_id(operon_id)
                self._operon = operon

        @property
        def genes_names(self):
            return self._genes_names

        @genes_names.setter
        def genes_names(self, genes_ids):
            self._genes_names = ""
            if len(genes_ids) > 0:
                for gene_id in genes_ids:
                    gene = multigenomic_api.genes.find_by_id(gene_id)
                    self._genes_names += f"{gene.name};"

        @property
        def promoter(self):
            return self._promoter

        @promoter.setter
        def promoter(self, promoter_id):
            self._promoter = {}
            if promoter_id:
                promoter = multigenomic_api.promoters.find_by_id(promoter_id)
                self._promoter = promoter
            else:
                self.promoter["id"] = ""
                self.promoter["name"] = ""

        @property
        def tfIds(self):
            return self._tfIds

        @tfIds.setter
        def tfIds(self, tu):
            self._tfIds = []
            ris = multigenomic_api.regulatory_interactions.find_by_regulated_entity_ids([tu.id, tu.promoters_id])
            for ri in ris:
                regulators = multigenomic_api.regulators.find_by_conformation(ri.regulator.id)
                for regulator in regulators:
                    if regulator.type == "transcriptionFactor":
                        self._tfIds.append(regulator.id)
            self._tfIds = ",".join(self._tfIds)

        @property
        def tf(self):
            return self._tf

        @tf.setter
        def tf(self, tu):
            self._tf = []
            ris = multigenomic_api.regulatory_interactions.find_by_regulated_entity_ids([tu.id, tu.promoters_id])
            for ri in ris:
                regulators = multigenomic_api.regulators.find_by_conformation(ri.regulator.id)
                for regulator in regulators:
                    if regulator.type == "transcriptionFactor":
                        self._tf.append(regulator.abbreviated_name)
            self._tf = ",".join(self._tf)

        @property
        def terminator_ids(self):
            return self._terminator_ids

        @terminator_ids.setter
        def terminator_ids(self, terminators_ids):
            self._terminator_ids = ""
            if len(terminators_ids) > 0:
                self._terminator_ids = ",".join(terminators_ids)

        @property
        def sigma(self):
            return self._sigma

        @sigma.setter
        def sigma(self, promoter):
            self._sigma = ""
            if promoter:
                if hasattr(promoter, "binds_sigma_factor"):
                    if promoter["binds_sigma_factor"]:
                        sigma = multigenomic_api.sigma_factors.find_by_id(promoter.binds_sigma_factor.sigma_factors_id)
                        self._sigma = sigma.abbreviated_name

        def to_row(self):
            return f"{self.tu.id}" \
                   f"\t{self.tu.name}" \
                   f"\t{self.operon.id}" \
                   f"\t{self.operon.name}" \
                   f"\t{",".join(self.tu.genes_ids)}" \
                   f"\t{self.genes_names}" \
                   f"\t{self.promoter["id"]}" \
                   f"\t{self.promoter["name"]}" \
                   f"\t{self.terminator_ids}" \
                   f"\t{self.sigma}" \
                   f"\t{self.tfIds}" \
                   f"\t{self.tf}"


def all_tus_rows(rdb_version, citation):
    trans_units = TranscriptionUnit()
    tus_content = ["1)tuId\t2)tuName\t3)operonId\t4)operonName\t5)tuGenesIds\t6)tuGenesNames\t7)promoterId\t8)promoterName\t9)terminatorIds\t10)sigmaFactor\t11)tfsIds\t12)tfsNames"]
    for tu in trans_units.objects:
        tus_content.append(tu.to_row())
    creation_date = datetime.now()
    tus_doc = {
        "_id": "RDBECOLIDLF00022",
        "fileName": "TuPromOperonTFSet",
        "title": "Transcription Unit with promoter, operon and TF Set",
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
        "columnsDetails": "# Columns:\n# "
                          "(1) tuId. Transcription Unit identifier assigned by RegulonDB\n# "
                          "(2) tuName. Transcription unit name \n# "
                          "(3) operonId. Operon id containing the transcription uni \n# "
                          "(4) operonName. Operon name containing the transcription unit\n# "
                          "(5) tuGenesIds. Ids of the gene(s) contained in the transcription unit\n# "
                          "(6) tuGenesNames. Names of the gene(s) contained in the transcription unit\n# "
                          "(7) promoterId. Promoter Id\n# "
                          "(8) promoterName. Promoter Name\n# "
                          "(9) terminatorIds. Ids of the terminators associated to transcription unit\n#"
                          "(10) sigmaFactor. Sigma Factor associated to \n#"
                          "(11) tfsIds. transcription Factors ids that regulates the transcription Unit\n#"
                          "(12) tfsNames. transcription Factors names that regulates the transcription Unit\n#",
        "content": " \n".join(tus_content),
        "rdbVersion": rdb_version,
        "description": "Transcription units with information of operon, promoter, terminator and tfs.",
        "group": "OPERON STRUCTURE"
    }
    return tus_doc
