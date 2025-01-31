import multigenomic_api
from datetime import datetime
from . import common_functions


class TranscriptionUnit:

    @property
    def objects(self):
        tu_objects = multigenomic_api.transcription_units.get_all()
        for tu_object in tu_objects:
            ri_row = TranscriptionUnit.TUDatamart(tu_object)
            yield ri_row
        del tu_objects

    class TUDatamart:
        def __init__(self, tu):
            self.tu = tu
            self.operon = tu.operons_id
            self.genes_names = tu.genes_ids
            self.promoter = tu.promoters_id
            self.promoter_pos = tu.promoters_id
            self.tfIds = tu
            self.tf = tu
            self.terminator_ids = tu.terminators_ids
            self.sigma = self.promoter
            self.genes_bnumber = tu.genes_ids
            self.first_gene_pos = tu
            self.last_gene_pos = tu
            self.terminator_pos = tu.terminators_ids
            self.dbxref = tu.external_cross_references[0].object_id

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
            genes = []
            if len(genes_ids) > 0:
                for gene_id in genes_ids:
                    gene = multigenomic_api.genes.find_by_id(gene_id)
                    genes.append(gene.name)
                self._genes_names = ";".join(genes)

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
        def promoter_pos(self):
            return self._promoter_pos

        @promoter_pos.setter
        def promoter_pos(self, promoter_id):
            self._promoter_pos = ""
            if promoter_id:
                promoter = multigenomic_api.promoters.find_by_id(promoter_id)
                if promoter.transcription_start_site:
                    self._promoter_pos = promoter.transcription_start_site.left_end_position

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

        @property
        def genes_bnumber(self):
            return self._genes_bnumber

        @genes_bnumber.setter
        def genes_bnumber(self, genes_ids):
            self._genes_bnumber = ""
            genes = []
            if len(genes_ids) > 0:
                for gene_id in genes_ids:
                    gene = multigenomic_api.genes.find_by_id(gene_id)
                    if gene.bnumber:
                        genes.append(gene.bnumber)
                self._genes_bnumber = ";".join(genes)

        @property
        def first_gene_pos(self):
            return self._first_gene_pos

        @first_gene_pos.setter
        def first_gene_pos(self, tu):
            self._first_gene_pos = ""
            if tu.operons_id:
                operon = multigenomic_api.operons.find_by_id(tu.operons_id)
                first_gene = get_gene_of_tu(tu.genes_ids, operon.strand, "first")
                if operon.strand == "forward":
                    self._first_gene_pos = first_gene.get("left_end_position")
                elif operon.strand == "reverse":
                    self._first_gene_pos = first_gene.get("right_end_position")

        @property
        def last_gene_pos(self):
            return self._last_gene_pos

        @last_gene_pos.setter
        def last_gene_pos(self, tu):
            self._last_gene_pos = ""
            if tu.operons_id:
                operon = multigenomic_api.operons.find_by_id(tu.operons_id)
                first_gene = get_gene_of_tu(tu.genes_ids, operon.strand, "last")
                if operon.strand == "forward":
                    self._last_gene_pos = first_gene.get("right_end_position")
                elif operon.strand == "reverse":
                    self._last_gene_pos = first_gene.get("left_end_position")

        @property
        def terminator_pos(self):
            return self._terminator_pos

        @terminator_pos.setter
        def terminator_pos(self, terminator_ids):
            self._terminator_pos = ""
            positions = []
            for terminator_id in terminator_ids:
                terminator = multigenomic_api.terminators.find_by_id(terminator_id)
                positions.append(f"{terminator.transcriptionTerminationSite.left_end_position}-{terminator.transcriptionTerminationSite.right_end_position}")
            self._terminator_pos = ";".join(positions)


        @property
        def dbxref(self):
            return self._dbxref

        @dbxref.setter
        def dbxref(self, ecocyc_id):
            self._dbxref = f"Ecocyc:{ecocyc_id}"

        def to_row(self):
            return f"{self.tu.id}" \
                   f"\t{self.tu.name}" \
                   f"\t{self.operon.id}" \
                   f"\t{self.operon.name}" \
                   f"\t{','.join(self.tu.genes_ids)}" \
                   f"\t{self.genes_names}" \
                   f"\t{self.promoter['id']}" \
                   f"\t{self.promoter['name']}" \
                   f"\t{self.promoter_pos}" \
                   f"\t{self.terminator_ids}" \
                   f"\t{self.sigma}" \
                   f"\t{self.tfIds}" \
                   f"\t{self.tf}" \
                   f"\t{self.genes_bnumber}" \
                   f"\t{self.first_gene_pos}" \
                   f"\t{self.last_gene_pos}" \
                   f"\t{self.terminator_pos}" \
                   f"\t{self.dbxref}"



def all_tus_rows(rdb_version, citation):
    trans_units = TranscriptionUnit()
    tus_content = ["1)tuId\t2)tuName\t3)operonId\t4)operonName\t5)tuGenesIds\t6)tuGenesNames\t7)promoterId\t8)promoterName\t9)promoterPos\t10)terminatorIds\t11)sigmaFactor\t12)tfsIds\t13)tfsNames\t14)tuGenesBnumber\t15)firstGenePos\t16)lastGenePos\t17)terminatorPositions\t18)DBXRef"]
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
                          "(9) promoterPos. Transcription start site of promoter"
                          "(10) terminatorIds. Ids of the terminators associated to transcription unit\n#"
                          "(11) sigmaFactor. Sigma Factor associated to \n#"
                          "(12) tfsIds. transcription Factors ids that regulates the transcription Unit\n#"
                          "(13) tfsNames. transcription Factors names that regulates the transcription Unit\n#"
                          "(14) tuGenesBnumber. bnumbers of the tu genes\n#"
                          "(15) firstGenePos. start position of the first tu gene\n#"
                          "(16) lastGenePos. end position of the first tu gene\n#"
                          "(17) terminatorPositions. positions of the terminator\n#"
                          "(18) DBXRef. id of the TU in Ecocyc\n#",
        "content": " \n".join(tus_content),
        "rdbVersion": rdb_version,
        "description": "Transcription units with information of operon, promoter, terminator and tfs.",
        "group": "OPERON STRUCTURE"
    }
    return tus_doc


def get_gene_of_tu(genes, strand, position):
    dict_genes = []
    selected_gene = None

    for gene in genes:
        gene_object = multigenomic_api.genes.find_by_id(gene)
        if gene_object.fragments:
            min_left_pos = min(gene_object.fragments, key=lambda x: x.left_end_position)
            max_right_pos = max(gene_object.fragments, key=lambda x: x.right_end_position)
            gene_object.left_end_position = min_left_pos.left_end_position
            gene_object.right_end_position = max_right_pos.right_end_position
        dict_genes.append({
            "id": gene_object.id,
            "name": gene_object.name,
            "left_end_position": gene_object.left_end_position,
            "right_end_position": gene_object.right_end_position
        })

    if len(dict_genes) > 0:
        if strand == "forward":
            if position == "first":
                selected_gene = min(dict_genes, key=lambda x: x["left_end_position"])
            elif position == "last":
                selected_gene = max(dict_genes, key=lambda x: x["right_end_position"])
        elif strand == "reverse":
            if position == "first":
                selected_gene = max(dict_genes, key=lambda x: x["right_end_position"])
            elif position == "last":
                selected_gene = min(dict_genes, key=lambda x: x["left_end_position"])

    return selected_gene

