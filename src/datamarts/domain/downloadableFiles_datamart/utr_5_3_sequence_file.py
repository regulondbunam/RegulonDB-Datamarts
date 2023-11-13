import multigenomic_api
from datetime import datetime


class UTR_Sequence:

    @property
    def objects(self):
        promoter_objects = multigenomic_api.promoters.get_all()
        for promoter in promoter_objects:
            # print(operon_object.id)
            promoter_row = UTR_Sequence.UTRSeq(promoter)
            yield promoter_row
        del promoter_objects

    class UTRSeq:
        def __init__(self, promoter):
            self.promoter = promoter
            self.trans_units = promoter.id

        @property
        def operon(self):
            return self._operon

        @operon.setter
        def operon(self, operon_id):
            self._operon = {
                "name": "",
                "strand": ""
            }
            if operon_id:
                self._operon = multigenomic_api.operons.find_by_id(operon_id)

        @property
        def trans_units(self):
            return self._trans_units

        @trans_units.setter
        def trans_units(self, promoter_id):
            self._trans_units = multigenomic_api.transcription_units.find_by_promoter_id(promoter_id)

        def to_row(self):
            response = []
            for tu in self.trans_units:
                coordinates = ""
                operon = multigenomic_api.operons.find_by_id(tu.operons_id)
                genes = get_tu_genes(tu)
                for terminator_id in tu.terminators.ids:
                    terminator = multigenomic_api.terminators.find_by_id(terminator_id)
                    if self.promoter.transcription_start_site is not None:
                        promoter_tss = self.promoter.transcription_start_site.left_end_position
                    else:
                        promoter_tss = ""

                    if self.promoter.transcription_start_site is not None and terminator.transcriptionTerminationSite is not None:
                        coordinates = f"{self.promoter.transcription_start_site.left_end_position}-{terminator.transcriptionTerminationSite.right_end_position}"

                    row = f"{operon.name}" \
                          f"\t{tu.name}" \
                          f"\t{self.promoter.name}" \
                          f"\t{promoter_tss}" \
                          f"\t{operon.strand}" \
                          f"\t{get_first_gene(genes)}" \
                          f"\t{get_last_gene(genes)}" \
                          f"\t{terminator.class_}({terminator.transcriptionTerminationSite.left_end_position},{terminator.transcriptionTerminationSite.right_end_position})" \
                          f"\t{coordinates}" \
                          f"\t{'Coordinates 5´ UTR'}" \
                          f"\t{'5´ UTR Sequence'}" \
                          f"\t{'Coordinates 3´ UTR'}" \
                          f"\t{'3´ UTR Sequence'}"
                    response.append(row)
            return response


def get_tu_genes(tu):
    genes = []
    if tu.genes_ids:
        for gene_id in tu.genes_ids:
            gene = multigenomic_api.genes.find_by_id(gene_id)
            if gene not in genes:
                genes.append(gene)
    genes = find_fragments(genes)
    return genes


def find_fragments(genes):
    for gene in genes:
        if gene.fragments:
            first_fragment = gene.fragments[0]
            for fragment in gene.fragments:
                if first_fragment.left_end_position < fragment.left_end_position:
                    first_fragment = fragment
            last_fragment = gene.fragments[0]
            for fragment in gene.fragments:
                if last_fragment.right_end_position < fragment.right_end_position:
                    last_fragment = fragment
            gene.left_end_position = first_fragment.left_end_position
            gene.right_end_position = last_fragment.right_end_position
    return genes


def get_first_gene(genes):
    first_gene = {
        "name": "",
        "rightEndPosition": "",
        "leftEndPosition": "",
    }
    if len(genes) == 1:
        first_gene = genes[0]
    elif len(genes) > 1:
        first_gene = genes[0]
        strand = genes[0].strand
        if strand == "forward":
            for gene in genes:
                if gene.left_end_position < first_gene.left_end_position:
                    first_gene = gene
        elif strand == "reverse":
            for gene in genes:
                if gene.right_end_position > first_gene.right_end_position:
                    first_gene = gene
    return f"{first_gene.name}({first_gene.left_end_position},{first_gene.right_end_position})"


def get_last_gene(genes):
    last_gene = {
        "name": "",
        "rightEndPosition": "",
        "leftEndPosition": "",
    }
    if len(genes) == 1:
        last_gene = genes[0]
    elif len(genes) > 1:
        last_gene = genes[0]
        strand = genes[0].strand
        if strand == "forward":
            for gene in genes:
                if gene.left_end_position > last_gene.left_end_position:
                    last_gene = gene
        elif strand == "reverse":
            for gene in genes:
                if gene.right_end_position < last_gene.right_end_position:
                    last_gene = gene
    return f"{last_gene.name}({last_gene.left_end_position},{last_gene.right_end_position})"


def all_utr_rows():
    operons = UTR_Sequence()
    operons_content = [
        "1)operonName\t2)tuName\t3)promoterName\t4)tss\t5)strand\t6)firstGene\t7)lastGene\t8)terminatorType\t9)utrCoordinates\t10)utrCoordinates5\t11)utrSequence5\t12)utrCoordinates3\t13)utrSequence3"]
    for operon in operons.objects:
        operons_content.extend(operon.to_row())
    creation_date = datetime.now()
    operons_doc = {
        "_id": "RDBECOLIDLF00020",
        "fileName": "UTR_5_3_sequence",
        "title": "Complete UTR_5_3_sequence Set",
        "fileFormat": "rif-version 1",
        "license": "# RegulonDB is free for academic/noncommercial use\n# User is not entitled to change or erase data sets of the RegulonDB\n# database or to eliminate copyright notices from RegulonDB. Furthermore,\n# User is not entitled to expand RegulonDB or to integrate RegulonDB partly\n# or as a whole into other databank systems, without prior written consent\n# from CCG-UNAM.\n# Please check the license at https://regulondb.ccg.unam.mx/manual/aboutUs/terms-conditions",
        "citation": "Salgado H., Gama-Castro S. et al (2023). RegulonDB 12.0: A Comprehensive resource of transcriptional regulation in E. coli K-12",
        "contact": {
            "person": "RegulonDB Team",
            "webPage": None,
            "email": "regulondb@ccg.unam.mx"
        },
        "version": "1.0",
        "creationDate": f"{creation_date.strftime('%m-%d-%Y')}",
        "columnsDetails": "# Columns:\n# (1) Operon Name\n# (2) Transcription Unit Name (TU)\n# (3) Promoter Name\n# (4) Transcription start site (+1)\n# (5) DNA strand  of TU\n# (6) First Gene Name of TU [Gene Position left, Gene Position right]\n# (7) Last Gene Name of TU [Gene Position left, Gene Position right]\n# (8) Terminator Type [Terminator Position left, Terminator Position right]\n# (9) Coordinates of UTR [+1 OR Start of the first gene] to [End of the terminator OR End of the last gene]*\n# (10) Coordinates 5' UTR\n# (11) 5' UTR Sequence	\n# (12) Coordinates 3' UTR\n# (13) 3' UTR Sequence\n#\n# * When the first one is NULL the other field is used",
        "content": " \n".join(operons_content),
        "rdbVersion": "12.0"
    }
    return operons_doc
