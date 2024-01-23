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
            self.seq = read_dna_full_seq_file()

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
                if self.promoter.transcription_start_site is not None:
                    promoter_tss = self.promoter.transcription_start_site.left_end_position
                else:
                    promoter_tss = ""
                if tu.terminators_ids:
                    for terminator_id in tu.terminators_ids:
                        terminator = multigenomic_api.terminators.find_by_id(terminator_id)
                        first_gene = get_first_gene(genes)
                        last_gene = get_last_gene(genes)
                        coordinates = get_coordinates(self.promoter.transcription_start_site, terminator.transcriptionTerminationSite, first_gene, last_gene, operon.strand)
                        five_utr_ccs = get_5_utr_coords(self.promoter.transcription_start_site, first_gene, last_gene, operon.strand)
                        three_utr_ccs = get_3_utr_coords(terminator.transcriptionTerminationSite, first_gene, last_gene, operon.strand)
                        row = f"{operon.name}" \
                              f"\t{tu.name}" \
                              f"\t{self.promoter.name}" \
                              f"\t{promoter_tss}" \
                              f"\t{operon.strand}" \
                              f"\t{first_gene['resume']}" \
                              f"\t{last_gene['resume']}" \
                              f"\t{terminator.class_}({terminator.transcriptionTerminationSite.left_end_position},{terminator.transcriptionTerminationSite.right_end_position})" \
                              f"\t{coordinates}" \
                              f"\t{five_utr_ccs}" \
                              f"\t{get_seq_from_coords(five_utr_ccs, self.seq, operon.strand)}" \
                              f"\t{three_utr_ccs}" \
                              f"\t{get_seq_from_coords(three_utr_ccs, self.seq, operon.strand)}"
                        response.append(row)
                else:
                    first_gene = get_first_gene(genes)
                    last_gene = get_last_gene(genes)
                    coordinates = get_coordinates(self.promoter.transcription_start_site, None, first_gene, last_gene, operon.strand)
                    five_utr_ccs = get_5_utr_coords(self.promoter.transcription_start_site, first_gene, last_gene, operon.strand)
                    row = f"{operon.name}" \
                          f"\t{tu.name}" \
                          f"\t{self.promoter.name}" \
                          f"\t{promoter_tss}" \
                          f"\t{operon.strand}" \
                          f"\t{first_gene['resume']}" \
                          f"\t{last_gene['resume']}" \
                          f"\t{''}" \
                          f"\t{coordinates}" \
                          f"\t{five_utr_ccs}" \
                          f"\t{get_seq_from_coords(five_utr_ccs, self.seq, operon.strand)}" \
                          f"\t{''}" \
                          f"\t{''}"
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
    first_gene["resume"] = f"{first_gene.name}({first_gene.left_end_position},{first_gene.right_end_position})"
    return first_gene


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
    last_gene["resume"] = f"{last_gene.name}({last_gene.left_end_position},{last_gene.right_end_position})"
    return last_gene


def get_coordinates(tss, tts, first_gene, last_gene, strand):
    coordinates = ""
    if tss is not None and tts is not None:
        if strand == "forward":
            coordinates = f"{tss.left_end_position}-{tts.right_end_position}"
        elif strand == "reverse":
            coordinates = f"{tts.left_end_position}-{tss.left_end_position}"
    elif tss is not None:
        if strand == "forward":
            coordinates = f"{tss.left_end_position}-{last_gene.right_end_position}"
        elif strand == "reverse":
            coordinates = f"{first_gene.right_end_position}-{tss.left_end_position}"
    elif tts is not None:
        if strand == "forward":
            coordinates = f"{tts.left_end_position}-{last_gene.right_end_position}"
        elif strand == "reverse":
            coordinates = f"{first_gene.right_end_position}-{tts.left_end_position}"
    else:
        if strand == "forward":
            coordinates = f"{first_gene.left_end_position}-{last_gene.right_end_position}"
        elif strand == "reverse":
            coordinates = f"{last_gene.left_end_position}-{first_gene.right_end_position}"
    return coordinates


def get_5_utr_coords(tss, first_gene, last_gene, strand):
    coordinates = ""
    if tss is not None:
        if strand == "forward":
            coordinates = f"{tss.left_end_position}-{first_gene.left_end_position}"
        elif strand == "reverse":
            coordinates = f"{last_gene.right_end_position}-{tss.left_end_position}"
    return coordinates


def get_3_utr_coords(tts, first_gene, last_gene, strand):
    coordinates = ""
    if tts is not None:
        if strand == "forward":
            coordinates = f"{first_gene.right_end_position}-{tts.right_end_position}"
        elif strand == "reverse":
            coordinates = f"{tts.left_end_position}-{last_gene.left_end_position}"
    return coordinates


def read_dna_full_seq_file():
    with open('docs/full-dna_sequence.txt', 'r') as archivo:
        sequence = archivo.read()
        return sequence.replace('\n', '')


def get_seq_from_coords(coords, seq, strand):
    ccs = coords.split("-")
    sequence = ""
    if len(ccs) > 1:
        sequence = seq[int(ccs[0])-1:int(ccs[1])-1]
        if strand == "reverse":
            sequence = reverse_complement(sequence)
    return sequence


def reverse_complement(sequence=None):
    if sequence:
        alt_map = {'ins': '0'}
        complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'a': 't', 'c': 'g', 'g': 'c', 't': 'a'}
        for k, v in alt_map.items():
            sequence = sequence.replace(k, v)
        bases = list(sequence)
        bases = reversed([complement.get(base, base) for base in bases])
        bases = ''.join(bases)
        for k, v in alt_map.items():
            bases = bases.replace(v, k)
        return bases


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
        "citation": "# Heladia Salgado, Socorro Gama-Castro, et al., RegulonDB v12.0: a comprehensive resource of transcriptional regulation in E. coli K-12,\n# Nucleic Acids Research, 2023;, gkad1072, https://doi.org/10.1093/nar/gkad1072",
        "contact": {
            "person": "RegulonDB Team",
            "webPage": None,
            "email": "regulondb@ccg.unam.mx"
        },
        "version": "1.0",
        "creationDate": f"{creation_date.strftime('%m-%d-%Y')}",
        "columnsDetails": "# Columns:\n# (1) Operon Name\n# (2) Transcription Unit Name (TU)\n# (3) Promoter Name\n# (4) Transcription start site (+1)\n# (5) DNA strand  of TU\n# (6) First Gene Name of TU [Gene Position left, Gene Position right]\n# (7) Last Gene Name of TU [Gene Position left, Gene Position right]\n# (8) Terminator Type [Terminator Position left, Terminator Position right]\n# (9) Coordinates of UTR [+1 OR Start of the first gene] to [End of the terminator OR End of the last gene]*\n# (10) Coordinates 5' UTR\n# (11) 5' UTR Sequence	\n# (12) Coordinates 3' UTR\n# (13) 3' UTR Sequence\n#\n# * When the first one is NULL the other field is used",
        "content": " \n".join(operons_content),
        "rdbVersion": "12.0",
        "description": "Untranslated Transcription unit Regions.",
        "group": "OPERON STRUCTURE"
    }
    return operons_doc
