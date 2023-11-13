import multigenomic_api
from datetime import datetime


class Operon:

    @property
    def objects(self):
        operon_objects = multigenomic_api.operons.get_all()
        for operon_object in operon_objects:
            # print(operon_object.id)
            ri_row = Operon.OperonDatamart(operon_object)
            yield ri_row
        del operon_objects

    class OperonDatamart:
        def __init__(self, operon):
            operon_genes = get_genes(operon.id)
            self.operon = operon
            self.genes = operon_genes
            self.first_gene = operon_genes
            self.last_gene = operon_genes
            self.tus_evidences = operon.id
            self.evidences = self.tus_evidences
            self.confidence_level = self.tus_evidences

        @property
        def genes(self):
            return self._genes

        @genes.setter
        def genes(self, genes):
            self._genes = ""
            for gene in genes:
                self._genes += f"{gene.name};"
            if len(self._genes) > 0:
                self._genes = self._genes[:-1]

        @property
        def first_gene(self):
            return self._first_gene

        @first_gene.setter
        def first_gene(self, genes):
            self._first_gene = {
                "leftEndPosition": "",
            }
            if len(genes) == 1:
                self._first_gene = genes[0]
            elif len(genes) > 1:
                self._first_gene = genes[0]
                strand = genes[0].strand
                if strand == "forward":
                    for gene in genes:
                        if gene.left_end_position < self._first_gene.left_end_position:
                            self._first_gene = gene
                elif strand == "reverse":
                    for gene in genes:
                        if gene.right_end_position > self._first_gene.right_end_position:
                            self._first_gene = gene

        @property
        def last_gene(self):
            return self._last_gene

        @last_gene.setter
        def last_gene(self, genes):
            self._last_gene = {
                "rightEndPosition": "",
            }
            if len(genes) == 1:
                self._last_gene = genes[0]
            elif len(genes) > 1:
                self._last_gene = genes[0]
                strand = genes[0].strand
                if strand == "forward":
                    for gene in genes:
                        if gene.left_end_position > self._last_gene.left_end_position:
                            self._last_gene = gene
                elif strand == "reverse":
                    for gene in genes:
                        if gene.right_end_position < self._last_gene.right_end_position:
                            self._last_gene = gene

        @property
        def tus_evidences(self):
            return self._tus_evidences

        @tus_evidences.setter
        def tus_evidences(self, operon_id):
            self._tus_evidences = []
            tus = multigenomic_api.transcription_units.find_by_operon_id(operon_id)
            for tu in tus:
                if tu.citations:
                    for citation in tu.citations:
                        if citation not in self._tus_evidences:
                            self._tus_evidences.append(citation)

        @property
        def evidences(self):
            return self._evidences

        @evidences.setter
        def evidences(self, citations):
            self._evidences = []
            for citation in citations:
                if citation.evidences_id:
                    citation_dict = multigenomic_api.evidences.find_by_id(citation.evidences_id)
                    citation_item = f"[{citation_dict.code}:{citation_dict.type}]"
                    if citation_item not in self._evidences:
                        self._evidences.append(citation_item)
            self._evidences = "".join(self._evidences)

        @property
        def confidence_level(self):
            return self._confidence_level

        @confidence_level.setter
        def confidence_level(self, citations):
            self._confidence_level = ""
            for citation in citations:
                if citation.evidences_id:
                    citation_dict = multigenomic_api.evidences.find_by_id(citation.evidences_id)
                    if citation_dict.type == "C":
                        self._confidence_level = "C"
            for citation in citations:
                if citation.evidences_id:
                    citation_dict = multigenomic_api.evidences.find_by_id(citation.evidences_id)
                    if citation_dict.type == "S":
                        self._confidence_level = "S"
            for citation in citations:
                if citation.evidences_id:
                    citation_dict = multigenomic_api.evidences.find_by_id(citation.evidences_id)
                    if citation_dict.type == "W":
                        self._confidence_level = "W"

        def to_row(self):
            return f"{self.operon.id}" \
                   f"\t{self.operon.name}" \
                   f"\t{self.first_gene['left_end_position']}" \
                   f"\t{self.last_gene['right_end_position']}" \
                   f"\t{self.operon.strand}" \
                   f"\t{len(self.genes.split(';'))}" \
                   f"\t{self.genes}" \
                   f"\t{self.evidences}" \
                   f"\t{self.confidence_level}" \



def get_genes(operon_id):
    genes = []
    tus = multigenomic_api.transcription_units.find_by_operon_id(operon_id)
    for tu in tus:
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


def all_operons_rows():
    operons = Operon()
    operons_content = ["1)operonId\t2)operonName\t3)firstGeneLeftPos\t4)lastGeneRightPos\t5)strand\t6)numberOfGenes\t7)operonGenes\t8)operonEvidence\t9)confidenceLevel"]
    for operon in operons.objects:
        operons_content.append(operon.to_row())
    creation_date = datetime.now()
    operons_doc = {
        "_id": "RDBECOLIDLF00009",
        "fileName": "OperonSet",
        "title": "Complete Operons Set",
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
        "columnsDetails": "# Columns:\n# (1) Operon name\n# (2) First gene-position left\n# (3) Last gene-position right\n# (4) DNA strand where the operon is coded\n# (5) Number of genes contained in the operon\n# (6) Name or Blattner number of the gene(s) contained in the operon\n# (7) Evidence that support the existence of the operon's TUs\n# (8) Evidence confidence level (Confirmed, Strong, Weak)",
        "content": " \n".join(operons_content),
        "rdbVersion": "12.0"
    }
    return operons_doc
