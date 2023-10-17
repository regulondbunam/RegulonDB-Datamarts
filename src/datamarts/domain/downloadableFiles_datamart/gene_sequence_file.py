import multigenomic_api
from datetime import datetime
import re


class GeneSequence:

    @property
    def objects(self):
        gene_objects = multigenomic_api.genes.get_all()
        for gene_object in gene_objects:
            # print(gene_object.id)
            ri_row = GeneSequence.GeneSequenceDatamart(gene_object)
            yield ri_row
        del gene_objects

    class GeneSequenceDatamart:
        def __init__(self, gene):
            self.gene = gene
            self.gene_lep = gene
            self.gene_rep = gene
            self.product = gene.id
            self.productType = gene.id
            self.startCodon = gene.sequence
            self.stopCodon = gene.sequence
            self.otherDBsIds = gene.external_cross_references
            self.relatedBNumbers = gene.synonyms

        @property
        def product(self):
            return self._product

        @product.setter
        def product(self, gene_id):
            self._product = ""
            products = multigenomic_api.products.find_by_gene_id(gene_id)
            for product in products:
                self._product += f"{product.name};"
            self._product = self._product[:-1]

        @property
        def productType(self):
            return self._productType

        @productType.setter
        def productType(self, gene_id):
            self._productType = ""
            products = multigenomic_api.products.find_by_gene_id(gene_id)
            for product in products:
                self._productType += f"{product.type};"
            self._productType = self._productType[:-1]

        @property
        def gene_lep(self):
            return self._gene_lep

        @gene_lep.setter
        def gene_lep(self, gene):
            if gene.fragments:
                self._gene_lep = gene.fragments[0].left_end_position
                for fragment in gene.fragments:
                    if fragment.left_end_position < self._gene_lep:
                        self._gene_lep = fragment.left_end_position
            else:
                self._gene_lep = gene.left_end_position

        @property
        def gene_rep(self):
            return self._gene_rep

        @gene_rep.setter
        def gene_rep(self, gene):
            if gene.fragments:
                self._gene_rep = gene.fragments[0].right_end_position
                for fragment in gene.fragments:
                    if fragment.right_end_position > self._gene_rep:
                        self._gene_rep = fragment.right_end_position
            else:
                self._gene_rep = gene.right_end_position

        @property
        def relatedBNumbers(self):
            return self._relatedBNumbers

        @relatedBNumbers.setter
        def relatedBNumbers(self, synonyms):
            self._relatedBNumbers = ""
            bnumber_pattern = r"^b[0-9]{4}$"
            for synonym in synonyms:
                if re.match(bnumber_pattern, synonym):
                    self._relatedBNumbers += f"{synonym};"
            if len(self._relatedBNumbers) > 0:
                self._relatedBNumbers = self._relatedBNumbers[:-1]

        @property
        def startCodon(self):
            return self._startCodon

        @startCodon.setter
        def startCodon(self, sequence):
            self._startCodon = ""
            if sequence:
                self._startCodon = sequence[0:3]

        @property
        def stopCodon(self):
            return self._stopCodon

        @stopCodon.setter
        def stopCodon(self, sequence):
            self._stopCodon = ""
            if sequence:
                self._stopCodon = sequence[-3:]

        @property
        def otherDBsIds(self):
            return self._otherDBsIds

        @otherDBsIds.setter
        def otherDBsIds(self, ext_cross_ref):
            self._otherDBsIds = []
            for ext_ref in ext_cross_ref:
                if ext_ref.external_cross_references_id:
                    ext_ref_dict = multigenomic_api.external_cross_references.find_by_id(
                        ext_ref.external_cross_references_id)
                    ext_ref_item = f"[{ext_ref_dict.name}:{ext_ref.object_id}]"
                    if ext_ref_item not in self._otherDBsIds:
                        self._otherDBsIds.append(ext_ref_item)
            self._otherDBsIds = "".join(self._otherDBsIds)

        def to_row(self):
            return f"{self.gene.id}" \
                   f"\t{self.gene.name}({self.gene.bnumber})" \
                   f"\t{self.gene_lep}" \
                   f"\t{self.gene_rep}" \
                   f"\t{self.gene.strand}" \
                   f"\t{self.productType}" \
                   f"\t{self.product}" \
                   f"\t{self.startCodon}" \
                   f"\t{self.stopCodon}" \
                   f"\t{self.gene.sequence}" \
                   f"\t{self.relatedBNumbers}" \
                   f"\t{self.otherDBsIds}"


def all_gene_rows():
    genes = GeneSequence()
    genes_content = [
        "1)geneId\t2)geneName\t3)leftEndPos\t4)rightEndPos\t5)strand\t6)productType\t7)productName\t8)startCodon\t9)stopCodon\t10)sequence\t11)relatedBnumbers\t12)otherDbsIds"]
    for gene in genes.objects:
        genes_content.append(gene.to_row())
    creation_date = datetime.now()
    genes_doc = {
        "_id": "RDBECOLIDLF00008",
        "fileName": "Gene_sequence",
        "title": "Complete Gene Sequence Set",
        "fileFormat": "rif-version 1",
        "license": "RegulonDB is free for academic/noncommercial use\t\tUser is not entitled to change or erase data sets of the RegulonDB\tdatabase or to eliminate copyright notices from RegulonDB. Furthermore,\tUser is not entitled to expand RegulonDB or to integrate RegulonDB partly\tor as a whole into other databank systems, without prior written consent\tfrom CCG-UNAM.\t\tPlease check the license at http://regulondb.ccg.unam.mx/menu/download/full_version/terms_and_conditions.jsp",
        "citation": "Tierrafría, V. H. et al. (2022). RegulonDB 11.0: Comprehensive high-throughput datasets on transcriptional regulation in Escherichia coli K-12,\tMicrob Genom. 2022 May;8(5). doi: 10.1099/mgen.0.000833. PMID: 35584008. https://doi.org/10.1099/mgen.0.000833",
        "contact": {
            "person": "RegulonDB Team",
            "webPage": "http://regulondb.ccg.unam.mx/menu/about_regulondb/contact_us/index.jsp",
            "email": "regulondb@ccg.unam.mx"
        },
        "version": "1.0",
        "creationDate": f"{creation_date.strftime('%m-%d-%Y')}",
        "columnsDetails": "Columns:\n(1) Gene identifier assigned by RegulonDB\n(2) Gene name (bnumber)\n(3) Gene left end position in the genome\n(4) Gene right end position in the genome\n(5) DNA strand where the gene is coded\n(6) Product type\n(7) Product name\n(8) Start codon sequence\n(9) Stop codon sequence\n(10) Gene sequence\n(11) All bnumber related to gene\n(12) Other database's id  related to gene",
        "content": " \n".join(genes_content),
        "rdbVersion": "12.0"
    }
    return genes_doc
