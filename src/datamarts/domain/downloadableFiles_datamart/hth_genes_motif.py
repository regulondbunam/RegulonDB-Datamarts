import multigenomic_api
from datetime import datetime
import re


class GeneProduct:

    @property
    def objects(self):
        gene_objects = multigenomic_api.genes.get_all()
        for gene_object in gene_objects:
            # print(gene_object.id)
            ri_row = GeneProduct.GeneProductDatamart(gene_object)
            yield ri_row
        del gene_objects

    class GeneProductDatamart:
        def __init__(self, gene):
            products = multigenomic_api.products.find_by_gene_id(gene.id)
            self.gene = gene
            self.gene_lep = gene
            self.gene_rep = gene
            self.products = gene.id
            self.rows = products

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
        def rows(self):
            return self._rows

        @rows.setter
        def rows(self, products):
            self._rows = []
            pattern = r"UniProt: H[-]?T[-]?H"
            regulators_name = []
            for product in products:
                product_size = len(product.sequence) if getattr(product, "sequence", None) else ""
                regulators = multigenomic_api.regulators.find_by_conformation(product.id)
                for regulator in regulators:
                    regulators_name.append(regulator.abbreviated_name)
                motifs = multigenomic_api.motifs.find_by_product_id(product.id)
                for motif in motifs:
                    if motif.note:
                        if re.search(pattern, motif.note):
                            self._rows.append(f"{self.gene.id}" 
                                               f"\t{self.gene.name}" 
                                               f"\t{self.gene_lep}" 
                                               f"\t{self.gene_rep}" 
                                               f"\t{self.gene.bnumber}" 
                                               f"\t{product.abbreviated_name}" 
                                               f"\t{",".join(regulators_name)}" 
                                               f"\t{product_size}" 
                                               f"\t{motif.id}" 
                                               f"\t{motif.left_end_position}" 
                                               f"\t{motif.right_end_position}" 
                                               f"\t{motif.sequence}")
        def to_row(self):
            return self.rows


def all_gene_rows(rdb_version, citation):
    genes = GeneProduct()
    genes_content = [
        "1)geneId\t2)geneName\t3)geneBnumber\t4)leftEndPos\t5)rightEndPos\t6)productName\t7)regulatorsName\t8)productSize\t9)motifId\t10)motif_leftEndPos\t11)motif_rightEndPos\t12)sequence"]
    for gene in genes.objects:
        if gene.to_row():
            genes_content.extend(gene.to_row())
    creation_date = datetime.now()
    genes_doc = {
        "_id": "RDBECOLIDLF00023",
        "fileName": "hth-genes-motif_internal",
        "title": "Complete Gene HTH Motif Set",
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
                          "# (1) Gene identifier assigned by RegulonDB\n"
                          "# (2) Gene name\n"
                          "# (3) Bnumber of the gene\n"
                          "# (4) Gene left end position in the genome\n"
                          "# (5) Gene right end position in the genome\n"
                          "# (6) Product Name\n"
                          "# (7) Regulator(s) Name\n"
                          "# (8) Product length\n"
                          "# (9) Id of the HTH motif\n"
                          "# (10) Motif left end position\n"
                          "# (11) Motif right end position\n"
                          "# (12) Sequence of the motif",
        "content": " \n".join(genes_content),
        "rdbVersion": rdb_version,
        "description": "Collection of genes with HTH motifs and its regulator that encodes.",
        "group": "GENE"
    }
    return genes_doc
