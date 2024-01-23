import multigenomic_api
from datetime import datetime
import re


class GeneProductIds:

    @property
    def objects(self):
        gene_objects = multigenomic_api.genes.get_all()
        for gene_object in gene_objects:
            # print(gene_object.id)
            ri_row = GeneProductIds.GeneProductIdsDatamart(gene_object)
            yield ri_row
        del gene_objects

    class GeneProductIdsDatamart:
        def __init__(self, gene):
            self.gene = gene
            self.gene_lep = gene
            self.gene_rep = gene
            self.products = gene.id
            self.otherDBsGeneIds = gene.external_cross_references

        @property
        def products(self):
            return self._products

        @products.setter
        def products(self, gene_id):
            self._products = multigenomic_api.products.find_by_gene_id(gene_id)

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
        def otherDBsGeneIds(self):
            return self._otherDBsGeneIds

        @otherDBsGeneIds.setter
        def otherDBsGeneIds(self, ext_cross_ref):
            self._otherDBsGeneIds = []
            for ext_ref in ext_cross_ref:
                if ext_ref.external_cross_references_id:
                    ext_ref_dict = multigenomic_api.external_cross_references.find_by_id(
                        ext_ref.external_cross_references_id)
                    ext_ref_item = f"[{ext_ref_dict.name}:{ext_ref.object_id}]"
                    if ext_ref_item not in self._otherDBsGeneIds:
                        self._otherDBsGeneIds.append(ext_ref_item)
            self._otherDBsGeneIds = "".join(self._otherDBsGeneIds)

        def to_row(self):
            products = []
            if len(self.products) > 0:
                for product in self.products:
                    products.append(f"{self.gene.id}"
                                    f"\t{self.gene.name}"
                                    f"\t{self.gene_lep}"
                                    f"\t{self.gene_rep}"
                                    f"\t{self.gene.strand}"
                                    f"\t{','.join(self.gene.synonyms)}"
                                    f"\t{self.otherDBsGeneIds}"
                                    f"\t{product.id}"
                                    f"\t{product.name}"
                                    f"\t{','.join(product.synonyms)}"
                                    f"\t{get_other_prod_ids(product.external_cross_references)}")
                return "\n".join(products)
            else:
                return f"{self.gene.id}" \
                       f"\t{self.gene.name}" \
                       f"\t{self.gene_lep}" \
                       f"\t{self.gene_rep}" \
                       f"\t{self.gene.strand}" \
                       f"\t{','.join(self.gene.synonyms)}" \
                       f"\t{self.otherDBsGeneIds}" \
                       f"\t{''}" \
                       f"\t{''}" \
                       f"\t{''}" \
                       f"\t{''}"


def all_gene_rows():
    genes = GeneProductIds()
    genes_content = [
        "1)geneId\t2)geneName\t3)leftEndPos\t4)rightEndPos\t5)strand\t6)geneSynonyms\t7)otherDbsGeneIds\t8)productId\t9)productName\t10)productSynonyms\t11)otherDbsProductsIds"]
    for gene in genes.objects:
        genes_content.append(gene.to_row())
    creation_date = datetime.now()
    genes_doc = {
        "_id": "RDBECOLIDLF00017",
        "fileName": "GeneProductAllIdentifiersSet",
        "title": "Complete Gene Product Identifiers Set",
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
        "columnsDetails": "# Columns:\n# (1) Gene identifier assigned by RegulonDB\n# (2) Gene name\n# (3) Gene left end position in the genome\n# (4) Gene right end position in the genome\n# (5) DNA strand where the gene is coded\n# (6) other gene synonyms\n# (7) Other database's id  related to gene\n# (8) Product identifier of the gene\n# (9) Product name of the gene\n# (10) Other products synonyms\n# (11) Other database's id  related to product",
        "content": " \n".join(genes_content),
        "rdbVersion": "12.0",
        "description": "Genes with information about their products, including synonyms and identifiers to other databases. Useful for database mapping and cross-referencing purposes.",
        "group": "GENE"
    }
    return genes_doc


def get_other_prod_ids(ext_cross_ref):
    otherDBsGeneIds = []
    for ext_ref in ext_cross_ref:
        if ext_ref.external_cross_references_id:
            ext_ref_dict = multigenomic_api.external_cross_references.find_by_id(
                ext_ref.external_cross_references_id)
            ext_ref_item = f"[{ext_ref_dict.name}:{ext_ref.object_id}]"
            if ext_ref_item not in otherDBsGeneIds:
                otherDBsGeneIds.append(ext_ref_item)
    return "".join(otherDBsGeneIds)
