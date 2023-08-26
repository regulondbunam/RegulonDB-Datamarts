import multigenomic_api
from datetime import datetime


class GeneProduct:

    @property
    def objects(self):
        gene_objects = multigenomic_api.genes.get_all()
        for gene_object in gene_objects:
            print(gene_object.id)
            ri_row = GeneProduct.GeneProductDatamart(gene_object)
            yield ri_row
        del gene_objects

    class GeneProductDatamart:
        def __init__(self, gene):
            self.gene = gene
            self.gene_lep = gene
            self.gene_rep = gene
            self.product = gene.id
            self.gene_pmids = gene.citations
            self.gene_evidences = gene.citations
            self.additive_evidences = gene.additive_evidences_ids

        @property
        def product(self):
            return self._product

        @product.setter
        def product(self, gene_id):
            self._product = None
            products = multigenomic_api.products.find_by_gene_id(gene_id)
            if len(products) > 0:
                self._product = ""
                for product in products:
                    self._product += f"{product.name},"
                self._product = self._product[:-1]

        @property
        def genes(self):
            return self._genes

        @genes.setter
        def genes(self, genes_ids):
            self._genes = None
            if len(genes_ids) > 0:
                self._genes = ""
                for gene_id in genes_ids:
                    gene = multigenomic_api.genes.find_by_id(gene_id)
                    self._genes += f"{gene.name},"
            if len(self._genes) > 0:
                self._genes = self._genes[:-1]

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
        def gene_pmids(self):
            return self._gene_pmids

        @gene_pmids.setter
        def gene_pmids(self, citations):
            self._gene_pmids = None
            if citations:
                self._gene_pmids = []
                for citation in citations:
                    if citation.publications_id:
                        publication = multigenomic_api.publications.find_by_id(citation.publications_id)
                        if publication.pmid not in self._gene_pmids:
                            self._gene_pmids.append(f"{publication.pmid},")
                self._gene_pmids = "".join(self._gene_pmids)[:-1]

        @property
        def gene_evidences(self):
            return self._gene_evidences

        @gene_evidences.setter
        def gene_evidences(self, citations):
            self._gene_evidences = None
            if len(citations) > 0:
                self._gene_evidences = ""
                for citation in citations:
                    if citation.evidences_id:
                        citation_dict = multigenomic_api.evidences.find_by_id(citation.evidences_id)
                        self._gene_evidences += f"[{citation_dict.code}:{citation_dict.type}]"
                if len(self._gene_evidences) == 0:
                    self._gene_evidences = None

        @property
        def additive_evidences(self):
            return self._additive_evidences

        @additive_evidences.setter
        def additive_evidences(self, additive_evs_ids):
            self._additive_evidences = None
            if len(additive_evs_ids) > 0:
                self._additive_evidences = ""
                for additive_evs_id in additive_evs_ids:
                    additive_evidence_dict = multigenomic_api.additive_evidences.find_by_id(additive_evs_id)
                    self._additive_evidences += f"[{additive_evidence_dict.code}:{additive_evidence_dict.confidence_level}]"

        def to_row(self):
            return f"{self.gene.id}" \
                   f"\t{self.gene.name}" \
                   f"\t{self.gene.bnumber}" \
                   f"\t{self.gene_lep}" \
                   f"\t{self.gene_rep}" \
                   f"\t{self.gene.strand}" \
                   f"\t{self.product}" \
                   f"\t{self.gene_evidences}" \
                   f"\t{self.gene_pmids}" \
                   f"\t{self.additive_evidences}" \
                   f"\t{self.gene.confidence_level}" \
                   f"\t{None}" \



def all_gene_rows():
    genes = GeneProduct()
    genes_content = ["1)geneId\t2)geneName\t3)bnumber\t4)leftEndPos\t5)rightEndPos\t6)strand\t7)productName\t8)geneEvidences\t9)genePMIDS\t10)additiveEvidences\t11)relatedPMIDS\t12)otherDbsIds"]
    for gene in genes.objects:
        genes_content.append(gene.to_row())
    creation_date = datetime.now()
    genes_doc = {
        "_id": "RDBECOLIDLF00006",
        "fileName": "GeneProductSet",
        "title": "Complete Gene Product Set",
        "fileFormat": "rif-version 1",
        "license": "RegulonDB is free for academic/noncommercial use\t\tUser is not entitled to change or erase data sets of the RegulonDB\tdatabase or to eliminate copyright notices from RegulonDB. Furthermore,\tUser is not entitled to expand RegulonDB or to integrate RegulonDB partly\tor as a whole into other databank systems, without prior written consent\tfrom CCG-UNAM.\t\tPlease check the license at http://regulondb.ccg.unam.mx/menu/download/full_version/terms_and_conditions.jsp",
        "citation": "Tierrafría, V. H. et al. (2022). RegulonDB 11.0: Comprehensive high-throughput datasets on transcriptional regulation in Escherichia coli K-12,\tMicrob Genom. 2022 May;8(5). doi: 10.1099/mgen.0.000833. PMID: 35584008. https://doi.org/10.1099/mgen.0.000833",
        "contact": {
            "person": "RegulonDB Team",
            "webPage": "http://regulondb.ccg.unam.mx/menu/about_regulondb/contact_us/index.jsp",
            "email": "regulondb@ccg.unam.mx"
        },
        "version": "",
        "creationDate": f"{creation_date.strftime('%m-%d-%Y')}",
        "columnsDetails": "Columns:\n(1) Gene identifier assigned by RegulonDB\n(2) Gene name\n(3) Blattner number (bnumber) of the gene\n(4) Gene left end position in the genome\n(5) Gene right end position in the genome\n(6) DNA strand where the gene is coded\n(7) Product name of the gene\n(8) Evidence that supports the existence of the gene\n(9) PMIDs list\n(10) Evidence confidence level (Confirmed, Strong, Weak)\n(11) All bnumber related to gene\n(12) Other database's id  related to gene",
        "content": " \n".join(genes_content)
    }
    return genes_doc
