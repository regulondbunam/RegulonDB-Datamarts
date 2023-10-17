import multigenomic_api
from datetime import datetime


class Promoters:

    @property
    def objects(self):
        promoter_objects = multigenomic_api.promoters.get_all()
        for promoter_object in promoter_objects:
            # print(promoter_object.id)
            ri_row = Promoters.PromoterDatamart(promoter_object)
            yield ri_row
        del promoter_objects

    class PromoterDatamart:
        def __init__(self, promoter):
            self.promoter = promoter
            self.tss = promoter
            self.first_gene = promoter
            self.prom_evidences = promoter.citations
            self.additive_evidences = promoter.additive_evidences_ids
            self.sigma_factor = promoter.binds_sigma_factor

        @property
        def tss(self):
            return self._tss

        @tss.setter
        def tss(self, promoter):
            self._tss = ""
            try:
                if promoter.transcription_start_site:
                    self._tss = promoter.transcription_start_site.left_end_position
            except AttributeError:
                pass

        @property
        def first_gene(self):
            return self._first_gene

        @first_gene.setter
        def first_gene(self, promoter):
            self._first_gene = {
                        "distanceToPromoter": "",
                        "name": ""
                    }
            dist = None
            trans_units = multigenomic_api.transcription_units.find_by_promoter_id(promoter.id)
            for tu in trans_units:
                first_gene = get_first_gene_of_tu(tu, promoter)
                if promoter.transcription_start_site:
                    if promoter.strand == "forward":
                        dist = promoter.transcription_start_site.left_end_position - first_gene["leftEndPosition"]
                    if promoter.strand == "reverse":
                        dist = first_gene["rightEndPosition"] - promoter.transcription_start_site.right_end_position
                    self._first_gene = {
                        "distanceToPromoter": dist,
                        "name": first_gene["name"]
                    }

        @property
        def sigma_factor(self):
            return self._sigma_factor

        @sigma_factor.setter
        def sigma_factor(self, sigma_factor):
            self._sigma_factor = ""
            if sigma_factor:
                sigma_factor = multigenomic_api.sigma_factors.find_by_id(sigma_factor.sigma_factors_id)
                self._sigma_factor = sigma_factor.abbreviated_name

        @property
        def prom_evidences(self):
            return self._prom_evidences

        @prom_evidences.setter
        def prom_evidences(self, citations):
            self._prom_evidences = []
            for citation in citations:
                if citation.evidences_id:
                    citation_dict = multigenomic_api.evidences.find_by_id(citation.evidences_id)
                    citation_item = f"[{citation_dict.code}:{citation_dict.type}]"
                    if citation_item not in self._prom_evidences:
                        self._prom_evidences.append(citation_item)
            self._prom_evidences = "".join(self._prom_evidences)

        @property
        def additive_evidences(self):
            return self._additive_evidences

        @additive_evidences.setter
        def additive_evidences(self, additive_evs_ids):
            self._additive_evidences = ""
            for additive_evs_id in additive_evs_ids:
                additive_evidence_dict = multigenomic_api.additive_evidences.find_by_id(additive_evs_id)
                self._additive_evidences += f"[{additive_evidence_dict.code}:{additive_evidence_dict.confidence_level}]"

        def to_row(self):
            return f"{self.promoter.id}" \
                   f"\t{self.promoter.name}" \
                   f"\t{self.promoter.strand}" \
                   f"\t{self.tss}" \
                   f"\t{self.sigma_factor}" \
                   f"\t{self.promoter.sequence}" \
                   f"\t{self.first_gene['name']}" \
                   f"\t{self.first_gene['distanceToPromoter']}" \
                   f"\t{self.prom_evidences}" \
                   f"\t{self.additive_evidences}" \
                   f"\t{self.promoter.confidence_level}"


def get_first_gene_of_tu(transcription_unit, promoter):
    dict_genes = []
    first_gene = {}
    for gene in transcription_unit.genes_ids:
        gene_object = multigenomic_api.genes.find_by_id(gene)
        if gene_object.fragments:
            min_left_pos = min(gene_object.fragments, key=lambda x: x.left_end_position)
            max_right_pos = max(gene_object.fragments, key=lambda x: x.right_end_position)
            gene_object.left_end_position = min_left_pos.left_end_position
            gene_object.right_end_position = max_right_pos.right_end_position
        dict_genes.append({
            "id": gene_object.id,
            "name": gene_object.name,
            "leftEndPosition": gene_object.left_end_position,
            "rightEndPosition": gene_object.right_end_position
        })
    if promoter.strand == "forward":
        first_gene = (min(dict_genes, key=lambda x: x["leftEndPosition"]))
    if promoter.strand == "reverse":
        first_gene = (max(dict_genes, key=lambda x: x["rightEndPosition"]))
    return first_gene



def all_promoters_rows():
    promoters = Promoters()
    promoters_content = ["1)pmId\t2)pmName\t3)strand\t4)posTSS\t5)sigmaF\t6)pmSequence\t7)firstGeneName\t8)distToFirstGene\t9)pmEvidence\t10)addEvidence\t11)confidenceLevel"]
    for promoter in promoters.objects:
        promoters_content.append(promoter.to_row())
    creation_date = datetime.now()
    promoters_doc = {
        "_id": "RDBECOLIDLF00011",
        "fileName": "PromoterSet",
        "title": "Complete Promoters Set",
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
        "columnsDetails": "Columns:\n(1) pmId. Promoter identifier assigned by RegulonDB\n(2) pmName. Promoter Name\n(3) strand. DNA strand where the promoter is located\n(4) posTSS. Genome map position of Transcription Start Site (+1)\n(5) sigmaF. Sigma Factor that recognize the promoter\n(6) pmSequence. Promoter Sequence (+1 upper case)\n(7) firstGeneName. Name of the first gene of promoter\n(8) distToFirstGene. distance to first gene of promoter\n(9) pmEvidence. Evidence that supports the existence of the promoter\n(10)addEvidence. Additive Evidence [CV(EvidenceCode1/EvidenceCodeN)|Confidence Level]\n(11) confidenceLevel. Promoter confidence level (Values: Confirmed, Strong, Weak)",
        "content": " \n".join(promoters_content),
        "rdbVersion": "12.0"
    }
    return promoters_doc
