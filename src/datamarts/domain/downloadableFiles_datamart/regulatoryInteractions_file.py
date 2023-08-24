import multigenomic_api
from datetime import datetime


class RegulatoryInteractions:

    @property
    def objects(self):
        ri_objects = multigenomic_api.regulatory_interactions.get_all()
        for ri_object in ri_objects:
            print(ri_object.id)
            ri_row = RegulatoryInteractions.RIDatamart(ri_object)
            yield ri_row
        del ri_objects

    class RIDatamart:
        def __init__(self, ri):
            self.ri = ri
            self.type = ri
            self.transcription_factor = ri.regulator
            self.regulatory_site = ri.regulatory_sites_id
            self.strand = ri.regulated_entity
            self.promoter = ri.regulated_entity
            self.tss = self.promoter
            self.target_tu_gene = ri.regulated_entity
            self.ri_evidences = ri.citations
            self.reg_site_evidences = self.regulatory_site
            self.additive_evidences = ri.additive_evidences_ids
            self.ri_ev_tech = ri.citations
            self.ri_ev_category = ri.citations

        @property
        def type(self):
            return self._type

        @type.setter
        def type(self, ri):
            self._type = None
            tf = multigenomic_api.transcription_factors.find_tf_id_by_conformation_id(ri.regulator.id)
            if len(tf) == 0:
                if ri.regulator.type == "regulatoryComplex":
                    reg_complex = multigenomic_api.regulatory_complexes.find_by_id(ri.regulator.id)
                    tf = multigenomic_api.transcription_factors.find_tf_id_by_conformation_name(reg_complex.abbreviated_name)

            if len(tf) > 0:
                regulator_type = "tf"
            elif ri.regulator.type == "product":
                regulator_type = "srna"
            else:
                regulator_type = "compound"

            if ri.regulated_entity.type == "gene":
                self._type = f"{regulator_type}-gene"
            elif ri.regulated_entity.type == "promoter":
                self._type = f"{regulator_type}-promoter"
            if ri.regulated_entity.type == "transcriptionUnit":
                self._type = f"{regulator_type}-tu"

        @property
        def transcription_factor(self):
            return self._transcription_factor

        @transcription_factor.setter
        def transcription_factor(self, regulator):
            self._transcription_factor = {
                "id": None,
                "name": None
            }
            tf = multigenomic_api.transcription_factors.find_tf_id_by_conformation_id(regulator.id)
            if len(tf) == 0:
                if regulator.type == "regulatoryComplex":
                    reg_complex = multigenomic_api.regulatory_complexes.find_by_id(regulator.id)
                    tf = multigenomic_api.transcription_factors.find_tf_id_by_conformation_name(reg_complex.abbreviated_name)
            if len(tf) > 0:
                self._transcription_factor = tf[0]

        @property
        def regulatory_site(self):
            return self._regulatory_site

        @regulatory_site.setter
        def regulatory_site(self, reg_site_id):
            self._regulatory_site = {
                "id": None,
                "left_end_position": None,
                "right_end_position": None,
                "strand": None,
                "sequence": None,
                "citations": []
            }
            if reg_site_id:
                self._regulatory_site = multigenomic_api.regulatory_sites.find_by_id(reg_site_id)

        @property
        def strand(self):
            return self._strand

        @strand.setter
        def strand(self, regulated_entity):
            self._strand = None
            if regulated_entity.type == "promoter":
                promoter = multigenomic_api.promoters.find_by_id(regulated_entity.id)
                self._strand = promoter.strand
            elif regulated_entity.type == "gene":
                gene = multigenomic_api.genes.find_by_id(regulated_entity.id)
                self._strand = gene.strand
            elif regulated_entity.type == "transcriptionUnit":
                trans_unit = multigenomic_api.transcription_units.find_by_id(regulated_entity.id)
                gene = multigenomic_api.genes.find_by_id(trans_unit.genes_ids[0])
                self._strand = gene.strand

        @property
        def promoter(self):
            return self._promoter

        @promoter.setter
        def promoter(self, reg_entity):
            self._promoter = {
                "id": None,
                "name": None
            }
            if reg_entity.type == "promoter":
                self._promoter = multigenomic_api.promoters.find_by_id(reg_entity.id)

        @property
        def tss(self):
            return self._tss

        @tss.setter
        def tss(self, promoter):
            self._tss = None
            try:
                if promoter.transcription_start_site:
                    self._tss = promoter.transcription_start_site.left_end_position
            except AttributeError:
                print("No promoter or tss")

        @property
        def target_tu_gene(self):
            return self._target_tu_gene

        @target_tu_gene.setter
        def target_tu_gene(self, reg_entity):
            self._target_tu_gene = None
            if reg_entity.type == "gene":
                self._target_tu_gene = f"{reg_entity.id}:{reg_entity.name}"
            elif reg_entity.type == "transcriptionUnit":
                self._target_tu_gene = f"{reg_entity.id}:{reg_entity.name}"
            if reg_entity.type == "promoter":
                tu = multigenomic_api.transcription_units.find_by_promoter_id(reg_entity.id)
                if tu:
                    self._target_tu_gene = f"{tu[0].id}:{tu[0].name}"

        @property
        def reg_site_evidences(self):
            return self._reg_site_evidences

        @reg_site_evidences.setter
        def reg_site_evidences(self, reg_site):
            self._reg_site_evidences = None
            if len(reg_site["citations"]) > 0:
                self._reg_site_evidences = ""
                for citation in reg_site.citations:
                    if citation.evidences_id:
                        citation_dict = multigenomic_api.evidences.find_by_id(citation.evidences_id)
                        self._reg_site_evidences += f"[{citation_dict.code}:{citation_dict.type}]"

        @property
        def ri_evidences(self):
            return self._ri_evidences

        @ri_evidences.setter
        def ri_evidences(self, citations):
            self._ri_evidences = None
            if len(citations) > 0:
                self._ri_evidences = ""
                for citation in citations:
                    if citation.evidences_id:
                        citation_dict = multigenomic_api.evidences.find_by_id(citation.evidences_id)
                        self._ri_evidences += f"[{citation_dict.code}:{citation_dict.type}]"
                if len(self._ri_evidences) == 0:
                    self._evidences = None

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

        @property
        def ri_ev_tech(self):
            return self._ri_ev_tech

        @ri_ev_tech.setter
        def ri_ev_tech(self, citations):
            self._ri_ev_tech = None
            if len(citations) > 0:
                self._ri_ev_tech = ""
                for citation in citations:
                    if citation.evidences_id:
                        citation_dict = multigenomic_api.evidences.find_by_id(citation.evidences_id)
                        self._ri_ev_tech += f"{citation_dict.approach}|"
                if len(self._ri_ev_tech) == 0:
                    self._ri_ev_tech = None
                else:
                    self._ri_ev_tech = self._ri_ev_tech[:-1]

        @property
        def ri_ev_category(self):
            return self._ri_ev_category

        @ri_ev_category.setter
        def ri_ev_category(self, citations):
            self._ri_ev_category = None
            if len(citations) > 0:
                self._ri_ev_category = ""
                for citation in citations:
                    if citation.evidences_id:
                        citation_dict = multigenomic_api.evidences.find_by_id(citation.evidences_id)
                        self._ri_ev_category += f"{citation_dict.category}|"
                if len(self._ri_ev_category) == 0:
                    self._ri_ev_category = None
                else:
                    self._ri_ev_category = self._ri_ev_category[:-1]

        def to_row(self):
            regulated_genes = get_regulated_genes(self.ri.regulated_entity)
            first_gene = get_first_gene(self.ri, regulated_genes)
            # TO CHECK: Strand on reg_site
            if self.regulatory_site:
                return f"{self.ri.id}" \
                       f"\t{self.type}" \
                       f"\t{self.transcription_factor['id']}" \
                       f"\t{self._transcription_factor['name']}" \
                       f"\t{self.ri.regulator.name}" \
                       f"\t{self.regulatory_site['id']}" \
                       f"\t{self.regulatory_site['left_end_position']}" \
                       f"\t{self.regulatory_site['right_end_position']}" \
                       f"\t{self.strand}" \
                       f"\t{self.regulatory_site['sequence']}" \
                       f"\t{self.ri.function}" \
                       f"\t{self.promoter['id']}" \
                       f"\t{self.promoter['name']}" \
                       f"\t{self.tss}" \
                       f"\t{self.ri.dist_site_promoter}" \
                       f"\t{first_gene['name']}" \
                       f"\t{first_gene['distanceTo']}" \
                       f"\t{self.target_tu_gene}" \
                       f"\t{self.ri.confidence_level}" \
                       f"\t{self.reg_site_evidences}" \
                       f"\t{self.ri_evidences}" \
                       f"\t{self.additive_evidences}" \
                       f"\t{self.ri_ev_tech}" \
                       f"\t{self.ri_ev_category}" \



def all_ris_rows():
    ris = RegulatoryInteractions()
    ris_content = []
    ris_content.append("1)riId	2)riType	3)tfId	4)tfName	5)cnfName	6)tfrsID	7)tfrsLeft	8)tfrsRight	9)strand	10)tfrsSeq	11)riFunction	12)promoterID	13)promoterName	14)tss	15)tfrsDistToPm	16)firstGene	17)tfrsDistTo1Gene	18)targetTuOrGene	19)confidenceLevel	20)tfrsEvidence	21)riEvidence	22)addEvidence	23)riEvTech	24)riEvCategory")
    for ri in ris.objects:
        ris_content.append(ri.to_row())
    creation_date = datetime.now()
    ri_doc = {
        "_id": "RDBECOLIDLF00001",
        "fileName": "RISet",
        "title": "Complete RIs Set",
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
        "columnsDetails": "# Columns:\t# (1) riId. Regulatory interaction (RI) identifier assigned by RegulonDB\t# (2) riType. Regulatory interaction type [tf-promoter, tf-tu, tf-gene]\t# (3) tfId. Transcription Factor (TF) identifier assigned by RegulonDB\t# (4) tfName. TF name\t# (5) cnfName. TF active conformation name\t# (6) tfrsId. TF regulatory site (TFRS) identifier assigned by RegulonDB\t# (7) tfrsLeft. TFRS left end position in the genome\t# (8) tfrsRight. TFRS right end position in the genome\t# (9) strand. DNA strand where the TFRS is located\t# (10) tfrsSeq. TFRS sequence (upper case)\t# (11) riFunction. Gene expression effect caused by the TF bound to the TFRS\t# (12) promoterId. Promoter Identifier assigned by RegulonDB\t# (13) promoterName. Promoter name\t# (14) tss. Transcription start site (+1) position in the genome\t# (15) tfrsDistToPm. Relative distance from the center position of TFRS to the Transcription Start Site\t# (16) firstGene. first transcribed gene name\t# (17) tfrsDistTo1Gene. Relative distance from center position of TFRS to the start of first gene\t# (18) targetTuOrGene. Transcription unit or gene (id:name) regulated by the TF\t# (19) confidenceLevel. RI confidence level (Values: Confirmed, Strong, Weak)\t# (20) tfrsEvidence. Evidence that supports the existence of the TFRS [EvidenceCode|EvidenceType(C:confirmed S:strong W:weak)]]\t# (21) riEvidence. Evidence that supports the RI function [EvidenceCode|EvidenceType(C:confirmed S:strong W:weak)]\t# (22) addEvidence. Additive Evidence [CV(EvidenceCode1/EvidenceCodeN)|Confidence Level]\t# (23) riEvTech. Evidence related to the type of technology used to determine the RI\t# (24) riEvCategory. Evidence  were categorized in classical, ht or non-experimental",
        "content": " \n".join(ris_content)
    }
    return ri_doc


def get_first_gene(reg_int, regulated_genes):
    distance_to_first_gene = None
    promoter = None
    first_gene = None
    if len(regulated_genes) > 0:
        if reg_int.regulatory_sites_id:
            reg_sites = multigenomic_api.regulatory_sites.find_by_id(reg_int.regulatory_sites_id)
            if reg_int.regulated_entity.type == "gene":
                first_gene = multigenomic_api.genes.find_by_id(reg_int.regulated_entity.id)
                if first_gene.strand:
                    if reg_sites.absolute_position:
                        if first_gene.strand == "forward":
                            distance_to_first_gene = reg_sites.absolute_position - first_gene.left_end_position
                        else:
                            distance_to_first_gene = first_gene.right_end_position - reg_sites.absolute_position
                        return {
                            "name": first_gene["name"],
                            "distanceTo": distance_to_first_gene
                        }
            else:
                if reg_int.regulated_entity.type == "promoter":
                    promoter = multigenomic_api.promoters.find_by_id(reg_int.regulated_entity.id)
                    first_gene = get_first_gene_of_tu(regulated_genes, promoter)
                elif reg_int.regulated_entity.type == "transcriptionUnit":
                    trans_unit = multigenomic_api.transcription_units.find_by_id(reg_int.regulated_entity.id)
                    if trans_unit.promoters_id:
                        promoter = multigenomic_api.promoters.find_by_id(trans_unit.promoters_id)
                    first_gene = get_first_gene_of_tu(regulated_genes, promoter)
                if promoter:
                    if reg_sites.absolute_position:
                        if promoter.strand == "forward":
                            distance_to_first_gene = reg_sites.absolute_position - first_gene["leftEndPosition"]
                        else:
                            distance_to_first_gene = first_gene["leftEndPosition"] - reg_sites.absolute_position
            return {
                "name": first_gene["name"],
                "distanceTo": distance_to_first_gene
            }
    return {
        "name": None,
        "distanceTo": None
    }


def get_first_gene_of_tu(genes, promoter):
    first_gene = None
    if len(genes) > 0:
        gene = genes[0]
        first_gene = multigenomic_api.genes.find_by_id(gene.get("_id"))
        if promoter:
            if promoter.strand == "reverse":
                first_gene.left_end_position = first_gene.right_end_position
            first_gene.left_end_position = first_gene.left_end_position or first_gene.fragments[0].left_end_position

            for gene in genes:
                current_gene = multigenomic_api.genes.find_by_id(gene.get("_id"))
                if promoter.strand == "forward":
                    if current_gene.left_end_position:
                        if current_gene.left_end_position < first_gene.left_end_position:
                            first_gene = current_gene
                    elif current_gene.fragments:
                        for fragment in current_gene.fragments:
                            if fragment.left_end_position < first_gene.left_end_position:
                                first_gene = current_gene
                                first_gene.left_end_position = fragment.left_end_position

                elif promoter.strand == "reverse":
                    if current_gene.left_end_position:
                        if current_gene.right_end_position > first_gene.left_end_position:
                            first_gene = current_gene
                            first_gene.left_end_position = first_gene.right_end_position
                    elif current_gene.fragments:
                        for fragment in current_gene.fragments:
                            if fragment.right_end_position > first_gene.left_end_position:
                                first_gene = current_gene
                                first_gene.left_end_position = fragment.right_end_position

                first_gene.left_end_position = first_gene.left_end_position or first_gene.fragments[0].left_end_position
        first_gene = {
            "id": first_gene.id,
            "name": first_gene.name,
            "leftEndPosition": first_gene.left_end_position
        }
    return first_gene


def get_regulated_genes(regulated_entity):
    regulated_genes = []
    transcription_units = []
    if regulated_entity.type == "promoter":
        transcription_units = multigenomic_api.transcription_units.find_by_promoter_id(regulated_entity.id)
    elif regulated_entity.type == "transcriptionUnit":
        trans_unit = multigenomic_api.transcription_units.find_by_id(regulated_entity.id)
        transcription_units.append(trans_unit)
    if transcription_units:
        for tu in transcription_units:
            for gene_id in tu.genes_ids:
                gene = multigenomic_api.genes.find_by_id(gene_id)
                gene_object = {
                    "_id": gene.id,
                    "name": gene.name,
                }
                if gene_object not in regulated_genes:
                    regulated_genes.append(gene_object)
    return regulated_genes
