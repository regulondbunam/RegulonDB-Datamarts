import multigenomic_api
from datetime import datetime


class RegulatoryInteractions:

    @property
    def objects(self):
        ri_objects = multigenomic_api.regulatory_interactions.get_all()
        for ri_object in ri_objects:
            # print(ri_object.id)
            ri_row = RegulatoryInteractions.RIDatamart(ri_object)
            yield ri_row
        del ri_objects

    class RIDatamart:
        def __init__(self, ri):
            self.ri = ri
            self.type = ri
            self.transcription_factor = ri.regulator
            self.regulatory_site = ri.regulatory_sites_id
            self.regulator_name = ri.regulator
            self.strand = ri.regulated_entity
            self.sequence = self.regulatory_site
            self.promoter = ri.regulated_entity
            self.tss = self.promoter
            self.sigma = self.promoter
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
            self._type = ""
            tf = multigenomic_api.transcription_factors.find_tf_id_by_conformation_id(ri.regulator.id)
            if tf is None:
                tf = multigenomic_api.transcription_factors.find_by_name(ri.regulator.name)
            if tf:
                regulator_type = "tf"
            elif ri.regulator.type == "product":
                regulator_type = "srna"
            elif ri.regulator.type == "regulatoryContinuant":
                regulator_type = "compound"
            else:
                regulator_type = "regulator"

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
                "id": "",
                "abbreviated_name": ""
            }
            tf = multigenomic_api.transcription_factors.find_tf_id_by_conformation_id(regulator.id)
            if tf is None:
                tf = multigenomic_api.transcription_factors.find_by_name(regulator.name)
                if tf is not None:
                    self._transcription_factor = tf
                else:
                    if regulator.type == "product":
                        product = multigenomic_api.products.find_by_id(regulator.id)
                        abb_name = product.abbreviated_name
                    else:
                        abb_name = regulator.name
                    self._transcription_factor = {
                        "id": regulator.id,
                        "abbreviated_name": abb_name
                    }
            else:
                self._transcription_factor = {
                    "id": tf.id,
                    "abbreviated_name": tf.abbreviated_name or tf.name
                }

        @property
        def regulator_name(self):
            return self._regulator_name

        @regulator_name.setter
        def regulator_name(self, regulator):
            self._regulator_name = ""
            if regulator.type == "product":
                reg = multigenomic_api.products.find_by_id(regulator.id)
                self._regulator_name = reg.abbreviated_name or reg.name
            if regulator.type == "regulatoryComplex":
                reg = multigenomic_api.regulatory_complexes.find_by_id(regulator.id)
                if reg.abbreviated_name or reg.name:
                    self._regulator_name = reg.abbreviated_name or reg.name
                else:
                    self._regulator_name = regulator.name
            if regulator.type == "regulatoryContinuant":
                reg = multigenomic_api.regulatory_continuants.find_by_id(regulator.id)
                self._regulator_name = reg.name

        @property
        def regulatory_site(self):
            return self._regulatory_site

        @regulatory_site.setter
        def regulatory_site(self, reg_site_id):
            self._regulatory_site = {
                "id": "",
                "left_end_position": "",
                "right_end_position": "",
                "strand": "",
                "sequence": "",
                "citations": []
            }
            if reg_site_id:
                self._regulatory_site = multigenomic_api.regulatory_sites.find_by_id(reg_site_id)

        @property
        def strand(self):
            return self._strand

        @strand.setter
        def strand(self, regulated_entity):
            self._strand = ""
            if regulated_entity.type == "promoter":
                promoter = multigenomic_api.promoters.find_by_id(regulated_entity.id)
                self._strand = promoter.strand
            elif regulated_entity.type == "gene":
                gene = multigenomic_api.genes.find_by_id(regulated_entity.id)
                self._strand = gene.strand
            elif regulated_entity.type == "transcriptionUnit":
                trans_unit = multigenomic_api.transcription_units.find_by_id(regulated_entity.id)
                operon = multigenomic_api.operons.find_by_id(trans_unit.operons_id)
                self._strand = operon.strand

        @property
        def sequence(self):
            return self._sequence

        @sequence.setter
        def sequence(self, reg_site):
            self._sequence = ""
            if reg_site['sequence']:
                if self.strand == "reverse":
                    self._sequence = reverse_complement(reg_site['sequence'])
                else:
                    self._sequence = reg_site['sequence']

        @property
        def promoter(self):
            return self._promoter

        @promoter.setter
        def promoter(self, reg_entity):
            self._promoter = {
                "id": "",
                "name": ""
            }
            if reg_entity.type == "promoter":
                self._promoter = multigenomic_api.promoters.find_by_id(reg_entity.id)

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
        def sigma(self):
            return self._sigma

        @sigma.setter
        def sigma(self, promoter):
            self._sigma = ""
            try:
                if promoter.binds_sigma_factor.sigma_factors_id:
                    sigma = multigenomic_api.sigma_factors.find_by_id(promoter.binds_sigma_factor.sigma_factors_id)
                    self._sigma = sigma.abbreviated_name
            except AttributeError:
                pass

        @property
        def target_tu_gene(self):
            return self._target_tu_gene

        @target_tu_gene.setter
        def target_tu_gene(self, reg_entity):
            self._target_tu_gene = ""
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
            self._reg_site_evidences = []
            if len(reg_site["citations"]) > 0:
                for citation in reg_site.citations:
                    if citation.evidences_id:
                        if citation.evidences_id:
                            citation_dict = multigenomic_api.evidences.find_by_id(citation.evidences_id)
                            citation_item = f"[{citation_dict.code}:{citation_dict.type}]"
                            if citation_item not in self._reg_site_evidences:
                                self._reg_site_evidences.append(citation_item)
            self._reg_site_evidences = "".join(self._reg_site_evidences)

        @property
        def ri_evidences(self):
            return self._ri_evidences

        @ri_evidences.setter
        def ri_evidences(self, citations):
            self._ri_evidences = []
            for citation in citations:
                if citation.evidences_id:
                    citation_dict = multigenomic_api.evidences.find_by_id(citation.evidences_id)
                    citation_item = f"[{citation_dict.code}:{citation_dict.type}]"
                    if citation_item not in self._ri_evidences:
                        self._ri_evidences.append(citation_item)
            self._ri_evidences = "".join(self._ri_evidences)

        @property
        def additive_evidences(self):
            return self._additive_evidences

        @additive_evidences.setter
        def additive_evidences(self, additive_evs_ids):
            self._additive_evidences = ""
            for additive_evs_id in additive_evs_ids:
                additive_evidence_dict = multigenomic_api.additive_evidences.find_by_id(additive_evs_id)
                self._additive_evidences += f"[{additive_evidence_dict.code}:{additive_evidence_dict.confidence_level}]"

        @property
        def ri_ev_tech(self):
            return self._ri_ev_tech

        @ri_ev_tech.setter
        def ri_ev_tech(self, citations):
            self._ri_ev_tech = ""
            for citation in citations:
                if citation.evidences_id:
                    citation_dict = multigenomic_api.evidences.find_by_id(citation.evidences_id)
                    self._ri_ev_tech += f"{citation_dict.approach}|"
            self._ri_ev_tech = self._ri_ev_tech[:-1]

        @property
        def ri_ev_category(self):
            return self._ri_ev_category

        @ri_ev_category.setter
        def ri_ev_category(self, citations):
            self._ri_ev_category = ""
            ev_categories = []
            for citation in citations:
                if citation.evidences_id:
                    citation_dict = multigenomic_api.evidences.find_by_id(citation.evidences_id)
                    if citation_dict.category:
                        if citation_dict.category not in ev_categories:
                            ev_categories.append(citation_dict.category)
            self._ri_ev_category = "|".join(ev_categories)

        def to_row(self):
            regulated_genes = get_regulated_genes(self.ri.regulated_entity)
            first_gene = get_first_gene(self.ri, regulated_genes)
            return f"{self.ri.id}" \
                   f"\t{self.type}" \
                   f"\t{self.transcription_factor['id']}" \
                   f"\t{self.transcription_factor['abbreviated_name']}" \
                   f"\t{self.regulator_name}" \
                   f"\t{self.regulatory_site['id']}" \
                   f"\t{self.regulatory_site['left_end_position'] or ''}" \
                   f"\t{self.regulatory_site['right_end_position'] or ''}" \
                   f"\t{self.strand}" \
                   f"\t{self.sequence}" \
                   f"\t{self.ri.function or ''}" \
                   f"\t{self.promoter['id']}" \
                   f"\t{self.promoter['name']}" \
                   f"\t{self.tss}" \
                   f"\t{self.sigma}" \
                   f"\t{self.ri.dist_site_promoter or ''}" \
                   f"\t{first_gene['name']}" \
                   f"\t{first_gene['distanceTo']}" \
                   f"\t{self.target_tu_gene}" \
                   f"\t{self.ri.confidence_level or '?'}" \
                   f"\t{self.reg_site_evidences}" \
                   f"\t{self.ri_evidences}" \
                   f"\t{self.additive_evidences}" \
                   f"\t{self.ri_ev_tech}" \
                   f"\t{self.ri_ev_category}"


def all_ris_rows():
    ris = RegulatoryInteractions()
    ris_content = []
    ris_content.append("1)riId\t2)riType\t3)regulatorId\t4)regulatorName\t5)cnfName\t6)tfrsID\t7)tfrsLeft\t8)tfrsRight\t9)strand\t10)tfrsSeq\t11)riFunction\t12)promoterID\t13)promoterName\t14)tss\t15)sigmaF\t16)tfrsDistToPm\t17)firstGene\t18)tfrsDistTo1Gene\t19)targetTuOrGene\t20)confidenceLevel\t21)tfrsEvidence\t22)riEvidence\t23)addEvidence\t24)riEvTech\t25)riEvCategory")
    for ri in ris.objects:
        if ri.type == "tf-gene" or ri.type == "tf-promoter" or ri.type == "tf-tu":
            ris_content.append(ri.to_row())
    creation_date = datetime.now()
    ri_doc = {
        "_id": "RDBECOLIDLF00002",
        "fileName": "TF-RISet",
        "title": "Complete TF RIs Set",
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
        "columnsDetails": "# Columns:\n# (1) riId. Regulatory interaction (RI) identifier assigned by RegulonDB\n# (2) riType. Regulatory interaction type [tf-promoter, tf-tu, tf-gene,srna-promoter, srna-tu, srna-gene, compound-promoter, compound-tu, compound-gene]\n# (3) tfId. Transcription Factor (TF) identifier assigned by RegulonDB\n# (4) regulatorName. Regulator name\n# (5) cnfName. regulator active conformation name\n# (6) tfrsId. regulator regulatory site (TFRS) identifier assigned by RegulonDB\n# (7) tfrsLeft. TFRS left end position in the genome\n# (8) tfrsRight. TFRS right end position in the genome\n# (9) strand. DNA strand where the TFRS is located\n# (10) tfrsSeq. TFRS sequence (upper case)\n# (11) riFunction. Gene expression effect caused by the TF bound to the TFRS\n# (12) promoterId. Promoter Identifier assigned by RegulonDB\n# (13) promoterName. Promoter name\n# (14) tss. Transcription start site (+1) position in the genome\n# (15) sigmaF. Sigma Factor that recognize the promotertfrs\n# (16)DistToPm. Relative distance from the center position of TFRS to the Transcription Start Site\n# (17) firstGene. first transcribed gene name\n# (18) tfrsDistTo1Gene. Relative distance from center position of TFRS to the start of first gene\n# (19) targetTuOrGene. Transcription unit or gene (id:name) regulated by the TF\n# (20) confidenceLevel. RI confidence level (Values: Confirmed, Strong, Weak)\n# (21) tfrsEvidence. Evidence that supports the existence of the TFRS [EvidenceCode|EvidenceType(C:confirmed S:strong W:weak)]]\n# (22) riEvidence. Evidence that supports the RI function [EvidenceCode|EvidenceType(C:confirmed S:strong W:weak)]\n# (23) addEvidence. Additive Evidence [CV(EvidenceCode1/EvidenceCodeN)|Confidence Level]\n# (24) riEvTech. Evidence related to the type of technology used to determine the RI\n# (25) riEvCategory. Evidence  were categorized in classical, ht or non-experimental",
        "content": " \n".join(ris_content),
        "rdbVersion": "12.0"
    }
    return ri_doc


def get_first_gene(reg_int, regulated_genes):
    strand = None
    distance_to_first_gene = ""
    first_gene = None
    if len(regulated_genes) > 0:
        if reg_int.regulated_entity.type == "gene":
            first_gene = multigenomic_api.genes.find_by_id(reg_int.regulated_entity.id)
            if first_gene.strand:
                strand = first_gene.strand
        else:
            if reg_int.regulated_entity.type == "promoter":
                promoter = multigenomic_api.promoters.find_by_id(reg_int.regulated_entity.id)
                strand = promoter.strand
                first_gene = get_first_gene_of_tu(regulated_genes, strand)
            elif reg_int.regulated_entity.type == "transcriptionUnit":
                trans_unit = multigenomic_api.transcription_units.find_by_id(reg_int.regulated_entity.id)
                operon = multigenomic_api.operons.find_by_id(trans_unit.operons_id)
                strand = operon.strand
                first_gene = get_first_gene_of_tu(regulated_genes, strand)
        if reg_int.regulatory_sites_id:
            reg_sites = multigenomic_api.regulatory_sites.find_by_id(reg_int.regulatory_sites_id)
            if reg_sites.absolute_position:
                if strand == "forward":
                    distance_to_first_gene = reg_sites.absolute_position - first_gene["left_end_position"]
                else:
                    distance_to_first_gene = first_gene["right_end_position"] - reg_sites.absolute_position
            elif reg_sites.left_end_position and reg_sites.right_end_position:
                abs_pos = (reg_sites.right_end_position - reg_sites.left_end_position) / 2 + reg_sites.left_end_position
                if first_gene:
                    if strand == "forward":
                        distance_to_first_gene = abs_pos - first_gene["left_end_position"]
                    else:
                        distance_to_first_gene = first_gene["right_end_position"] - abs_pos
            return {
                "name": first_gene["name"],
                "distanceTo": distance_to_first_gene
            }
        return {
            "name": first_gene["name"],
            "distanceTo": ""
        }
    return {
        "name": "",
        "distanceTo": ""
    }


def get_first_gene_of_tu(genes, strand):
    dict_genes = []
    first_gene = ""
    for gene in genes:
        gene_object = multigenomic_api.genes.find_by_id(gene.get("_id"))
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
            first_gene = (min(dict_genes, key=lambda x: x["left_end_position"]))
        elif strand == "reverse":
            first_gene = (max(dict_genes, key=lambda x: x["right_end_position"]))
    return first_gene


def get_regulated_genes(regulated_entity):
    regulated_genes = []
    transcription_units = []
    if regulated_entity.type == "gene":
        gene = multigenomic_api.genes.find_by_id(regulated_entity.id)
        gene_object = {
            "_id": gene.id,
            "name": gene.name,
        }
        regulated_genes.append(gene_object)
    elif regulated_entity.type == "promoter":
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

