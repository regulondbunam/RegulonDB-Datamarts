import multigenomic_api
from src.datamarts.domain.general.biological_base import BiologicalBase
from src.datamarts.domain.general.additiveEvidences import AdditiveEvidences


class RegulatoryInteractions(BiologicalBase):
    def __init__(self, reg_interaction):
        super().__init__([], reg_interaction.citations, "")
        self.active_conformation = reg_interaction.regulator
        self.regulatory_interaction = reg_interaction
        self.regulated_entity = reg_interaction.regulated_entity
        self.regulated_genes = reg_interaction.regulated_entity
        self.regulatory_binding_sites = reg_interaction

    def to_dict(self):
        reg_genes = self.regulated_genes
        citations = self.citations
        reg_bind_sites = self.regulatory_binding_sites
        additive_evs = AdditiveEvidences(citations + reg_bind_sites.get("citations", []))
        regulatory_interactions = {
            "_id": self.regulatory_interaction.id,
            "function": self.regulatory_interaction.function,
            "regulatedEntity": self.regulated_entity,
            "activeConformation": self.active_conformation,
            "distanceToFirstGene": get_distance_to_first_gene(self.regulatory_interaction, reg_genes),
            "distanceToPromoter": self.regulatory_interaction.dist_site_promoter,
            "regulatedGenes": self.regulated_genes,
            "regulatoryBindingSites":  reg_bind_sites,
            "citations": citations,
            # TODO: Check if this field is correctly obtained
            "mechanism": self.regulatory_interaction.mechanism,
            "additiveEvidences": additive_evs.to_dict(),
            "confidenceLevel": additive_evs.get_confidence_level()
        }
        return regulatory_interactions

    @property
    def active_conformation(self):
        return self._active_conformation

    @active_conformation.setter
    def active_conformation(self, regulator):
        name = regulator.name
        if regulator.type == "regulatoryComplex":
            reg = multigenomic_api.regulatory_complexes.find_by_id(regulator.id)
            if reg.abbreviated_name:
                name = reg.abbreviated_name
            else:
                tf = multigenomic_api.transcription_factors.find_by_name(regulator.name)
                if tf:
                    name = tf.abbreviated_name
        elif regulator.type == "product":
            reg = multigenomic_api.products.find_by_id(regulator.id)
            if reg.abbreviated_name:
                name = reg.abbreviated_name
        self._active_conformation = {
            "_id": regulator.id,
            "type": regulator.type,
            "name": name
        }

    @property
    def regulated_entity(self):
        return self._regulated_entity

    @regulated_entity.setter
    def regulated_entity(self, regulated_entity):
        self._regulated_entity = {
            "_id": regulated_entity.id,
            "type": regulated_entity.type,
            "name": regulated_entity.name
        }

    @property
    def regulated_genes(self):
        return self._regulated_genes

    @regulated_genes.setter
    def regulated_genes(self, regulated_entity):
        self._regulated_genes = []
        transcription_units = []
        strand = None
        if regulated_entity.type == "gene":
            self._regulated_genes.append({
                        "_id": regulated_entity.id,
                        "name": regulated_entity.name,
                    })
        elif regulated_entity.type == "promoter":
            transcription_units = multigenomic_api.transcription_units.find_by_promoter_id(regulated_entity.id)
            promoter = multigenomic_api.promoters.find_by_id(regulated_entity.id)
            strand = promoter.strand
        elif regulated_entity.type == "transcriptionUnit":
            trans_unit = multigenomic_api.transcription_units.find_by_id(regulated_entity.id)
            transcription_units.append(trans_unit)
            operon = multigenomic_api.operons.find_by_id(trans_unit.operons_id)
            strand = operon.strand
        if transcription_units:
            for tu in transcription_units:
                for gene_id in tu.genes_ids:
                    min_left_pos = None
                    max_right_pos = None
                    gene = multigenomic_api.genes.find_by_id(gene_id)
                    if gene.fragments:
                        min_left_pos = min(gene.fragments, key=lambda x: x.left_end_position)
                        max_right_pos = max(gene.fragments, key=lambda x: x.right_end_position)
                    gene_object = {
                        "_id": gene.id,
                        "name": gene.name,
                        "leftEndPosition": gene.left_end_position or min_left_pos.left_end_position,
                        "rightEndPosition": gene.right_end_position or max_right_pos.right_end_position
                    }
                    if gene_object not in self._regulated_genes:
                        self._regulated_genes.append(gene_object)
            if strand == "forward":
                self._regulated_genes = sorted(self._regulated_genes, key=lambda x: x["leftEndPosition"])
            else:
                self._regulated_genes = sorted(self._regulated_genes, key=lambda x: x["rightEndPosition"])

    @property
    def regulatory_binding_sites(self):
        return self._regulatory_binding_sites

    @regulatory_binding_sites.setter
    def regulatory_binding_sites(self, reg_interaction):
        self._regulatory_binding_sites = {}
        strand = ""
        if reg_interaction.regulatory_sites_id:
            if reg_interaction.regulated_entity.type == "promoter":
                promoter = multigenomic_api.promoters.find_by_id(reg_interaction.regulated_entity.id)
                strand = promoter.strand
            elif reg_interaction.regulated_entity.type == "transcriptionUnit":
                trans_unit = multigenomic_api.transcription_units.find_by_id(reg_interaction.regulated_entity.id)
                if trans_unit.promoters_id:
                    promoter = multigenomic_api.promoters.find_by_id(trans_unit.promoters_id)
                    strand = promoter.strand
            elif reg_interaction.regulated_entity.type == "gene":
                gene = multigenomic_api.genes.find_by_id(reg_interaction.regulated_entity.id)
                strand = gene.strand
            regulatory_sites = multigenomic_api.regulatory_sites.find_by_id(reg_interaction.regulatory_sites_id)
            reg_binding_sites_obj = RegulatoryBindingSites(regulatory_sites).to_dict()
            if strand == "reverse":
                reg_binding_sites_obj["sequence"] = reverse_complement(reg_binding_sites_obj["sequence"])
            reg_binding_sites_obj["strand"] = strand
            self._regulatory_binding_sites = reg_binding_sites_obj


def get_distance_to_first_gene(reg_int, regulated_genes):
    strand = None
    first_gene = None
    if reg_int.regulatory_sites_id:
        reg_sites = multigenomic_api.regulatory_sites.find_by_id(reg_int.regulatory_sites_id)
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
        if reg_sites.absolute_position:
            if first_gene:
                if strand == "forward":
                    return reg_sites.absolute_position - first_gene["left_end_position"]
                else:
                    return first_gene["right_end_position"] - reg_sites.absolute_position
        elif reg_sites.left_end_position and reg_sites.right_end_position:
            abs_pos = (reg_sites.right_end_position - reg_sites.left_end_position)/2 + reg_sites.left_end_position
            if first_gene:
                if strand == "forward":
                    return abs_pos - first_gene["left_end_position"]
                else:
                    return first_gene["right_end_position"] - abs_pos
    return None


class RegulatoryBindingSites(BiologicalBase):
    def __init__(self, reg_sites):
        super().__init__(reg_sites.external_cross_references, reg_sites.citations, reg_sites.note)
        self.reg_sites = reg_sites

    def to_dict(self):
        regulatory_binding_sites = {
            "_id": self.reg_sites.id,
            "absolutePosition": self.reg_sites.absolute_position,
            "leftEndPosition": self.reg_sites.left_end_position,
            "rightEndPosition": self.reg_sites.right_end_position,
            "sequence": self.reg_sites.sequence,
            "strand": self.reg_sites.strand,
            "citations": self.citations
        }
        return regulatory_binding_sites


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


def get_first_gene_of_tu(genes, strand):
    dict_genes = []
    first_gene = None
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
