import multigenomic_api
from src.datamarts.domain.general.biological_base import BiologicalBase
from src.datamarts.domain.general.additiveEvidences import AdditiveEvidences


class RegulatoryInteractions(BiologicalBase):
    def __init__(self, reg_interaction):
        super().__init__([], reg_interaction.citations, "")
        self.regulatory_interaction = reg_interaction
        self.regulator = reg_interaction.regulator
        self.regulated_entity = reg_interaction.regulated_entity
        self.regulated_genes = reg_interaction.regulated_entity
        self.regulatory_binding_sites = reg_interaction

    def to_dict(self):
        reg_genes = self.regulated_genes
        distance_to_first_gene = get_distance_to_first_gene(self.regulatory_interaction, reg_genes)
        citations = self.citations
        reg_bind_sites = self.regulatory_binding_sites
        additive_evs = AdditiveEvidences(citations + reg_bind_sites.get("citations", []))
        regulatory_interactions = {
            "_id": self.regulatory_interaction.id,
            "regulator": self.regulator,
            "function": self.regulatory_interaction.function,
            "regulatedEntity": self.regulated_entity,
            "distanceToFirstGene": distance_to_first_gene,
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
    def regulator(self):
        return self._regulator

    @regulator.setter
    def regulator(self, regulator):
        self._regulator = {
            "_id": regulator.id,
            "type": regulator.type,
            "name": regulator.name,
        }

    @property
    def regulated_entity(self):
        return self._regulated_entity

    @regulated_entity.setter
    def regulated_entity(self, regulated_entity):
        self._regulated_entity = {}
        reg_entity = {}
        if regulated_entity.type == "promoter":
            reg_entity = multigenomic_api.promoters.find_by_id(regulated_entity.id)
        elif regulated_entity.type == "transcriptionUnit":
            reg_entity = multigenomic_api.transcription_units.find_by_id(regulated_entity.id)
        if reg_entity:
            self._regulated_entity = {
                "_id": regulated_entity.id,
                "type": regulated_entity.type,
                "name": reg_entity["name"]
            }

    @property
    def regulated_genes(self):
        return self._regulated_genes

    @regulated_genes.setter
    def regulated_genes(self, regulated_entity):
        self._regulated_genes = []
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
                    if gene_object not in self._regulated_genes:
                        self._regulated_genes.append(gene_object)

    @property
    def regulatory_binding_sites(self):
        return self._regulatory_binding_sites

    @regulatory_binding_sites.setter
    def regulatory_binding_sites(self, reg_interaction):
        self._regulatory_binding_sites = {}
        promoter = None
        if reg_interaction.regulatory_sites_id:
            if reg_interaction.regulated_entity.type == "promoter":
                promoter = multigenomic_api.promoters.find_by_id(reg_interaction.regulated_entity.id)
            elif reg_interaction.regulatory_sites_id == "transcriptionUnit":
                trans_unit = multigenomic_api.transcription_units.find_by_id(reg_interaction.regulated_entity_id)
                if trans_unit.promoters_id:
                    promoter = multigenomic_api.promoters.find_by_id(trans_unit.promoters_id)
            regulatory_sites = multigenomic_api.regulatory_sites.find_by_id(reg_interaction.regulatory_sites_id)
            reg_binding_sites_obj = RegulatoryBindingSites(regulatory_sites).to_dict()
            if promoter:
                if promoter.strand == "reverse":
                    reg_binding_sites_obj["sequence"] = reverse_complement(reg_binding_sites_obj["sequence"])
                reg_binding_sites_obj["strand"] = promoter.strand
            self._regulatory_binding_sites = reg_binding_sites_obj


def get_distance_to_first_gene(reg_int, regulated_genes):
    distance_to_first_gene = None
    promoter = None
    if reg_int.regulated_entity.type == "promoter":
        promoter = multigenomic_api.promoters.find_by_id(reg_int.regulated_entity.id)
    elif reg_int.regulated_entity.type == "transcriptionUnit":
        trans_unit = multigenomic_api.transcription_units.find_by_id(reg_int.regulated_entity.id)
        if trans_unit.promoters_id:
            promoter = multigenomic_api.promoters.find_by_id(trans_unit.promoters_id)
    first_gene = get_first_gene_of_tu(regulated_genes, promoter)
    if reg_int.regulatory_sites_id:
        reg_sites = multigenomic_api.regulatory_sites.find_by_id(reg_int.regulatory_sites_id)
        if reg_sites and promoter:
            if reg_sites.absolute_position:
                if promoter.strand == "forward":
                    distance_to_first_gene = reg_sites.absolute_position - first_gene["leftEndPosition"]
                else:
                    distance_to_first_gene = first_gene["leftEndPosition"] - reg_sites.absolute_position

    return distance_to_first_gene


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
            "leftEndPosition": first_gene.left_end_position
        }
    return first_gene
