import multigenomic_api
from src.datamarts.domain.general.biological_base import BiologicalBase


class RegulatoryInteractions(BiologicalBase):
    def __init__(self, reg_int):
        super().__init__(reg_int.external_cross_references, reg_int.citations, reg_int.note)
        self.reg_int = reg_int
        self.regulated_entity = reg_int.regulated_entity
        self.regulatory_binding_sites = reg_int.regulatory_sites_id
        self.distance_to_gene = reg_int

    @property
    def regulated_entity(self):
        return self._regulated_entity

    @regulated_entity.setter
    def regulated_entity(self, reg_ent):
        reg_entity_obj = {}
        if reg_ent.type == "gene":
            reg_entity_obj = multigenomic_api.genes.find_by_id(reg_ent.id)
        if reg_ent.type == "transcriptionUnit":
            reg_entity_obj = multigenomic_api.transcription_units.find_by_id(reg_ent.id)
        if reg_ent.type == "promoter":
            reg_entity_obj = multigenomic_api.promoters.find_by_id(reg_ent.id)
        self._regulated_entity = {
            "_id": reg_ent.id,
            "name": reg_entity_obj.name,
            "type": reg_ent.type
        }

    @property
    def regulatory_binding_sites(self):
        return self._regulatory_binding_sites

    @regulatory_binding_sites.setter
    def regulatory_binding_sites(self, reg_site_id):
        self._regulatory_binding_sites = {}
        if reg_site_id:
            reg_site = multigenomic_api.regulatory_sites.find_by_id(reg_site_id)
            self._regulatory_binding_sites = {
                "absolutePosition": reg_site.absolute_position,
                "leftEndPosition": reg_site.left_end_position,
                "rightEndPosition": reg_site.right_end_position,
                "sequence": dna_seq_to_rna_seq(reg_site.sequence)
            }

    @property
    def distance_to_gene(self):
        return self._distance_to_gene

    @distance_to_gene.setter
    def distance_to_gene(self, reg_int):
        self._distance_to_gene = 0
        if reg_int.regulatory_sites_id:
            reg_site = multigenomic_api.regulatory_sites.find_by_id(reg_int.regulatory_sites_id)
            if reg_int.regulated_entity.type == "gene":
                gene = multigenomic_api.genes.find_by_id(reg_int.regulated_entity.id)
                self._distance_to_gene = ((reg_site.left_end_position + reg_site.right_end_position)/2) - gene.left_end_position
            elif reg_int.regulated_entity.type == "transcriptionUnit":
                trans_unit = multigenomic_api.transcription_units.find_by_id(reg_int.regulated_entity.id)
                if len(trans_unit.genes_ids) > 1:
                    if trans_unit.promoters_id:
                        promoter = multigenomic_api.promoters.find_by_id(trans_unit.promoters_id)
                        first_gene = get_first_gene_of_tu(trans_unit.genes_ids, promoter)
                        self._distance_to_gene = ((reg_site.left_end_position + reg_site.right_end_position)/2) - first_gene["leftEndPosition"]
                else:
                    gene = multigenomic_api.genes.find_by_id(trans_unit.genes_ids[0])
                    # TODO: add a way to identify when gene has fragments instead left and right end positions
                    if gene.left_end_position:
                        self._distance_to_gene = ((reg_site.left_end_position + reg_site.right_end_position) / 2) - gene.left_end_position

    def to_dict(self):
        regulatory_interactions = {
            "_id": self.reg_int.id,
            "function": self.reg_int.function,
            "regulatedEntity": self.regulated_entity,
            "distanceToGene": self.distance_to_gene,
            "regulatoryBindingSites": self.regulatory_binding_sites,
            "mechanism": self.reg_int.mechanism,
            "citations": self.citations
        }
        return regulatory_interactions


def dna_seq_to_rna_seq(dna_seq):
    if dna_seq:
        alt_map = {'ins': '0'}
        complement = {'T': 'U', 't': 'u'}
        for k, v in alt_map.items():
            dna_seq = dna_seq.replace(k, v)
        bases = list(dna_seq)
        bases = [complement.get(base, base) for base in bases]
        bases = ''.join(bases)
        for k, v in alt_map.items():
            bases = bases.replace(v, k)
        return bases


def get_first_gene_of_tu(genes, promoter):
    first_gene = multigenomic_api.genes.find_by_id(genes[0])
    if promoter.strand == "reverse":
        first_gene.left_end_position = first_gene.right_end_position
    first_gene.left_end_position = first_gene.left_end_position or first_gene.fragments[0].left_end_position

    for gene in genes:
        current_gene = multigenomic_api.genes.find_by_id(gene)
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
