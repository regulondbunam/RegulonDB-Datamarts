import multigenomic_api

from src.datamarts.domain.general.biological_base import BiologicalBase

from src.datamarts.domain.operon_datamart.transcription_unit.promoters import Promoters
from src.datamarts.domain.operon_datamart.transcription_unit.terminators import Terminators

from src.datamarts.domain.operon_datamart.regulator_binding_sites import Regulator_Binding_Sites


class TranscriptionUnit(BiologicalBase):

    def __init__(self, transcription_unit):
        super().__init__(transcription_unit.external_cross_references, transcription_unit.citations, transcription_unit.note)
        regulatory_int = multigenomic_api.regulatory_interactions.find_regulatory_interactions_by_reg_entity_id(transcription_unit.id)
        self.transcription_unit = transcription_unit
        self.first_gene = transcription_unit
        self.genes = transcription_unit.genes_ids
        self.promoter = transcription_unit
        self.regulator_binding_sites = transcription_unit.id
        self.sites = regulatory_int
        self.terminators = transcription_unit.terminators_ids
        self.transcription_factors = regulatory_int

    @property
    def first_gene(self):
        return self._first_gene

    @first_gene.setter
    def first_gene(self, transcription_unit):
        self._first_gene = {}
        if transcription_unit.promoters_id:
            promoter = multigenomic_api.promoters.find_by_id(transcription_unit.promoters_id)
            first_gene = get_first_gene_of_tu(transcription_unit, promoter)
            if promoter.transcription_start_site:
                self._first_gene = {
                    "distanceToPromoter": abs(first_gene["pos"] - promoter.transcription_start_site.left_end_position),
                    "gene_id": first_gene["id"],
                    "gene_name": first_gene["name"]
                }

    @property
    def genes(self):
        return self._genes

    @genes.setter
    def genes(self, gene_ids):
        self._genes = []
        for gene_id in gene_ids:
            gene = multigenomic_api.genes.find_by_id(gene_id)
            tf_binding_sites = Regulator_Binding_Sites(gene_id)
            gene = {
                "id": gene.id,
                "name": gene.name,
                "regulatorBindingSites": tf_binding_sites.to_dict()
            }
            self._genes.append(gene.copy())

    @property
    def promoter(self):
        return self._promoter

    @promoter.setter
    def promoter(self, transcription_unit):
        self._promoter = {}
        promoter_id = transcription_unit.promoters_id
        if promoter_id:
            promoter = Promoters(multigenomic_api.promoters.find_by_id(promoter_id))
            self._promoter = promoter.to_dict().copy()

    @property
    def sites(self):
        return self._sites

    @sites.setter
    def sites(self, reg_int):
        self._sites = []
        for ri in reg_int:
            if ri.regulatory_sites_id:
                self._sites.append(ri.regulatory_sites_id)


    @property
    def terminators(self):
        return self._terminators

    @terminators.setter
    def terminators(self, terminator_ids):
        self._terminators = []
        for terminator_id in terminator_ids:
            terminator = multigenomic_api.terminators.find_by_id(terminator_id)
            terminator_dict = Terminators(terminator).to_dict()
            self._terminators.append(terminator_dict.copy())

    @property
    def transcription_factors(self):
        return self._transcription_factors

    @transcription_factors.setter
    def transcription_factors(self, reg_int):
        self._transcription_factors = []
        for ri in reg_int:
            if ri.regulator:
                trans_factor = multigenomic_api.transcription_factors.find_tf_id_by_active_conformation_id(ri.regulator.id)
                self._transcription_factors.extend(trans_factor)
        self._transcription_factors = set(list(self._transcription_factors))


    @property
    def regulator_binding_sites(self):
        return self._regulator_binding_sites

    @regulator_binding_sites.setter
    def regulator_binding_sites(self, transcription_unit_id):
        self._regulator_binding_sites = []
        tf_binding_sites_dict = Regulator_Binding_Sites(transcription_unit_id)
        self._regulator_binding_sites = tf_binding_sites_dict.to_dict()

    def to_dict(self):
        transcription_unit = {
            "id": self.transcription_unit.id,
            "name": self.transcription_unit.name,
            "citations": self.citations,
            "firstGene": self.first_gene,
            "genes": self.genes,
            "note": self.formatted_note,
            "synonyms": self.transcription_unit.synonyms,
            "promoter": self.promoter,
            "statistics": {
                "genes": len(self.genes),
                "sites": len(self.sites),
                "transcriptionFactors": len(self.transcription_factors)
            },
            "terminators": self.terminators,
            "regulatorBindingSites": self._regulator_binding_sites
        }
        return transcription_unit


def get_first_gene_of_tu(transcription_unit, promoter):
    first_gene = multigenomic_api.genes.find_by_id(transcription_unit.genes_ids[0])
    if promoter.strand == "reverse":
        first_gene.left_end_position = first_gene.right_end_position
    first_gene.left_end_position = first_gene.left_end_position or first_gene.fragments[0].left_end_position

    for gene in transcription_unit.genes_ids:
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
        "name": first_gene.name,
        "pos": first_gene.left_end_position
    }
    return first_gene

