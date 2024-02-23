import multigenomic_api

from src.datamarts.domain.general.biological_base import BiologicalBase

from src.datamarts.domain.operon_datamart.transcription_unit.promoters import Promoters
from src.datamarts.domain.operon_datamart.transcription_unit.terminators import Terminators

from src.datamarts.domain.operon_datamart.regulator_binding_sites import RegulatoryBindingSites
from src.datamarts.domain.general.additiveEvidences import AdditiveEvidences


class TranscriptionUnit(BiologicalBase):

    def __init__(self, trans_unit):
        super().__init__(trans_unit.external_cross_references, trans_unit.citations, trans_unit.note)
        regulatory_int = get_all_ri_from_tu(trans_unit)
        self.transcription_unit = trans_unit
        self.first_gene = trans_unit
        self.genes = trans_unit.genes_ids
        self.promoter = trans_unit
        self.regulator_binding_sites = trans_unit.id
        self.sites = regulatory_int
        self.terminators = trans_unit.terminators_ids
        self.transcription_factors = regulatory_int

    @property
    def first_gene(self):
        return self._first_gene

    @first_gene.setter
    def first_gene(self, transcription_unit):
        self._first_gene = {}
        dist = None
        if transcription_unit.promoters_id:
            promoter = multigenomic_api.promoters.find_by_id(transcription_unit.promoters_id)
            first_gene = get_first_gene_of_tu(transcription_unit, promoter)
            if promoter.transcription_start_site:
                if promoter.strand == "forward":
                    dist = first_gene["leftEndPosition"] - promoter.transcription_start_site.left_end_position
                if promoter.strand == "reverse":
                    dist = first_gene["rightEndPosition"] - promoter.transcription_start_site.right_end_position
                self._first_gene = {
                    "distanceToPromoter": dist,
                    "_id": first_gene["id"],
                    "name": first_gene["name"]
                }

    @property
    def genes(self):
        return self._genes

    @genes.setter
    def genes(self, gene_ids):
        self._genes = []
        for gene_id in gene_ids:
            gene = multigenomic_api.genes.find_by_id(gene_id)
            tf_binding_sites = RegulatoryBindingSites(gene_id)
            gene = {
                "_id": gene.id,
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
                trans_factor = multigenomic_api.transcription_factors.find_tf_id_by_conformation_id(ri.regulator.id)
                self._transcription_factors.append(trans_factor)
                trans_factor = multigenomic_api.transcription_factors.find_by_name(ri.regulator.name)
                self._transcription_factors.append(trans_factor)
        self._transcription_factors = set(list(self._transcription_factors))

    @property
    def regulator_binding_sites(self):
        return self._regulator_binding_sites

    @regulator_binding_sites.setter
    def regulator_binding_sites(self, transcription_unit_id):
        self._regulator_binding_sites = []
        tf_binding_sites_dict = RegulatoryBindingSites(transcription_unit_id)
        self._regulator_binding_sites = tf_binding_sites_dict.to_dict()

    def to_dict(self):
        citations = self.citations
        additive_evs = AdditiveEvidences(citations)
        transcription_unit = {
            "_id": self.transcription_unit.id,
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
            "regulatorBindingSites": self._regulator_binding_sites,
            "additiveEvidences": additive_evs.to_dict(),
            "confidenceLevel": additive_evs.get_confidence_level()
        }
        return transcription_unit


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


def get_all_ri_from_tu(trans_unit):
    ris = []
    ris.extend(
        multigenomic_api.regulatory_interactions.find_regulatory_interactions_by_reg_entity_id(trans_unit.id))
    if trans_unit.promoters_id:
        ris.extend(
            multigenomic_api.regulatory_interactions.find_regulatory_interactions_by_reg_entity_id(trans_unit.promoters_id))
    return ris
