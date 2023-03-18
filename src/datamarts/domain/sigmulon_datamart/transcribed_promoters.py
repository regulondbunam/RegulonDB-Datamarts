from src.datamarts.domain.general.biological_base import BiologicalBase
import multigenomic_api


class TranscribedPromoters(BiologicalBase):

    def __init__(self, promoter):
        super().__init__(promoter.external_cross_references, promoter.citations, promoter.note)
        self.promoter = promoter
        self.transcribed_genes = promoter
        self.boxes = promoter.boxes
        self.operon_id = promoter.id

    def to_dict(self):
        tss_pos = None
        if self.promoter.transcription_start_site:
            tss_pos = self.promoter.transcription_start_site.left_end_position
        transcribed_promoters = {
            "_id": self.promoter.id,
            "name": self.promoter.name,
            "transcribedGenes": self.transcribed_genes,
            "operonId": self.operon_id,
            "sequence": self.promoter.sequence,
            "TSSPosition": tss_pos,
            "boxes": self.boxes,
            "citations": self.citations
        }
        return transcribed_promoters

    @property
    def transcribed_genes(self):
        return self._transcribed_genes

    @transcribed_genes.setter
    def transcribed_genes(self, promoter):
        self._transcribed_genes = []
        trans_units = multigenomic_api.transcription_units.find_by_promoter_id(promoter.id)
        for tu in trans_units:
            if tu.genes_ids:
                for gene_id in tu.genes_ids:
                    gene = multigenomic_api.genes.find_by_id(gene_id)
                    gene_obj = {
                        "_id": gene.id,
                        "name": gene.name,
                        "distanceFromTSS": abs(get_distance_from_tss(promoter, gene))
                    }
                    if gene_obj not in self._transcribed_genes:
                        self._transcribed_genes.append(gene_obj)

    @property
    def boxes(self):
        return self._boxes

    @boxes.setter
    def boxes(self, boxes):
        self._boxes = []
        for box in boxes:
            self._boxes.append({
                "leftEndPosition": box.left_end_position,
                "rightEndPosition": box.right_end_position,
                "sequence": box.sequence,
                "type": box.type
            })

    @property
    def operon_id(self):
        return self._operon_id

    @operon_id.setter
    def operon_id(self, promoter_id):
        self._operon_id = None
        operons_ids = []
        trans_units = multigenomic_api.transcription_units.find_by_promoter_id(promoter_id)
        for tu in trans_units:
            if tu.operons_id not in operons_ids:
                operons_ids.append(tu.operons_id)
        if len(operons_ids) == 1:
            self._operon_id = operons_ids[0]


def get_distance_from_tss(promoter, gene):
    distance = 0
    if not promoter.transcription_start_site:
        return distance
    if gene.strand == "forward":
        if gene.left_end_position:
            distance = promoter.transcription_start_site.left_end_position - gene.left_end_position
    else:
        if gene.right_end_position:
            distance = gene.right_end_position - promoter.transcription_start_site.right_end_position
    return distance
