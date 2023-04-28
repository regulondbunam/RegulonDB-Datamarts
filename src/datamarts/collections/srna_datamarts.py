import multigenomic_api

from src.datamarts.domain.srna_datamart.summary import Summary

from src.datamarts.domain.srna_datamart.regulatory_interactions import RegulatoryInteractions
from src.datamarts.domain.general.biological_base import BiologicalBase

from src.datamarts.domain.general.remove_items import remove_empty_items


class SrnaDatamarts:

    @property
    def objects(self):
        srna_objects = multigenomic_api.products.get_all_srnas()
        for srna in srna_objects:
            print(srna.id)
            if srna.type == "small RNA":
                srna_datamart = SrnaDatamarts.SrnaDatamart(srna)
                yield srna_datamart
        del srna_objects

    class SrnaDatamart(BiologicalBase):
        def __init__(self, srna):
            super().__init__(srna.external_cross_references, srna.citations, srna.note)
            self.id = srna.id
            self.product = srna
            self.regulatory_interactions = srna
            self.summary = srna

        @property
        def product(self):
            return self._small_rna

        @product.setter
        def product(self, srna):
            gene = multigenomic_api.genes.find_by_id(srna.genes_id)
            genes = {
                "_id": gene.id,
                "name": gene.name,
                "leftEndPosition": gene.left_end_position,
                "rightEndPosition": gene.right_end_position,
                "strand": gene.strand,
                "gcContent": gene.gc_content
            }
            self._small_rna = {
                "_id": srna.id,
                "gene": genes,
                "name": srna.name,
                "synonyms": srna.synonyms,
                "note": self.formatted_note,
                "sequence": srna.sequence,
                "citations": self.citations,
                "externalCrossReferences": self.external_cross_references
            }

        @property
        def regulatory_interactions(self):
            return self._regulatory_interactions

        @regulatory_interactions.setter
        def regulatory_interactions(self, srna):
            self._regulatory_interactions = []
            reg_ints = multigenomic_api.regulatory_interactions.find_by_regulator_id(srna.id)
            for ri in reg_ints:
                self._regulatory_interactions.append(RegulatoryInteractions(ri).to_dict())

        @property
        def summary(self):
            return self._summary

        @summary.setter
        def summary(self, srna):
            self._summary = {}
            summary = Summary(srna, self.regulatory_interactions).to_dict()
            self._summary = summary

        def to_dict(self):
            srna_datamart = {
                "_id": self.id,
                "product": self.product,
                "regulatoryInteractions": self.regulatory_interactions,
                "summary": self.summary,
                "allCitations": BiologicalBase.get_all_citations()
            }
            return srna_datamart


def all_srna_datamarts():
    srna_items = SrnaDatamarts()
    json_srna = []
    for srna in srna_items.objects:
        srna_dict = remove_empty_items(srna.to_dict().copy())
        json_srna.append(remove_empty_items(srna_dict))
    return json_srna
