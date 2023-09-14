import multigenomic_api
import numpy as np
from src.datamarts.domain.general.biological_base import BiologicalBase


class TerminatorDNAFeatures:

    @property
    def objects(self):
        # Terminator Objects
        terminator_objects = multigenomic_api.terminators.get_all()
        for terminator in terminator_objects:
            dtt_datamart = TerminatorDNAFeatures.DTTDatamart(terminator)
            yield dtt_datamart
        del terminator_objects

    class DTTDatamart(BiologicalBase):

        def __init__(self, entity):
            super().__init__(entity.external_cross_references, entity.citations, entity.note)
            self.entity = entity
            self.line_type = entity.citations
            self.strand = entity

        @property
        def line_type(self):
            return self._line_type

        @line_type.setter
        def line_type(self, citations):
            self._line_type = 1
            if citations:
                for citation in citations:
                    if citation.evidences_id:
                        evidence = multigenomic_api.evidences.find_by_id(citation.evidences_id)
                        if evidence.type == "S":
                            self._line_type = 2
                            break
                        elif evidence.type == "C":
                            self._line_type = 3
                            break

        @property
        def strand(self):
            return self._strand

        @strand.setter
        def strand(self, entity):
            self._strand = ""
            trans_units = multigenomic_api.transcription_units.find_by_terminator_id(entity.id)
            if len(trans_units) > 1:
                print(entity.id)
                gene = get_common_gene_of_tu_list(trans_units)
                if len(gene) > 0:
                    gene = gene.pop()
                else:
                    gene = trans_units[0].genes_ids[0]
            else:
                gene = trans_units[0].genes_ids[0]
            gene = multigenomic_api.genes.find_by_id(gene)
            self._strand = gene.strand

        def to_dict(self):
            dttDatamart = {
                "_id": self.entity.id,
                "labelFont": "arial",
                "labelRGBColor": "0,0,0",
                "labelSize": 12,
                "citations": self.citations,
                "leftEndPosition": self.entity.transcriptionTerminationSite.left_end_position,
                "lineRGBColor": "0,0,0",
                "lineType": self._line_type,
                "lineWidth": 1,
                "objectType": "terminator",
                "objectRGBColor": "0,0,0",
                "organism": {
                    "_id": self.entity.organisms_id
                },
                "strand": self.strand,
                "rightEndPosition": self.entity.transcriptionTerminationSite.right_end_position,
                "tooltip": f"{self.entity.class_}; \n"
                           f"Genome Position: "
                           f"{self.entity.transcriptionTerminationSite.left_end_position}-"
                           f"{self.entity.transcriptionTerminationSite.right_end_position}; "
            }
            return dttDatamart


def get_common_gene_of_tu_list(trans_units_list):
    list_of_genes = []
    for tu in trans_units_list:
        list_of_genes.append(tu.genes_ids)
    return set(list_of_genes[0]).intersection(*list_of_genes)