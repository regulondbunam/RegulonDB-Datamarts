import multigenomic_api
from src.datamarts.domain.general.biological_base import BiologicalBase


class PromoterDnaFeatures(object):
    def __init__(self, dict_colors):
        self.dict_colors = dict_colors

    @property
    def objects(self):
        # Promoter Objects
        promoter_objects = multigenomic_api.promoters.get_all()
        for promoter in promoter_objects:
            dtt_datamart = PromoterDnaFeatures.DTTDatamart(promoter, self.dict_colors)
            yield dtt_datamart
        del promoter_objects

    class DTTDatamart(BiologicalBase):

        def __init__(self, entity, dict_colors):
            super().__init__(entity.external_cross_references, entity.citations, entity.note)
            self.entity = entity
            self.product = entity.id
            self.positions = entity.transcription_start_site
            self.dict_colors = dict_colors
            self.linked_object_when_no_positions = entity
            self.related_genes = entity
            self.line_type = entity.citations

        @property
        def product(self):
            return self._product

        @product.setter
        def product(self, entity_id):
            self._product = ""
            products = multigenomic_api.products.find_by_gene_id(entity_id)
            if products:
                product_tooltip = products[0].name
                if len(products) > 1:
                    for product in products:
                        product_tooltip = f"{product_tooltip}, {product.name}"
                self._product = product_tooltip

        @property
        def positions(self):
            return self._positions

        @positions.setter
        def positions(self, transcription_start_site):
            self._positions = {
                    "leftEndPosition": None,
                    "rightEndPosition": None
                }
            if transcription_start_site:
                self._positions = {
                    "leftEndPosition": transcription_start_site.left_end_position,
                    "rightEndPosition": transcription_start_site.right_end_position
                }

        @property
        def linked_object_when_no_positions(self):
            return self._linked_object_when_no_positions

        @linked_object_when_no_positions.setter
        def linked_object_when_no_positions(self, entity):
            self._linked_object_when_no_positions = {}
            if not entity.transcription_start_site:
                trans_units = multigenomic_api.transcription_units.find_by_promoter_id(entity.id)
                gene = get_first_gene_of_tu(trans_units, entity.strand)
                if gene.fragments:
                    gene = get_first_fragment_position(gene)
                self._linked_object_when_no_positions = {
                    "_id": gene.id,
                    "leftEndPosition": gene.left_end_position,
                    "name": gene.name,
                    "rightEndPosition": gene.right_end_position,
                    "strand": gene.strand,
                    "type": "gene"
                }

        @property
        def related_genes(self):
            return self._related_genes

        @related_genes.setter
        def related_genes(self, entity):
            self._related_genes = []
            trans_units = multigenomic_api.transcription_units.find_by_promoter_id(entity.id)
            for tu in trans_units:
                for gene_id in tu.genes_ids:
                    gene = multigenomic_api.genes.find_by_id(gene_id)
                    terms = multigenomic_api.terms.get_terms_by_gene_member_id(gene_id)
                    color = ""
                    if terms:
                        color = self.dict_colors[terms[0].id]
                    location = ""
                    if gene.strand == "forward":
                        location = f"{gene.left_end_position}->{gene.right_end_position}"
                    else:
                        location = f"{gene.left_end_position}<-{gene.right_end_position}"
                    gene_object = {
                        "_id": gene_id,
                        "effect": "",
                        "name": gene.name,
                        "objectRGBColor": color,
                        "strand": gene.strand,
                        "tooltip": f"Gene: {gene.name};\n"
                                   f"Synonyms: {gene.synonyms}\n"
                                   f"Location: {location};\n"
                    }
                    if gene_object not in self._related_genes:
                        self._related_genes.append(gene_object)

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

        def to_dict(self):
            dttDatamart = {
                "_id": self.entity.id,
                "labelFont": "arial",
                "labelRGBColor": "0,0,0",
                "labelSize": 12,
                "labelName": self.entity.name,
                "citations": self.citations,
                "leftEndPosition": self.positions["leftEndPosition"],
                "lineRGBColor": "0,0,0",
                # TODO: this gonna be defined by evidence, if is weak or strong
                "lineType": self.line_type,
                "lineWidth": 1,
                "linkedObjectWhenNoPositions": self.linked_object_when_no_positions,
                "objectType": "promoter",
                "objectRGBColor": "0,0,0",
                "organism": {
                    "_id": self.entity.organisms_id
                },
                "relatedGenes": self.related_genes,
                "rightEndPosition": self.positions["rightEndPosition"],
                "strand": self.entity.strand,
                # TODO: ask what will be showed on this tooltip
                "tooltip": f"Promoter: {self.entity.name}; \n"
                           f"Tr. Start site: {self.positions['leftEndPosition']}->{self.positions['rightEndPosition']};"
            }
            return dttDatamart


def get_first_fragment_position(gene):
    gene.left_end_position = gene.fragments[0].left_end_position
    for fragment in gene.fragments:
        if fragment.left_end_position < gene.left_end_position:
            gene.left_end_position = fragment.left_end_position
            gene.right_end_position = fragment.right_end_position
    return gene


def get_first_gene_of_tu(trans_units, strand):
    first_gene = multigenomic_api.genes.find_by_id(trans_units[0].genes_ids[0])
    if first_gene.fragments:
        first_gene = get_first_fragment_position(first_gene)
    first_gene = {
        "gene_id": first_gene.id,
        "leftEndPosition": first_gene.left_end_position
    }
    for tu in trans_units:
        for gene_id in tu.genes_ids:
            gene = multigenomic_api.genes.find_by_id(gene_id)
            if gene.fragments:
                gene = get_first_fragment_position(gene)
            if strand == "forward":
                if first_gene["leftEndPosition"] > gene.left_end_position:
                    first_gene = {
                        "gene_id": gene.id,
                        "leftEndPosition": gene.left_end_position
                    }
            elif strand == "reverse":
                if first_gene["leftEndPosition"] < gene.left_end_position:
                    first_gene = {
                        "gene_id": gene.id,
                        "leftEndPosition": gene.left_end_position
                    }
    return multigenomic_api.genes.find_by_id(first_gene["gene_id"])
