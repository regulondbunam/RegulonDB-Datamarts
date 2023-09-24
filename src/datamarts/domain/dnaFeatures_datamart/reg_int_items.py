import multigenomic_api
import re
from src.datamarts.domain.general.biological_base import BiologicalBase


class RegIntDnaFeatures(object):
    def __init__(self, dict_colors):
        self.dict_colors = dict_colors

    @property
    def objects(self):
        # Regulatory Interactions Objects
        reg_int_objects = multigenomic_api.regulatory_interactions.get_all()
        reg_int_objects = find_dual_reg_ints(reg_int_objects)
        for reg_int in reg_int_objects:
            if reg_int.regulator:
                print(reg_int.id)
                dtt_datamart = RegIntDnaFeatures.DTTDatamart(reg_int, self.dict_colors)
                yield dtt_datamart
        del reg_int_objects

    class DTTDatamart(BiologicalBase):

        def __init__(self, entity, dict_colors):
            super().__init__(entity.external_cross_references, entity.citations, entity.note)
            self.entity = entity
            self.dict_colors = dict_colors
            self.object_rgb_color = entity.function
            self.positions = entity
            self.related_genes = entity
            self.linked_object_when_no_positions = entity
            self.regulator = entity
            self.tooltip = entity
            self.strand = entity.regulated_entity
            self.line_type = entity.citations
            self.object_type = entity
            self.name = self.object_type

        @property
        def object_rgb_color(self):
            return self._color

        @object_rgb_color.setter
        def object_rgb_color(self, entity_function):
            self._color = ""
            if entity_function == "activator":
                self._color = "10, 255, 5"
            elif entity_function == "repressor":
                self._color = "243, 10, 2"
            elif entity_function == "dual":
                self._color = "163, 167, 167"
            else:
                self._color = "165, 166, 166"

        @property
        def positions(self):
            return self._positions

        @positions.setter
        def positions(self, entity):
            self._positions = {
                "leftEndPosition": None,
                "rightEndPosition": None
            }
            if entity.regulatory_sites_id:
                reg_sites = multigenomic_api.regulatory_sites.find_by_id(entity.regulatory_sites_id)
                self._positions = {
                    "leftEndPosition": reg_sites.left_end_position,
                    "rightEndPosition": reg_sites.right_end_position
                }

        @property
        def regulator(self):
            return self._regulator

        @regulator.setter
        def regulator(self, entity):
            if entity.regulator.type == "product":
                self._regulator = multigenomic_api.products.find_by_id(entity.regulator.id)
            elif entity.regulator.type == "regulatoryComplex":
                self._regulator = multigenomic_api.regulatory_complexes.find_by_id(entity.regulator.id)
            elif entity.regulator.type == "regulatoryContinuant":
                self._regulator = entity.regulator

        @property
        def tooltip(self):
            return self._tooltip

        @tooltip.setter
        def tooltip(self, entity):
            self._tooltip = f"{self.regulator.name,}; \n"\
                            f"Distance to transcription start site: {entity.dist_site_promoter} \n"\
                            f"Evidence:"  # TODO: add evidences

        @property
        def related_genes(self):
            return self._related_genes

        @related_genes.setter
        def related_genes(self, entity):
            self._related_genes = []
            if entity.regulated_entity.type == "promoter":
                trans_units = multigenomic_api.transcription_units.find_by_promoter_id(entity.regulated_entity.id)
                for tu in trans_units:
                    for gene_id in tu.genes_ids:
                        gene = multigenomic_api.genes.find_by_id(gene_id)
                        terms = multigenomic_api.terms.get_terms_by_gene_member_id(gene_id)
                        color = ""
                        if terms:
                            color = self.dict_colors[terms[0].id]
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
        def linked_object_when_no_positions(self):
            return self._linked_object_when_no_positions

        @linked_object_when_no_positions.setter
        def linked_object_when_no_positions(self, entity):
            self._linked_object_when_no_positions = {}
            if self.positions["leftEndPosition"] is None and self.positions["rightEndPosition"] is None:
                if entity.regulated_entity.type == "promoter":
                    promoter = multigenomic_api.promoters.find_by_id(entity.regulated_entity.id)
                    if promoter.transcription_start_site:
                        self._linked_object_when_no_positions = {
                            "_id": promoter.id,
                            "leftEndPosition": promoter.transcription_start_site.left_end_position or None,
                            "name": promoter.name,
                            "rightEndPosition": promoter.transcription_start_site.right_end_position or None,
                            "strand": promoter.strand,
                            "type": "promoter"
                        }
                    else:
                        self._linked_object_when_no_positions = linked_gene_to_promoter(promoter)
                elif entity.regulated_entity.type == "transcriptionUnit":
                    tu = multigenomic_api.transcription_units.find_by_id(entity.regulated_entity.id)
                    gene = multigenomic_api.genes.find_by_id(tu.genes_ids[0])
                    self._linked_object_when_no_positions = {
                        "_id": gene.id,
                        "leftEndPosition": gene.left_end_position or None,
                        "name": gene.name,
                        "rightEndPosition": gene.right_end_position or None,
                        "strand": gene.strand,
                        "type": "gene"
                    }
                elif entity.regulated_entity.type == "gene":
                    gene = multigenomic_api.genes.find_by_id(entity.regulated_entity.id)
                    self._linked_object_when_no_positions = {
                        "_id": gene.id,
                        "leftEndPosition": gene.left_end_position or None,
                        "name": gene.name,
                        "rightEndPosition": gene.right_end_position or None,
                        "strand": gene.strand,
                        "type": "gene"
                    }

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
        def object_type(self):
            return self._object_type

        @object_type.setter
        def object_type(self, reg_int):
            self._object_type = "tf_binding_site"
            if reg_int.regulation_type == "sRNA-Regulation":
                self._object_type = "srna"
            elif reg_int.regulation_type == "Protein-Regulation":
                self._object_type = "translational_tf_binding_site"
            elif reg_int.regulator.type == "regulatoryContinuant" or reg_int.regulator.name == "DksA-ppGpp":
                self._object_type = "ppGpp"

        @property
        def name(self):
            return self._name

        @name.setter
        def name(self, obj_type):
            self._name = self.entity.regulator.name
            if obj_type == "ppGpp":
                self._name = self.regulator.name
            elif obj_type == 'translational_tf_binding_site':
                if self.entity.regulator.type == "regulatoryComplex":
                    reg_complex = multigenomic_api.regulatory_complexes.find_by_id(self.regulator.id)
                    product = multigenomic_api.products.find_by_id(reg_complex.products[0].products_id)
                    self._name = reg_complex.abbreviated_name or product.abbreviated_name
                elif self.entity.regulator.type == "product":
                    product = multigenomic_api.products.find_by_id(self.regulator.id)
                    self._name = product.abbreviated_name or product.name
            elif obj_type == "srna":
                product = multigenomic_api.products.find_by_id(self.regulator.id)
                self._name = product.abbreviated_name or product.name
            else:
                if self.regulator.abbreviated_name:
                    self._name = self.regulator.abbreviated_name
                tf = multigenomic_api.transcription_factors.find_tf_id_by_conformation_id(self.regulator.id)
                if len(tf) > 0:
                    self._name = tf[0].abbreviated_name
                else:
                    tf = multigenomic_api.transcription_factors.find_by_name(self.regulator.name)
                    if tf:
                        self._name = tf[0].abbreviated_name

        @property
        def strand(self):
            return self._strand

        @strand.setter
        def strand(self, regulated_entity):
            if regulated_entity.type == "promoter":
                promoter = multigenomic_api.promoters.find_by_id(regulated_entity.id)
                self._strand = promoter.strand
            elif regulated_entity.type == "transcriptionUnit":
                trans_unit = multigenomic_api.transcription_units.find_by_id(regulated_entity.id)
                if trans_unit.operons_id:
                    operon = multigenomic_api.operons.find_by_id(trans_unit.operons_id)
                    self._strand = operon.strand
                else:
                    gene = multigenomic_api.genes.find_by_id(trans_unit.genes_ids)
                    self._strand = gene.strand
            elif regulated_entity.type == "gene":
                gene = multigenomic_api.genes.find_by_id(regulated_entity.id)
                self._strand = gene.strand

        def to_dict(self):
            dtt_datamart = {
                "_id": self.entity.id,
                "labelFont": "arial",
                "labelRGBColor": "0,0,0",
                "labelSize": 12,
                "labelName": self.name,
                "citations": self.citations,
                "leftEndPosition": self.positions["leftEndPosition"],
                "lineRGBColor": "0,0,0",
                "lineType": self.line_type,
                "lineWidth": 1,
                "linkedObjectWhenNoPositions": self.linked_object_when_no_positions,
                "objectType": self.object_type,
                "objectRGBColor": self._color,
                "organism": {
                    "_id": self.entity.organisms_id
                },
                "relatedGenes": self.related_genes,
                "rightEndPosition": self.positions["rightEndPosition"],
                "strand": self.strand,
                # TODO: ask what will be showed on this tooltip
                "tooltip": self.tooltip
            }
            return dtt_datamart


def find_dual_reg_ints(reg_ints_items):
    final_reg_items = []
    for reg_int in reg_ints_items:
        final_reg_items.append(reg_int)
        for rev_reg_int in final_reg_items:
            if reg_int.absolute_center_position == rev_reg_int.absolute_center_position and \
                    reg_int.regulatory_sites_id == rev_reg_int.regulatory_sites_id:
                if reg_int.function == "repressor" and rev_reg_int.function == "activator":
                    rev_reg_int.function = "dual"
                    final_reg_items.remove(reg_int)
                    break
                elif reg_int.function == "activator" and rev_reg_int.function == "repressor":
                    rev_reg_int.function = "dual"
                    final_reg_items.remove(reg_int)
                    break
    return final_reg_items


def linked_gene_to_promoter(promoter):
    trans_units = multigenomic_api.transcription_units.find_by_promoter_id(promoter.id)
    gene = get_first_gene_of_tu(trans_units, promoter.strand)
    if gene.fragments:
        gene = get_first_fragment_position(gene)
    return {
        "_id": gene.id,
        "leftEndPosition": gene.left_end_position,
        "name": gene.name,
        "rightEndPosition": gene.right_end_position,
        "strand": gene.strand,
        "type": "gene"
    }


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
