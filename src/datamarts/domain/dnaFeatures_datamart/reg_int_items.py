import multigenomic_api


class RegIntDnaFeatures(object):
    def __init__(self, dict_colors):
        self.dict_colors = dict_colors

    @property
    def objects(self):
        # Regulatory Interactions Objects
        reg_int_objects = multigenomic_api.regulatory_interactions.get_all()
        for reg_int in reg_int_objects[0:10]:
            if reg_int.regulator.type == "product" or reg_int.regulator.type == "regulatoryComplex":
                dtt_datamart = RegIntDnaFeatures.DTTDatamart(reg_int, self.dict_colors)
                yield dtt_datamart
        del reg_int_objects

    class DTTDatamart:

        def __init__(self, entity, dict_colors):
            self.entity = entity
            self.dict_colors = dict_colors
            self.objectRGBColor = entity.function
            self.positions = entity
            self.related_genes = entity
            self.linked_object_when_no_positions = entity
            self.regulator = entity
            self.tooltip = entity

        @property
        def objectRGBColor(self):
            return self._color

        @objectRGBColor.setter
        def objectRGBColor(self, entity_function):
            self._color = ""
            if entity_function == "activator":
                self._color = "10, 255, 5"
            elif entity_function == "repressor":
                self._color = "243, 10, 2"

        @property
        def positions(self):
            return self._positions

        @positions.setter
        def positions(self, entity):
            self._positions = {
                "leftEndPosition": None,
                "rightEndPostion": None
            }
            reg_sites = multigenomic_api.regulatory_sites.find_by_id(entity.regulatory_sites_id)
            if reg_sites:
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

        @property
        def tooltip(self):
            return self._tooltip

        @tooltip.setter
        def tooltip(self, entity):
            self._tooltip = f"{self.regulator.name,}; \n"\
                            f"Distance to transcription start site: {self.entity.absolute_center_position}" \
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
                        location = ""
                        if gene.strand == "forward":
                            location = f"{gene.left_end_position}->{gene.right_end_position}"
                        else:
                            location = f"{gene.left_end_position}<-{gene.right_end_position}"
                        gene_object = {
                            "gene_id": gene_id,
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

        def to_dict(self):
            dttDatamart = {
                "_id": self.entity.id,
                "labelFont": "arial",
                "labelRGBColor": "0,0,0",
                "labelSize": 12,
                "labelName": self.regulator.abbreviated_name or self.regulator.name,
                "leftEndPosition": self.positions["leftEndPosition"],
                "lineRGBColor": "0,0,0",
                # TODO: this gonna be defined by evidence, if is weak or strong
                "lineType": 1,
                "lineWidth": 1,
                "linkedObjectWhenNoPositions": [],
                "objectType": "tf_binding_site",
                "objectRGBColor": self._color,
                "organism": {
                    "organism_id": self.entity.organisms_id
                },
                "relatedGenes": self.related_genes,
                "rightEndPosition": self.positions["rightEndPosition"],
                "strand": self.entity.strand,
                # TODO: ask what will be showed on this tooltip
                "tooltip": self.tooltip
            }
            return dttDatamart