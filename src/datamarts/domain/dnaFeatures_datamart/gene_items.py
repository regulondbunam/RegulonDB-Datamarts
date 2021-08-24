import multigenomic_api


class GeneDnaFeatures(object):
    def __init__(self, dict_colors):
        self.dict_colors = dict_colors
    @property
    def objects(self):
        # Gene Objects
        gene_objects = multigenomic_api.genes.get_all()
        for gene in gene_objects:
            if gene.left_end_position:
                dtt_datamart = GeneDnaFeatures.DTTDatamart(gene, self.dict_colors)
                yield dtt_datamart
        del gene_objects

    class DTTDatamart:

        def __init__(self, entity, dict_colors):
            self.dict_colors = dict_colors
            self.entity = entity
            self.color = entity
            self.product = entity.id

        @property
        def color(self):
            return self._color

        @color.setter
        def color(self, entity):
            self._color = ""
            terms = multigenomic_api.terms.get_terms_by_gene_member_id(entity.id)
            if terms:
                self._color = self.dict_colors[terms[0].id]

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

        def to_dict(self):
            location = ""
            if self.entity.strand == "forward":
                location = f"{self.entity.left_end_position}->{self.entity.right_end_position}"
            else:
                location = f"{self.entity.left_end_position}<-{self.entity.right_end_position}"
            dttDatamart = {
                "_id": self.entity.id,
                "labelFont": "arial",
                "labelRGBColor": "0,0,0",
                "labelSize": 12,
                "labelName": self.entity.name,
                "leftEndPosition": self.entity.left_end_position,
                "lineRGBColor": "0,0,0",
                # TODO: this gonna be defined by evidence, if is weak or strong
                "lineType": 1,
                "lineWidth": 1,
                "objectType": "gene",
                "objectRGBColor": self.color,
                "organism": {
                    "organism_id": self.entity.organisms_id
                },
                "rightEndPosition": self.entity.right_end_position,
                "strand": self.entity.strand,
                "tooltip": f"Gene: {self.entity.name};\n"
                           f"Synonyms: {self.entity.synonyms}\n"
                           f"Location: {location};\n"
                           f"Product: {self.product}"
            }
            return dttDatamart

