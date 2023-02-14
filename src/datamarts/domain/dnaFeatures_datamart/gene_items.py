import multigenomic_api
from string import digits


class GeneDnaFeatures(object):
    def __init__(self, dict_colors):
        self.dict_colors = dict_colors

    @property
    def objects(self):
        # Gene Objects
        gene_objects = multigenomic_api.genes.get_all()
        for gene in gene_objects:
            print(gene.id)
            if gene.left_end_position:
                dtt_datamart = GeneDnaFeatures.DTTDatamart(gene, self.dict_colors, False)
                yield dtt_datamart
            if gene.fragments:
                frag_count = 0
                for fragment in gene.fragments:
                    gene.left_end_position = fragment.left_end_position
                    gene.right_end_position = fragment.right_end_position
                    gene.name = fragment.name if fragment.name is not None else fragment.id
                    dtt_datamart = GeneDnaFeatures.DTTDatamart(gene, self.dict_colors, True)
                    yield dtt_datamart
                    frag_count += 1
        del gene_objects

    class DTTDatamart:

        def __init__(self, entity, dict_colors, is_fragment):
            self.is_fragment = is_fragment
            self.dict_colors = dict_colors
            self.entity = entity
            self.color = entity
            self.product = entity.id
            self.location = entity

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

        @property
        def location(self):
            return self._location

        @location.setter
        def location(self, entity):
            self._location = ""
            if self.entity.strand == "forward":
                self._location = f"{self.entity.left_end_position}->{self.entity.right_end_position}"
            else:
                self._location = f"{self.entity.left_end_position}<-{self.entity.right_end_position}"

        def to_dict(self):
            dttDatamart = {
                "_id": f"{self.entity.id}_{self.entity.name}" if self.is_fragment else self.entity.id,
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
                           f"Location: {self.location};\n"
                           f"Product: {self.product}"
            }
            return dttDatamart
