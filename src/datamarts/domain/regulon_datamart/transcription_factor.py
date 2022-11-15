import multigenomic_api
from src.datamarts.domain.general.biological_base import BiologicalBase
from src.datamarts.domain.general.additiveEvidences import AdditiveEvidences


class TranscriptionFactor(BiologicalBase):
    def __init__(self, trans_factor):
        super().__init__(trans_factor.external_cross_references, trans_factor.citations, trans_factor.note)
        self.transcription_factor = trans_factor
        self.conformations = trans_factor.active_conformations + trans_factor.inactive_conformations
        self.genes = trans_factor.products_ids
        self.operons = trans_factor.products_ids
        self.products = trans_factor.products_ids

    def to_dict(self):
        citations = self.citations
        additive_evs = AdditiveEvidences(citations)
        transcription_factor = {
            "_id": self.transcription_factor.id,
            "citations": self.citations,
            "name": self.transcription_factor.name,
            "synonyms": self.transcription_factor.synonyms,
            "note": self.formatted_note,
            "conformations": self.conformations,
            "encodedFrom": {
                "genes": self.genes,
                "operon": self.operons
            },
            # TODO: This will be added later by local process
            # "sensingClass": self.regulator.sensing_class,
            # "connectivityClass": self.regulator.connectivity_class,
            "products": self.products,
            "symmetry": self.transcription_factor.symmetry,
            "siteLength": self.transcription_factor.site_length,
            "family": self.transcription_factor.family,
            "additiveEvidences": additive_evs.to_dict(),
            "confidenceLevel": additive_evs.get_confidence_level()
        }
        return transcription_factor

    @property
    def conformations(self):
        return self._conformations

    @conformations.setter
    def conformations(self, conformations):
        product = None
        self._conformations = []
        for conformation in conformations:
            if conformation.type == "product":
                product = multigenomic_api.products.find_by_id(conformation.id)
            elif conformation.type == "regulatoryComplex":
                product = multigenomic_api.regulatory_complexes.find_by_id(conformation.id)
            conformation_object = Conformation(product, conformation.type)
            self._conformations.append(conformation_object.to_dict().copy())

    @property
    def genes(self):
        return self._genes

    @genes.setter
    def genes(self, products_ids):
        self._genes = []
        for product_id in products_ids:
            prod = multigenomic_api.products.find_by_id(product_id)
            gene = multigenomic_api.genes.find_by_id(prod.genes_id)
            gene_properties = self.get_gene_properties(multigenomic_api.genes.find_by_id(prod.genes_id))
            gene = {
                "gene_id": prod.genes_id,
                "gene_name": gene.name,
                "genomePosition": gene_properties[0],
                "length": gene_properties[1]
            }
            self.genes.append(gene.copy())

    @property
    def operons(self):
        return self._operons

    @operons.setter
    def operons(self, products_ids):
        self._operons = []
        for product_id in products_ids:
            prod = multigenomic_api.products.find_by_id(product_id)
            transcription_units = multigenomic_api.transcription_units.find_by_gene_id(prod.genes_id)
            operons_id = multigenomic_api.transcription_units.get_operons_id_by_gene_id(prod.genes_id)
            # operons_id = list(set(operons_id))
            operon = multigenomic_api.operons.find_by_id(operons_id)
            tus_encoding_reg = []
            promoter = None
            for transcription_unit in transcription_units:
                tu_dict = {
                    "transcriptionUnitName": transcription_unit.name
                }
                if transcription_unit.promoters_id:
                    promoter = multigenomic_api.promoters.find_by_id(transcription_unit.promoters_id)
                    tu_dict["promoterName"]: promoter.name
                if tu_dict not in tus_encoding_reg:
                    tus_encoding_reg.append(tu_dict)
            operon_dict = {
                "operon_id": operon.id,
                "name": operon.name,
                "tusEncodingRegulator": tus_encoding_reg
            }
            if operon_dict not in self._operons:
                self._operons.append(operon_dict.copy())

    # TODO: Check if this is correctly obtained
    @property
    def products(self):
        return self._products

    @products.setter
    def products(self, products_ids):
        self._products = []
        for product_id in products_ids:
            prod = multigenomic_api.products.find_by_id(product_id)
            prod_dict = {
                "_id": prod.id,
                "name": prod.name
            }
            if prod_dict not in self._products:
                self._products.append(prod_dict.copy())

    @staticmethod
    def get_gene_properties(gene):
        genome_pos = ""
        length = 0
        if gene.fragments:
            for fragment in gene.fragments:
                if fragment.strand == "forward":
                    genome_pos = genome_pos + f"{fragment.left_end_position} -> {fragment.right_end_position};"
                elif fragment.strand == "reverse":
                    genome_pos = genome_pos + f"{fragment.left_end_position} <- {fragment.right_end_position};"
                if gene.right_end_position and gene.left_end_position:
                    length = abs(gene.right_end_position - gene.left_end_position) + 1
        else:
            if gene.strand == "forward":
                genome_pos = f"{gene.left_end_position} -> {gene.right_end_position}"
            elif gene.strand == "reverse":
                genome_pos = f"{gene.left_end_position} <- {gene.right_end_position}"
            length = abs(gene.right_end_position - gene.left_end_position) + 1
        return [genome_pos, length]


class Conformation(BiologicalBase):
    def __init__(self, product, conformation_type):
        super().__init__([], product.citations, [])
        self.product = product
        self.type = conformation_type

    def to_dict(self):
        citations = self.citations
        additive_evs = AdditiveEvidences(citations)
        conformation = {
            "id": self.product.id,
            "name": self.product.abbreviated_name or self.product.name,
            "type": self.type,
            "citations": citations,
            "additiveEvidences": additive_evs.to_dict(),
            "confidenceLevel": additive_evs.get_confidence_level(),
            # TODO: This will be added later by local process
            # "effectorInteractionType": None,
            # TODO: This will be added later by local process
            # "functionalType": None
        }
        return conformation
