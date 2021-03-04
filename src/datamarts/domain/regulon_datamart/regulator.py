import multigenomic_api
from src.datamarts.domain.general.biological_base import BiologicalBase

class Regulator(BiologicalBase):
    def __init__(self, regulator):
        super().__init__(regulator.external_cross_references, regulator.citations, regulator.note)
        self.regulator = regulator
        self.conformations = regulator.active_conformations + regulator.inactive_conformations
        self.genes = regulator.products_ids
        self.operons = regulator.products_ids

    def to_dict(self):
        regulator = {
            "name": self.regulator.name,
            "synonyms": self.regulator.synonyms,
            "note": self.regulator.note,
            "conformations": self.conformations,
            "genes": self.genes,
            "operons": self.operons,
            # TODO Defined on multigenomic/models but not found
            # "sensingClass": self.regulator.sensing_class,
            # "connectivityClass": self.regulator.connectivity_class
        }
        return regulator

    @property
    def conformations(self):
        return self._conformations

    @conformations.setter
    def conformations(self, conformations):
        global product
        self._conformations = []
        for conformation in conformations:
            if conformation.type == "product":
                product = multigenomic_api.products.find_by_id(conformation.id)
            elif conformation.type == "regulatoryComplex":
                product = multigenomic_api.regulatory_complexes.find_by_id(conformation.id)

            super().__init__(product.external_cross_references, product.citations, product.note)
            conformation = {
                "id": product.id,
                "name": product.abbreviated_name or product.name,
                "type": conformation.type,
                "citations": self.citations,
                # TODO: This will be added later
                # "effectorInteractionType": None,
                # TODO: This will be added later
                # "functionalType": None
            }
            self._conformations.append(conformation.copy())

    @property
    def genes(self):
        return self._genes

    @genes.setter
    def genes(self, products_ids):
        self._genes = []
        for product_id in products_ids:
            product = multigenomic_api.products.find_by_id(product_id)
            gene_properties = self.get_gene_properties(multigenomic_api.genes.find_by_id(product.genes_id))
            gene = {
                "id": product.genes_id,
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
            product = multigenomic_api.products.find_by_id(product_id)
            transcription_units = multigenomic_api.transcription_units.find_by_gene_id(product.genes_id)
            operons_id = multigenomic_api.transcription_units.get_operons_id_by_gene_id(product.genes_id)
            # operons_id = list(set(operons_id))
            operon = multigenomic_api.operons.find_by_id(operons_id)
            tus_encoding_reg = []
            global promoter
            for transcription_unit in transcription_units:
                if transcription_unit.promoters_ids:
                    promoter = multigenomic_api.promoters.find_by_id(transcription_unit.promoters_ids[0])
                tu_dict = {
                    "transcriptionUnitName": transcription_unit.name,
                    "promoterName": promoter.name
                }
                if tu_dict not in tus_encoding_reg:
                    tus_encoding_reg.append(tu_dict)
            operon_dict = {
                "operon_id": operon.id,
                "name": operon.name,
                "tusEncodingRegulator": tus_encoding_reg
            }
            self._operons.append(operon_dict.copy())

    @staticmethod
    def get_gene_properties(gene):
        global genome_pos, length
        if gene.fragments:
            for fragment in gene.fragments:
                if fragment.strand == "forward":
                    genome_pos = genome_pos + f"{fragment.left_end_position} -> {fragment.right_end_position};"
                elif fragment.strand == "reverse":
                    genome_pos = genome_pos + f"{fragment.left_end_position} <- {fragment.right_end_position};"
                length = abs(gene.right_end_position - gene.left_end_position) + 1
        else:
            if gene.strand == "forward":
                genome_pos = f"{gene.left_end_position} -> {gene.right_end_position}"
            elif gene.strand == "reverse":
                genome_pos = f"{gene.left_end_position} <- {gene.right_end_position}"
            length = abs(gene.right_end_position - gene.left_end_position) + 1
        return [genome_pos, length]

