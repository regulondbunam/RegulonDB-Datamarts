import multigenomic_api
from mongoengine.errors import DoesNotExist

from src.datamarts.domain.general.biological_base import BiologicalBase
from src.datamarts.domain.general.additiveEvidences import AdditiveEvidences


class Regulator(BiologicalBase):
    def __init__(self, regulator):
        super().__init__(regulator.external_cross_references, regulator.citations, regulator.note)
        self.regulator = regulator
        self.conformations = regulator
        self.genes = regulator
        self.operons = regulator
        self.products = regulator

    def to_dict(self):
        citations = self.citations
        additive_evs = AdditiveEvidences(citations)
        note = self.formatted_note
        if not note:
            note = self.get_longest_note()
        if self.regulator.regulator_type == "transcriptionFactor":
            regulator = {
                "_id": self.regulator.id,
                "citations": self.citations,
                "abbreviatedName": self.regulator.abbreviated_name,
                "name": self.regulator.name,
                "synonyms": self.regulator.synonyms,
                "note": note,
                "conformations": self.conformations,
                "encodedBy": {
                    "genes": self.genes,
                    "operon": self.operons
                },
                # TODO: This will be added later by local process
                # "sensingClass": self.regulator.sensing_class,
                # "connectivityClass": self.regulator.connectivity_class,
                "products": self.products,
                "symmetry": self.regulator.symmetry,
                "siteLength": self.regulator.site_length,
                "family": self.regulator.family,
                "additiveEvidences": additive_evs.to_dict(),
                "confidenceLevel": additive_evs.get_confidence_level(),
                "type": self.regulator.regulator_type
            }
        else:
            abbreviated_name = None
            if self.regulator.regulator_type == "sRNA":
                abbreviated_name = self.regulator.abbreviated_name
            regulator = {
                "_id": self.regulator.id,
                "citations": self.citations,
                "name": self.regulator.name,
                "encodedBy": {
                    "genes": self.genes
                },
                "abbreviatedName": abbreviated_name,
                "synonyms": self.regulator.synonyms,
                "note": self.formatted_note,
                "additiveEvidences": additive_evs.to_dict(),
                "confidenceLevel": additive_evs.get_confidence_level(),
                "type": self.regulator.regulator_type
            }
        return regulator

    @property
    def conformations(self):
        return self._conformations

    @conformations.setter
    def conformations(self, regulator):
        self._conformations = []
        reg_complex = None
        product = None
        if regulator.regulator_type == "transcriptionFactor":
            try:
                reg_complex = multigenomic_api.regulatory_complexes.find_by_name(regulator.name)
            except DoesNotExist:
                pass
            if reg_complex is not None:
                reg_ints = multigenomic_api.regulatory_interactions.find_by_regulator_id(
                    reg_complex.id)
                if len(reg_ints) > 0:
                    self._conformations.append(Conformation(reg_complex, "regulatoryComplex", "active").to_dict().copy())
            try:
                product = multigenomic_api.products.find_by_name(regulator.name)
            except DoesNotExist:
                pass
            if product is not None:
                reg_ints = multigenomic_api.regulatory_interactions.find_by_regulator_id(product.id)
                if len(reg_ints) > 0:
                    self._conformations.append(Conformation(product, "product", "active").to_dict().copy())
            self._conformations = insert_conf(regulator.active_conformations, self._conformations, "active")
            self._conformations = insert_conf(regulator.inactive_conformations, self._conformations, "inactive")

    @property
    def genes(self):
        return self._genes

    @genes.setter
    def genes(self, regulator):
        self._genes = []
        if regulator.regulator_type == "transcriptionFactor":
            for product_id in regulator.products_ids:
                prod = multigenomic_api.products.find_by_id(product_id)
                gene = self.get_gene_properties(prod.genes_id)
                self.genes.append(gene.copy())
        elif regulator.regulator_type == "sRNA":
            prod = multigenomic_api.products.find_by_id(regulator.id)
            gene = self.get_gene_properties(prod.genes_id)
            self.genes.append(gene.copy())

    @property
    def operons(self):
        return self._operons

    @operons.setter
    def operons(self, regulator):
        self._operons = []
        if regulator.regulator_type == "transcriptionFactor":
            for product_id in regulator.products_ids:
                prod = multigenomic_api.products.find_by_id(product_id)
                transcription_units = multigenomic_api.transcription_units.find_by_gene_id(prod.genes_id)
                operons_id = multigenomic_api.transcription_units.get_operons_id_by_gene_id(prod.genes_id)
                operon = multigenomic_api.operons.find_by_id(operons_id)
                tus_encoding_reg = []
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
                    "_id": operon.id,
                    "name": operon.name,
                    "tusEncodingRegulator": tus_encoding_reg
                }
                if operon_dict not in self._operons:
                    self._operons.append(operon_dict.copy())

    @property
    def products(self):
        return self._products

    @products.setter
    def products(self, regulator):
        self._products = []
        if regulator.regulator_type == "transcriptionFactor":
            for product_id in regulator.products_ids:
                prod = multigenomic_api.products.find_by_id(product_id)
                prod_dict = {
                    "_id": prod.id,
                    "name": prod.name,
                    "abbreviatedName": prod.abbreviated_name
                }
                if prod_dict not in self._products:
                    self._products.append(prod_dict.copy())

    @staticmethod
    def get_gene_properties(gene_id):
        gene_obj = multigenomic_api.genes.find_by_id(gene_id)
        length = None
        if gene_obj.left_end_position and gene_obj.right_end_position:
            length = abs(gene_obj.right_end_position - gene_obj.left_end_position) + 1
        gene = {
            "_id": gene_id,
            "name": gene_obj.name,
            "leftEndPosition": gene_obj.left_end_position,
            "rightEndPosition": gene_obj.right_end_position,
            "length": length
        }
        return gene

    def get_longest_note(self):
        longest_note = ""
        for conf in self.conformations:
            if conf["note"]:
                if len(conf["note"]) > len(longest_note):
                    longest_note = conf["note"]
        return longest_note


def insert_conf(conformations, final_list, conf_class):
    conformation_object = None
    for conformation in conformations:
        if conformation.type == "product":
            prod = multigenomic_api.products.find_by_id(conformation.id)
            conformation_object = Conformation(prod, conformation.type, conf_class)
        elif conformation.type == "regulatoryComplex":
            complx = multigenomic_api.regulatory_complexes.find_by_id(conformation.id)
            conformation_object = Conformation(complx, conformation.type, conf_class)
        if conformation_object.to_dict() not in final_list:
            final_list.append(conformation_object.to_dict().copy())
    return final_list


class Conformation(BiologicalBase):
    def __init__(self, conf, conformation_type, conf_class):
        super().__init__([], conf.citations, conf.note)
        self.conf = conf
        self.type = conformation_type
        self.name = None
        self.effector = conf
        self.conf_class = conf_class

    @property
    def effector(self):
        return self._effector

    @effector.setter
    def effector(self, regulator):
        self._effector = []
        if self.type == "regulatoryComplex":
            for continuant_id in regulator.regulatory_continuants_ids:
                continuant = multigenomic_api.regulatory_continuants.find_by_id(continuant_id)
                effector = {
                    "_id": continuant.id,
                    "name": continuant.name
                }
                self._effector.append(effector)

    def to_dict(self):
        citations = self.citations
        additive_evs = AdditiveEvidences(citations)
        if self.type != "regulatoryContinuant":
            self.name = self.conf.abbreviated_name or self.conf.name
        else:
            self.name = self.conf.name
        conformation = {
            "_id": self.conf.id,
            "name": self.name,
            "type": self.type,
            "effector": self.effector,
            "note": self.formatted_note,
            "citations": self.citations,
            "additiveEvidences": additive_evs.to_dict(),
            "confidenceLevel": additive_evs.get_confidence_level(),
            "class": self.conf_class,
            # TODO: This will be added later by local process
            # "effectorInteractionType": None,
            # TODO: This will be added later by local process
            # "functionalType": None
        }
        return conformation

