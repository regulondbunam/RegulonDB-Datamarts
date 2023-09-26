import multigenomic_api
from mongoengine.errors import DoesNotExist

from src.datamarts.domain.general.biological_base import BiologicalBase


class Regulates(BiologicalBase):
    def __init__(self, regulator):
        super().__init__(regulator.external_cross_references, regulator.citations, regulator.note)
        trans_units = get_all_transcription_units(regulator)
        regulated_genes = get_all_regulated_genes(regulator)
        self.genes = regulated_genes
        self.transcription_factors = trans_units
        self.transcription_units = trans_units
        self.operons = trans_units
        self.sigma_factors = trans_units

    def to_dict(self):
        regulates = {
            "genes": self.genes,
            "transcriptionFactors": self.transcription_factors,
            "transcriptionUnits": self.transcription_units,
            "operons": self.operons,
            "sigmaFactors": self.sigma_factors
        }
        return regulates

    @property
    def genes(self):
        return self._genes

    @genes.setter
    def genes(self, reg_genes):
        self._genes = []
        for tu_object in reg_genes:
            genes_ids = tu_object["genes"]
            for gene_id in genes_ids:
                gene = multigenomic_api.genes.find_by_id(gene_id)
                terms = get_gene_terms(gene)
                products = multigenomic_api.products.find_by_gene_id(gene_id)
                gene_object = {
                    "_id": gene.id,
                    "name": gene.name,
                    "terms": {
                        "multifun": terms,
                        "geneOntology": gene_ontology_extrac(products)
                    },
                    "function": tu_object["function"]
                }
                self._genes = insert_element_on(gene_object, self._genes)

    @property
    def transcription_factors(self):
        return self._transcription_factors

    @transcription_factors.setter
    def transcription_factors(self, trans_units):
        self._transcription_factors = []
        for tu_object in trans_units:
            transcription_units = tu_object["transcription_units"]
            if transcription_units:
                for tu in transcription_units:
                    for gene_id in tu.genes_ids:
                        products = multigenomic_api.products.find_by_gene_id(gene_id)
                        for product in products:
                            trans_factors = multigenomic_api.transcription_factors.\
                                find_tf_id_by_conformation_id(product.id)
                            for tf in trans_factors:
                                genes = []
                                for product_id in tf.products_ids:
                                    prod = multigenomic_api.products.find_by_id(product_id)
                                    gene = multigenomic_api.genes.find_by_id(prod.genes_id)
                                    terms = get_gene_terms(gene)
                                    products = multigenomic_api.products.find_by_gene_id(gene_id)
                                    gene = {
                                        "_id": prod.genes_id,
                                        "name": gene.name,
                                        "terms": {
                                            "multifun": terms,
                                            "geneOntology": gene_ontology_extrac(products)
                                        }
                                    }
                                    if gene.copy() not in genes:
                                        genes.append(gene.copy())
                                tf_object = {
                                    "_id": tf.id,
                                    "name": tf.abbreviated_name,
                                    "function": tu_object["function"],
                                    "genes": genes
                                }
                                self._transcription_factors = insert_element_on(tf_object, self._transcription_factors)

    @property
    def transcription_units(self):
        return self._transcription_units

    @transcription_units.setter
    def transcription_units(self, trans_units):
        self._transcription_units = []
        for tu_object in trans_units:
            transcription_units = tu_object["transcription_units"]
            if transcription_units:
                for tu in transcription_units:
                    if tu.promoters_id:
                        promoter = multigenomic_api.promoters.find_by_id(tu.promoters_id)
                        first_gene = get_first_gene_of_tu(tu.genes_ids, promoter)
                        tu_object = {
                            "_id": tu.id,
                            "name": tu.name,
                            "firstGene": first_gene,
                            "function": tu_object["function"],
                            "promoter": {
                                "_id": promoter.id,
                                "name": promoter.name
                            }
                        }
                        self._transcription_units = insert_element_on(tu_object, self._transcription_units)

    @property
    def operons(self):
        return self._operons

    @operons.setter
    def operons(self, trans_units):
        self._operons = []
        for tu_object in trans_units:
            transcription_units = tu_object["transcription_units"]
            if transcription_units:
                genes_ids = []
                for tu in transcription_units:
                    genes_ids.extend(tu.genes_ids)
                for tu in transcription_units:
                    if tu.promoters_id:
                        promoter = multigenomic_api.promoters.find_by_id(tu.promoters_id)
                        first_gene = get_first_gene_of_tu(genes_ids, promoter)
                        operon = multigenomic_api.operons.find_by_id(tu.operons_id)
                        operon_object = {
                            "_id": operon.id,
                            "name": operon.name,
                            "firstGene": first_gene,
                            "function": tu_object["function"]
                        }
                        self._operons = insert_element_on(operon_object, self._operons)

    @property
    def sigma_factors(self):
        return self._sigma_factors

    @sigma_factors.setter
    def sigma_factors(self, trans_units):
        self._sigma_factors = []
        for tu_object in trans_units:
            transcription_units = tu_object["transcription_units"]
            if transcription_units:
                for tu in transcription_units:
                    if tu.promoters_id:
                        promoter = multigenomic_api.promoters.find_by_id(tu.promoters_id)
                        if promoter.binds_sigma_factor:
                            sigma_factor = multigenomic_api.\
                                sigma_factors.find_by_id(promoter.binds_sigma_factor.sigma_factors_id)
                            gene = multigenomic_api.genes.find_by_id(sigma_factor.genes_id)
                            sigma_factor_object = {
                                "_id": sigma_factor.id,
                                "name": sigma_factor.abbreviated_name,
                                "function": tu_object["function"],
                                "gene": {
                                    "_id": gene.id,
                                    "name": gene.name
                                }
                            }
                            if sigma_factor_object not in self._sigma_factors:
                                self._sigma_factors = insert_element_on(sigma_factor_object, self._sigma_factors)


def get_all_transcription_units(regulator):
    all_transcription_units = []
    conformations_ids = []
    if regulator.regulator_type == "transcriptionFactor":
        reg_complex = None
        try:
            reg_complex = multigenomic_api.regulatory_complexes.find_by_name(regulator.name)
        except DoesNotExist:
            pass
        if reg_complex is not None:
            conformations_ids.append(reg_complex.id)
        for active_conf in regulator.active_conformations:
            conformations_ids.append(active_conf.id)
        conformations_ids.extend(regulator.products_ids)
    else:
        conformations_ids.append(regulator.id)
    for conformation_id in conformations_ids:
        regulatory_interactions = multigenomic_api.regulatory_interactions.find_by_regulator_id(conformation_id)
        for ri in regulatory_interactions:
            transcription_units = []
            if ri.regulated_entity.type == "promoter":
                transcription_units = multigenomic_api.transcription_units.find_by_promoter_id(
                    ri.regulated_entity.id)
            elif ri.regulated_entity.type == "transcriptionUnit":
                tu = multigenomic_api.transcription_units.find_by_id(ri.regulated_entity.id)
                transcription_units = [tu]
            tu_object = {
                "transcription_units": transcription_units,
                "function": ri.function
            }
            if tu_object not in all_transcription_units:
                all_transcription_units.append(tu_object)
    return all_transcription_units


def get_all_regulated_genes(regulator):
    conformations_ids = []
    all_reg_genes = []
    if regulator.regulator_type == "transcriptionFactor":
        reg_complex = None
        try:
            reg_complex = multigenomic_api.regulatory_complexes.find_by_name(regulator.name)
        except DoesNotExist:
            pass
        if reg_complex is not None:
            conformations_ids.append(reg_complex.id)
        for active_conf in regulator.active_conformations:
            conformations_ids.append(active_conf.id)
        conformations_ids.extend(regulator.products_ids)
    else:
        conformations_ids.append(regulator.id)
    for conformation_id in conformations_ids:
        regulatory_interactions = multigenomic_api.regulatory_interactions.find_by_regulator_id(conformation_id)
        for ri in regulatory_interactions:
            regulated_genes = []
            if ri.regulated_entity.type == "promoter":
                transcription_units = multigenomic_api.transcription_units.find_by_promoter_id(
                    ri.regulated_entity.id)
                for tu in transcription_units:
                    regulated_genes.extend(tu.genes_ids)
            elif ri.regulated_entity.type == "transcriptionUnit":
                tu = multigenomic_api.transcription_units.find_by_id(ri.regulated_entity.id)
                regulated_genes.extend(tu.genes_ids)
            if ri.regulated_entity.type == "gene":
                regulated_genes.append(ri.regulated_entity.id)
            gene_object = {
                "genes": list(set(regulated_genes)),
                "function": ri.function
            }
            if gene_object not in all_reg_genes:
                all_reg_genes.append(gene_object)
    return all_reg_genes


def gene_ontology_extrac(products):
    terms_dict = {
        'biologicalProcess': [],
        'cellularComponent': [],
        'molecularFunction': []
    }
    for product in products:
        terms = product.terms
        if terms:
            for term in terms.biological_process:
                term = Term(term)
                term_dict = term.to_dict()
                if term_dict not in terms_dict['biologicalProcess']:
                    terms_dict['biologicalProcess'].append(term_dict)

            for term in terms.cellular_component:
                term = Term(term)
                term_dict = term.to_dict()
                if term_dict not in terms_dict['cellularComponent']:
                    terms_dict['cellularComponent'].append(term_dict)

            for term in terms.molecular_function:
                term = Term(term)
                term_dict = term.to_dict()
                if term_dict not in terms_dict['molecularFunction']:
                    terms_dict['molecularFunction'].append(term_dict)
    return terms_dict


class Term(BiologicalBase):

    def __init__(self, term):
        super().__init__([], [], [])
        self.term = term

    def to_dict(self):
        term = {
            '_id': self.term.terms_id,
            'name': self.term.terms_name
        }
        return term


def get_first_gene_of_tu(genes_ids, promoter):
    first_gene = multigenomic_api.genes.find_by_id(genes_ids[0])
    if promoter.strand == "reverse":
        first_gene.left_end_position = first_gene.right_end_position
    first_gene.left_end_position = first_gene.left_end_position or first_gene.fragments[0].left_end_position

    for gene in genes_ids:
        current_gene = multigenomic_api.genes.find_by_id(gene)
        if promoter.strand == "forward":
            if current_gene.left_end_position:
                if current_gene.left_end_position < first_gene.left_end_position:
                    first_gene = current_gene
            elif current_gene.fragments:
                for fragment in current_gene.fragments:
                    if fragment.left_end_position < first_gene.left_end_position:
                        first_gene = current_gene
                        first_gene.left_end_position = fragment.left_end_position

        elif promoter.strand == "reverse":
            if current_gene.left_end_position:
                if current_gene.right_end_position > first_gene.left_end_position:
                    first_gene = current_gene
                    first_gene.left_end_position = first_gene.right_end_position
            elif current_gene.fragments:
                for fragment in current_gene.fragments:
                    if fragment.right_end_position > first_gene.left_end_position:
                        first_gene = current_gene
                        first_gene.left_end_position = fragment.right_end_position

        first_gene.left_end_position = first_gene.left_end_position or first_gene.fragments[0].left_end_position
    first_gene = {
        "_id": first_gene.id,
        "name": first_gene.name
    }
    return first_gene


def insert_element_on(reg_object, regulated_list):
    for regulated in regulated_list:
        if reg_object["_id"] == regulated["_id"]:
            if reg_object['function'] != regulated['function']:
                reg_object["function"] = "dual"
                regulated["function"] = "dual"
    if reg_object not in regulated_list:
        regulated_list.append(reg_object)
    return regulated_list


def get_gene_terms(gene):
    terms = gene.terms
    terms_list = []
    for term in terms:
        term = {
            '_id': term.terms_id,
            'name': term.terms_name
        }
        terms_list.append(term.copy())
    return terms_list

