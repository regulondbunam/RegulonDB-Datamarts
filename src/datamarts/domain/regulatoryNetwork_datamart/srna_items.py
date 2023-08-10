import multigenomic_api
from src.datamarts.domain.general.biological_base import BiologicalBase


class RegulatoryNetworkSRNA:

    @property
    def objects(self):
        srna_objects = get_all_regulators(multigenomic_api.regulatory_interactions.get_all())
        for srna in srna_objects:
            print(srna.id)
            srna_dict = multigenomic_api.products.find_by_id(srna.id)
            reg_network_node = RegulatoryNetworkSRNA.NodeItem(srna_dict)
            yield reg_network_node
        del srna_objects

    class NodeItem(BiologicalBase):

        def __init__(self, node_object):
            super().__init__(node_object.external_cross_references, node_object.citations, node_object.note)
            self.id = node_object.id
            self.node = node_object
            self.outdegree = node_object
            self.indegree = node_object

        @property
        def outdegree(self):
            return self._outdegree

        @outdegree.setter
        def outdegree(self, node_object):
            self._outdegree = []
            genes_ids = []
            reg_ints = multigenomic_api.regulatory_interactions.find_by_regulator_id(node_object.id)
            for ri in reg_ints:
                if ri.regulated_entity.type == "promoter":
                    tus = multigenomic_api.transcription_units.find_by_promoter_id(ri.regulated_entity.id)
                    for tu in tus:
                        for gene_id in tu.genes_ids:
                            genes_ids.append(gene_id)
                elif ri.regulated_entity.type == "transcriptionUnit":
                    tu = multigenomic_api.transcription_units.find_by_id(ri.regulated_entity.id)
                    for gene_id in tu.genes_ids:
                        genes_ids.append(gene_id)
                elif ri.regulated_entity.type == "gene":
                    genes_ids.append(ri.regulated_entity.id)
                for gene_id in genes_ids:
                    gene_outdegree_item = outdegree_gene(gene_id, ri.function, node_object.name)
                    if gene_outdegree_item not in self._outdegree:
                        self._outdegree.append(gene_outdegree_item.copy())

        @property
        def indegree(self):
            return self._indegree

        @indegree.setter
        def indegree(self, node_object):
            self._indegree = []
            product = multigenomic_api.products.find_by_id(node_object.id)
            reg_ints = multigenomic_api.regulatory_interactions.find_regulatory_interactions_by_reg_entity_id(
                product.genes_id)
            if len(reg_ints) > 0:
                self._indegree = indegree_gene(reg_ints, self._indegree, node_object)
            trans_units = multigenomic_api.transcription_units.find_by_gene_id(product.genes_id)
            if len(trans_units):
                for tu in trans_units:
                    reg_ints = multigenomic_api.regulatory_interactions. \
                        find_regulatory_interactions_by_reg_entity_id(tu.id)
                    if len(reg_ints) > 0:
                        self._indegree = indegree_gene(reg_ints, self._indegree, node_object)
                    reg_ints = multigenomic_api.regulatory_interactions. \
                        find_regulatory_interactions_by_reg_entity_id(tu.promoters_id)
                    if len(reg_ints) > 0:
                        self._indegree = indegree_gene(reg_ints, self._indegree, node_object)

        def to_dict(self):
            reg_network_node = {
                "_id": self.id,
                "name": self.node.name,
                "type": "Small RNA",
                "indegree": self.indegree,
                "outdegree": self.outdegree,
                # "citations": self.citations
            }
            return reg_network_node


def outdegree_gene(gene_id, reg_int_function, object_name):
    gene = multigenomic_api.genes.find_by_id(gene_id)
    tooltip = define_tooltip(reg_int_function, f"sRNA {object_name}", f"Gene {gene.name}")
    gene_outdegree_item = BuildDict(gene, "Gene", reg_int_function, tooltip, "sRNA-Gene").to_dict()
    return gene_outdegree_item


def indegree_gene(reg_ints, indegree_list, node_object):
    for ri in reg_ints:
        if ri.regulator:
            products = []
            if ri.regulator.type == "product":
                product = multigenomic_api.products.find_by_id(ri.regulator.id)
                products.append(product)
            elif ri.regulator.type == "regulatoryComplex":
                reg_complex = multigenomic_api.regulatory_complexes.find_by_id(ri.regulator.id)
                if reg_complex.products:
                    for product in reg_complex.products:
                        products.append(multigenomic_api.products.find_by_id(product.products_id))
            for product in products:
                gene = multigenomic_api.genes.find_by_id(product.genes_id)
                tooltip = define_tooltip(ri.function, f"Gene {gene.name}",
                                         f"sRNA {node_object.name}")
                gene_indegree_item = BuildDict(gene, "Gene", ri.function, tooltip, "Gene-sRNA").to_dict()
                if gene_indegree_item not in indegree_list:
                    indegree_list.append(gene_indegree_item.copy())
    return indegree_list


def define_tooltip(function, regulator, regulated):
    if function == "repressor":
        tooltip = f"{regulator} represses {regulated}"
    elif function == "activator":
        tooltip = f"{regulator} activates {regulated}"
    else:
        tooltip = f"{regulator} function to {regulated} is unknown"
    return tooltip


def get_all_regulators(reg_ints):
    regulators = []
    for ri in reg_ints:
        if ri.regulator.type == "product":
            product = multigenomic_api.products.find_by_id(ri.regulator.id)
            if product.type == "small RNA":
                if ri.regulator not in regulators:
                    regulators.append(ri.regulator)
    return regulators


class BuildDict(BiologicalBase):
    def __init__(self, item, item_type, reg_int_function, tooltip, network_type):
        super().__init__([], item.citations, [])
        self.item = item
        self.item_type = item_type
        self.reg_int_function = reg_int_function,
        self.tooltip = tooltip
        self.network_type = network_type

    def to_dict(self):
        item_dict = {
            "_id": self.item.id,
            "name": self.item.name,
            "type": self.item_type,
            "regulatoryEffect": self.reg_int_function[0] or "unknown",
            "citations": self.citations,
            "tooltip": self.tooltip,
            "networkType": self.network_type
        }
        return item_dict
