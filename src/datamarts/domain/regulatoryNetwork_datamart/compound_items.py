import multigenomic_api
from src.datamarts.domain.general.biological_base import BiologicalBase


class RegulatoryNetworkCompound:

    @property
    def objects(self):
        compound_objects = get_all_regulators(multigenomic_api.regulatory_interactions.get_all())
        for compound in compound_objects:
            print(compound.id)
            srna_dict = multigenomic_api.regulatory_continuants.find_by_id(compound.id)
            reg_network_node = RegulatoryNetworkCompound.NodeItem(srna_dict)
            yield reg_network_node
        del compound_objects

    class NodeItem(BiologicalBase):

        def __init__(self, node_object):
            super().__init__(node_object.external_cross_references, node_object.citations, node_object.note)
            self.id = node_object.id
            self.node = node_object
            self.outdegree = node_object

        @property
        def outdegree(self):
            return self._outdegree

        @outdegree.setter
        def outdegree(self, node_object):
            self._outdegree = []
            reg_ints = multigenomic_api.regulatory_interactions.find_by_regulator_id(node_object.id)
            for ri in reg_ints:
                genes_ids = []
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
                genes_ids = list(set(genes_ids))
                for gene_id in genes_ids:
                    self._outdegree = outdegree_tf(gene_id, ri.function, node_object.name, self._outdegree)
                    gene_outdegree_item = outdegree_gene(gene_id, ri.function, node_object.name)
                    if gene_outdegree_item not in self._outdegree:
                        self._outdegree.append(gene_outdegree_item.copy())

        def to_dict(self):
            reg_network_node = {
                "_id": self.id,
                "name": self.node.name,
                "type": "Compound",
                "outdegree": self.outdegree
                # "citations": self.citations
            }
            return reg_network_node


def outdegree_gene(gene_id, reg_int_function, object_name):
    gene = multigenomic_api.genes.find_by_id(gene_id)
    tooltip = define_tooltip(reg_int_function, f"Compound {object_name}", f"Gene {gene.name}")
    gene_outdegree_item = BuildDict(gene, "Gene", reg_int_function, tooltip, "Compound-Gene").to_dict()
    return gene_outdegree_item


def outdegree_tf(gene_id, reg_int_function, object_name, outdegree_list):
    trans_factors = []
    products = multigenomic_api.products.find_by_gene_id(gene_id)
    for product in products:
        if product.type == "small RNA":
            tooltip = define_tooltip(reg_int_function, f"Compound {object_name}", f"sRNA {product.abbreviated_name}")
            tf_outdegree_item = BuildDict(product, "sRNA", reg_int_function, tooltip, "Compound-sRNA").to_dict()
            if tf_outdegree_item not in outdegree_list:
                outdegree_list.append(tf_outdegree_item.copy())
        trans_factors.extend(multigenomic_api.transcription_factors.find_tf_id_by_conformation_id(product.id))
        trans_factors.extend(multigenomic_api.transcription_factors.find_tf_id_by_product_id(product.id))
        for tf in trans_factors:
            tooltip = define_tooltip(reg_int_function, f"Compound {object_name}", f"Transcription Factor {tf.abbreviated_name}")
            tf_outdegree_item = BuildDict(tf, "Transcription Factor", reg_int_function, tooltip, "Compound-TF").to_dict()
            if tf_outdegree_item not in outdegree_list:
                outdegree_list.append(tf_outdegree_item)
    return outdegree_list


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
        if ri.regulator.type == "regulatoryContinuant":
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
        if self.item_type == "Gene":
            name = self.item.name
        else:
            name = self.item.abbreviated_name or self.item.name
        item_dict = {
            "_id": self.item.id,
            "name": name,
            "type": self.item_type,
            "regulatoryEffect": self.reg_int_function[0] or "unknown",
            "citations": self.citations,
            "tooltip": self.tooltip,
            "networkType": self.network_type
        }
        return item_dict