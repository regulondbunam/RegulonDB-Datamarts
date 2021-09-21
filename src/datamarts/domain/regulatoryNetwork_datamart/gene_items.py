import multigenomic_api


class RegulatoryNetworkGene:

    @property
    def objects(self):
        genes = multigenomic_api.genes.get_all()
        for gene in genes:
            print(gene.id)
            reg_network_node = RegulatoryNetworkGene.NodeItem(gene)
            yield reg_network_node
        del genes

    class NodeItem:

        def __init__(self, node_object):
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
            self._outdegree = outdegree_genes(node_object, self._outdegree)
            self._outdegree = outdegree_tfs(node_object, self._outdegree)

        @property
        def indegree(self):
            return self._indegree

        @indegree.setter
        def indegree(self, node_object):
            self._indegree = []
            self._indegree = indegree_items(node_object, self._indegree)

        def to_dict(self):
            reg_network_node = {
                "_id": self.id,
                "name": self.node.name,
                "type": "Gene",
                "indegree": self.indegree,
                "outdegree": self.outdegree,
            }
            return reg_network_node


def outdegree_genes(node_object, outdegree):
    products = multigenomic_api.products.find_by_gene_id(node_object.id)
    for conformation in products:
        reg_ints = multigenomic_api.regulatory_interactions.find_by_regulator_id(conformation.id)
        for ri in reg_ints:
            if ri.regulated_entity.type == "promoter":
                tus = multigenomic_api.transcription_units.find_by_promoter_id(ri.regulated_entity.id)
                for tu in tus:
                    for gene_id in tu.genes_ids:
                        gene_outdegree_item = build_outdegree_gene_dict(gene_id, ri.function, node_object.name)
                        if gene_outdegree_item not in outdegree:
                            outdegree.append(gene_outdegree_item.copy())
            elif ri.regulated_entity.type == "transcriptionUnit":
                tu = multigenomic_api.transcription_units.find_by_id(ri.regulated_entity.id)
                for gene_id in tu.genes_ids:
                    gene_outdegree_item = build_outdegree_gene_dict(gene_id, ri.function, node_object.name)
                    if gene_outdegree_item not in outdegree:
                        outdegree.append(gene_outdegree_item.copy())
            elif ri.regulated_entity.type == "gene":
                gene_outdegree_item = build_outdegree_gene_dict(ri.regulated_entity.id, ri.function, node_object.name)
                if gene_outdegree_item not in outdegree:
                    outdegree.append(gene_outdegree_item.copy())
    return outdegree


def outdegree_tfs(node_object, outdegree):
    products = multigenomic_api.products.find_by_gene_id(node_object.id)
    for conformation in products:
        reg_ints = multigenomic_api.regulatory_interactions.find_by_regulator_id(conformation.id)
        for ri in reg_ints:
            if ri.regulated_entity.type == "promoter":
                tus = multigenomic_api.transcription_units.find_by_promoter_id(ri.regulated_entity.id)
                for tu in tus:
                    for gene_id in tu.genes_ids:
                        outdegree = build_outdegree_tf_dict(gene_id, ri.function, node_object.name, outdegree)
            elif ri.regulated_entity.type == "transcriptionUnit":
                tu = multigenomic_api.transcription_units.find_by_id(ri.regulated_entity.id)
                for gene_id in tu.genes_ids:
                    outdegree = build_outdegree_tf_dict(gene_id, ri.function, node_object.name, outdegree)
            elif ri.regulated_entity.type == "gene":
                outdegree = build_outdegree_tf_dict(ri.regulated_entity.id, ri.function, node_object.name, outdegree)
    return outdegree


def build_outdegree_gene_dict(gene_id, reg_int_function, object_name):
    tooltip = ""
    gene = multigenomic_api.genes.find_by_id(gene_id)
    if reg_int_function == "repressor":
        tooltip = f"Gene {object_name} represses gene {gene.name}"
    elif reg_int_function == "activator":
        tooltip = f"Gene {object_name} activates gene {gene.name}"
    else:
        tooltip = f"Gene {object_name} function to gene {gene.name} is unknown"
    gene_outdegree_item = {
        "_id": gene.id,
        "name": gene.name,
        "type": "Gene",
        "regulatoryEffect": reg_int_function or "unknown",
        "citations": [],
        "tooltip": tooltip,
        "networkType": "Gene-Gene"
    }
    return gene_outdegree_item


def build_outdegree_tf_dict(gene_id, reg_int_function, object_name, outdegree_lits):
    tooltip = ""
    products = multigenomic_api.products.find_by_gene_id(gene_id)
    for product in products:
        trans_factors = multigenomic_api.transcription_factors.find_tf_id_by_active_conformation_id(product.id)
        for tf in trans_factors:
            if reg_int_function == "repressor":
                tooltip = f"Gene {object_name} represses Transcription factor {tf.name}"
            elif reg_int_function == "activator":
                tooltip = f"Gene {object_name} activates Transcription factor {tf.name}"
            else:
                tooltip = f"Gene {object_name} function to Transcription factor {tf.name} is unknown"
            tf_outdegree_item = {
                "_id": tf.id,
                "name": tf.name,
                "type": "Transcription Factor",
                "regulatoryEffect": reg_int_function or "unknown",
                "citations": [],
                "tooltip": tooltip,
                "networkType": "Gene-TF"
            }
            if tf_outdegree_item not in outdegree_lits:
                outdegree_lits.append(tf_outdegree_item)
    return outdegree_lits


def build_tf_dict(gene_id, reg_int_function, object_name, outdegree_lits):
    tooltip = ""
    gene = multigenomic_api.genes.find_by_id(gene_id)
    if reg_int_function == "repressor":
        tooltip = f"Gene {object_name} represses Gene {gene.name}"
    elif reg_int_function == "activator":
        tooltip = f"Gene {object_name} activates Gene {gene.name}"
    else:
        tooltip = f"Gene {object_name} function to Gene {gene.name} is unknown"
    tf_outdegree_item = {
        "_id": gene.id,
        "name": gene.name,
        "type": "Gene",
        "regulatoryEffect": reg_int_function or "unknown",
        "citations": [],
        "tooltip": tooltip
    }
    if tf_outdegree_item not in outdegree_lits:
        outdegree_lits.append(tf_outdegree_item)
    return outdegree_lits


def indegree_items(node_object, indegree):
    products = multigenomic_api.products.find_by_gene_id(node_object.id)
    if products:
        for product in products:
            reg_ints = multigenomic_api.regulatory_interactions.find_regulatory_interactions_by_reg_entity_id(
                product.genes_id)
            if len(reg_ints) > 0:
                indegree = build_indegree_tf_dict(reg_ints, indegree, node_object)
                indegree = build_indegree_gene_dict(reg_ints, indegree, node_object)
            trans_units = multigenomic_api.transcription_units.find_by_gene_id(product.genes_id)
            if len(trans_units):
                for tu in trans_units:
                    reg_ints = multigenomic_api.regulatory_interactions.find_regulatory_interactions_by_reg_entity_id(
                        tu.id)
                    if len(reg_ints) > 0:
                        indegree = build_indegree_tf_dict(reg_ints, indegree, node_object)
                        indegree = build_indegree_gene_dict(reg_ints, indegree, node_object)
                    reg_ints = multigenomic_api.regulatory_interactions.find_regulatory_interactions_by_reg_entity_id(
                        tu.promoters_id)
                    if len(reg_ints) > 0:
                        indegree = build_indegree_tf_dict(reg_ints, indegree, node_object)
                        indegree = build_indegree_gene_dict(reg_ints, indegree, node_object)
    return indegree


def remove_none_fields_empty_lists(gene_object):
    if isinstance(gene_object, dict):
        return {property: remove_none_fields_empty_lists(property_value) for property, property_value in gene_object.items() if property_value}
    elif isinstance(gene_object, list):
        if len(gene_object) != 0:
            return [remove_none_fields_empty_lists(v) for v in gene_object]
    else:
        return gene_object


def build_indegree_gene_dict(reg_ints, indegree_list, node_object):
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
                tooltip = ""
                if ri.function == "repressor":
                    tooltip = f"Gene {gene.name} represses Gene {node_object.name}"
                elif ri.function == "activator":
                    tooltip = f"Gene {gene.name} activates Gene {node_object.name}"
                else:
                    tooltip = f"Gene {gene.name} function to Gene {node_object.name} is unknown"
                gene_indegree_item = {
                    "_id": gene.id,
                    "name": gene.name,
                    "type": "Gene",
                    "regulatoryEffect": ri.function or "unknown",
                    "citations": [],
                    "tooltip": tooltip,
                    "networkType": "Gene-Gene"
                }
                if gene_indegree_item not in indegree_list:
                    indegree_list.append(gene_indegree_item.copy())
    return indegree_list


def build_indegree_tf_dict(reg_ints, indegree_list, node_object):
    for ri in reg_ints:
        if ri.regulator:
            trans_factors = multigenomic_api.transcription_factors.find_tf_id_by_active_conformation_id(ri.regulator.id)
            for tf in trans_factors:
                tooltip = ""
                if ri.function == "repressor":
                    tooltip = f"Transcription factor {tf.name} represses Gene {node_object.name}"
                elif ri.function == "activator":
                    tooltip = f"Transcription factor {tf.name} activates Gene {node_object.name}"
                else:
                    tooltip = f"Transcription factor {tf.name} function to Gene {node_object.name} is unknown"
                gene_indegree_item = {
                    "_id": tf.id,
                    "name": tf.name,
                    "type": "Transcription Factor",
                    "regulatoryEffect": ri.function or "unknown",
                    "citations": [],
                    "tooltip": tooltip,
                    "networkType": "TF-Gene"
                }
                if gene_indegree_item not in indegree_list:
                    indegree_list.append(gene_indegree_item.copy())
    return indegree_list

