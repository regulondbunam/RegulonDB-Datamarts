from src.datamarts.domain.regulatoryNetwork_datamart.tf_items import RegulatoryNetworkTF
from src.datamarts.domain.regulatoryNetwork_datamart.gene_items import RegulatoryNetworkGene
from src.datamarts.domain.general.biological_base import BiologicalBase


def all_regulatory_network_nodes():
    network_nodes_tf = RegulatoryNetworkTF()
    network_nodes_gene = RegulatoryNetworkGene()
    json_nodes = []
    for item in network_nodes_tf.objects:
        item = locate_dual_items(item)
        item = remove_none_fields_empty_lists(item.to_dict())
        json_nodes.append(remove_none_fields_empty_lists(item))
    for item in network_nodes_gene.objects:
        item = locate_dual_items(item)
        item = remove_none_fields_empty_lists(item.to_dict())
        json_nodes.append(remove_none_fields_empty_lists(item))
    return json_nodes


def remove_none_fields_empty_lists(object_item):
    if isinstance(object_item, dict):
        return {property: remove_none_fields_empty_lists(property_value) for property, property_value in
                object_item.items() if property_value}
    elif isinstance(object_item, list):
        if len(object_item) != 0:
            return [remove_none_fields_empty_lists(v) for v in object_item]
    else:
        return object_item


def locate_dual_items(item):
    if item.outdegree:
        for first_node in item.outdegree:
            for node in item.outdegree:
                if first_node["_id"] == node["_id"]:
                    if first_node["regulatoryEffect"] == "activator" and node["regulatoryEffect"] == "repressor":
                        node["regulatoryEffect"] = "dual"
                        node["tooltip"] = node["tooltip"].replace("represses", "activates and represses")
                        item.outdegree.remove(first_node)
                    elif first_node["regulatoryEffect"] == "repressor" and node["regulatoryEffect"] == "activator":
                        node["regulatoryEffect"] = "dual"
                        node["tooltip"] = node["tooltip"].replace("activates", "activates and represses")
                        item.outdegree.remove(first_node)
    if item.indegree:
        for first_node in item.indegree:
            for node in item.indegree:
                if first_node["_id"] == node["_id"]:
                    if first_node["regulatoryEffect"] == "activator" and node["regulatoryEffect"] == "repressor":
                        node["regulatoryEffect"] = "dual"
                        node["tooltip"] = node["tooltip"].replace("represses", "activates and represses")
                        item.indegree.remove(first_node)
                    elif first_node["regulatoryEffect"] == "repressor" and node["regulatoryEffect"] == "activator":
                        node["regulatoryEffect"] = "dual"
                        node["tooltip"] = node["tooltip"].replace("activates", "activates and represses")
                        item.indegree.remove(first_node)
    return item

