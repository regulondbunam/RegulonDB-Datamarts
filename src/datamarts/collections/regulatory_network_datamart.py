from src.datamarts.domain.regulatoryNetwork_datamart.tf_items import RegulatoryNetworkTF
from src.datamarts.domain.regulatoryNetwork_datamart.gene_items import RegulatoryNetworkGene
from src.datamarts.domain.regulatoryNetwork_datamart.srna_items import RegulatoryNetworkSRNA
from src.datamarts.domain.regulatoryNetwork_datamart.compound_items import RegulatoryNetworkCompound

from src.datamarts.domain.general.remove_items import remove_empty_items


def all_regulatory_network_nodes():
    network_nodes_tf = RegulatoryNetworkTF()
    network_nodes_gene = RegulatoryNetworkGene()
    network_nodes_srna = RegulatoryNetworkSRNA()
    network_nodes_compound = RegulatoryNetworkCompound()
    json_nodes = []
    for item in network_nodes_tf.objects:
        item = locate_dual_items(item)
        item = remove_empty_items(item.to_dict())
        json_nodes.append(remove_empty_items(item))
    for item in network_nodes_gene.objects:
        item = locate_dual_items(item)
        item = remove_empty_items(item.to_dict())
        json_nodes.append(remove_empty_items(item))
    for item in network_nodes_srna.objects:
        item = locate_dual_items(item)
        item = remove_empty_items(item.to_dict())
        json_nodes.append(remove_empty_items(item))
    for item in network_nodes_compound.objects:
        item = locate_dual_items(item)
        item = remove_empty_items(item.to_dict())
        json_nodes.append(remove_empty_items(item))
    return json_nodes


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
    try:
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
    except AttributeError:
        print("There is no indegree items")
    return item
