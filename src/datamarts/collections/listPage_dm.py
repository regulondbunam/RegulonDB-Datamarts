from src.datamarts.domain.general.remove_items import remove_empty_items
from src.datamarts.domain.listPage_datamart import regulonList, geneList, operonList, sigmulonList


def get_all_list_page_docs():
    list_dm = []

    regulons = regulonList.list_all_regulon_datamarts()
    for regulon in regulons:
        regulon_dict = remove_empty_items(regulon.copy())
        list_dm.append(remove_empty_items(regulon_dict))

    genes = geneList.list_all_gene_datamarts()
    for gene in genes:
        gene_dict = remove_empty_items(gene.copy())
        list_dm.append(remove_empty_items(gene_dict))

    operons = operonList.list_all_operon_datamarts()
    for operon in operons:
        operon_dict = remove_empty_items(operon.copy())
        list_dm.append(remove_empty_items(operon_dict))

    sigmulons = sigmulonList.list_all_sigma_factor_datamarts()
    for sigmulon in sigmulons:
        sigmulon_dict = remove_empty_items(sigmulon.copy())
        list_dm.append(remove_empty_items(sigmulon_dict))

    return list_dm
