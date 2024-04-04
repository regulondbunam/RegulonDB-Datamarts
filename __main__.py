import argparse
import json
import os
from datetime import datetime

import multigenomic_api

from src.datamarts.collections import \
    gene_datamarts, \
    operon_datamarts, \
    regulon_datamart, \
    sigmulon_datamarts, \
    dttDatamart, \
    regulatory_network_datamart, \
    listPage_dm, \
    gensorUnit_datamarts, \
    downloadable_files_dm, \
    termTree_datamart


def load_arguments_parser():
    parser = argparse.ArgumentParser(description="Extract and Transforms Datamarts from Multigenomic DB.",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter
                                     )
    parser.add_argument(
        "-db", "--database",
        help="Database where the multigenomic data is been stored",
        metavar="multigenomic",
        default="regulondbmultigenomic",
        required=False
    )

    parser.add_argument(
        "-u", "--url",
        help="URL to establish a connection between the process and MongoDB",
        metavar="mongodb://user:pass@localhost:27017/regulondbmultigenomic",
        default="mongodb://andresloal:15091052@localhost:27017/regulondbmultigenomic",
        required=False
    )

    parser.add_argument(
        "-rdbv", "--rdbversion",
        help="RegulonDB Version",
        metavar="12.5",
        default="12.5",
        required=False
    )

    parser.add_argument(
        "-c", "--citation",
        help="RegulonDB citation for downloadable files",
        metavar="# Heladia Salgado, Socorro Gama-Castro, et al., RegulonDB v12.0: a comprehensive resource of transcriptional regulation in E. coli K-12,\n# Nucleic Acids Research, 2023;, gkad1072, https://doi.org/10.1093/nar/gkad1072",
        default="# Heladia Salgado, Socorro Gama-Castro, et al., RegulonDB v12.0: a comprehensive resource of transcriptional regulation in E. coli K-12,\n# Nucleic Acids Research, 2023;, gkad1072, https://doi.org/10.1093/nar/gkad1072",
        required=False
    )

    args = parser.parse_args()

    return args


def create_json(objects_collections, filename, output):
    filename = os.path.join(output, filename)
    with open("{}.json".format(filename), 'w') as json_file:
        json.dump(objects_collections, json_file, indent=2, sort_keys=True)


def write_summary(datamarts_info, started_time):
    with open("DatamartExtractorSummary.md", 'w') as out_file:
        text = "# DatamartsExtractorSummary \n" \
               "Creation date: " + datetime.today().strftime('%Y-%m-%d') + "\n \n" \
                + started_time + \
                "Creation time: " + datetime.today().strftime('%H:%M:%S') + "\n \n" \
                "RegulonDB Version: " + arguments.rdbversion + "\n" \
                "\n" \
                "## RegulonDB Datamarts Summary \n" \
                "" + datamarts_info
        out_file.write(text)


if __name__ == '__main__':
    arguments = load_arguments_parser()
    multigenomic_api.connect(arguments.database, arguments.url)

    started_time = "Started time: " + datetime.today().strftime('%H:%M:%S') + "\n \n"

    datamart_files = dict()

    datamart_files["geneDatamart"] = gene_datamarts.all_genes_datamarts()
    datamart_files["operonDatamart"] = operon_datamarts.all_operon_datamarts()
    datamart_files["regulonDatamart"] = regulon_datamart.all_regulon_datamarts()
    datamart_files["sigmulonDatamart"] = sigmulon_datamarts.all_sigmulon_datamarts()
    datamart_files["dnaFeatures"] = dttDatamart.all_dtt_datamarts()
    datamart_files["regulatoryNetworkDatamart"] = regulatory_network_datamart.all_regulatory_network_nodes()
    datamart_files["listPage"] = listPage_dm.get_all_list_page_docs()
    datamart_files["gensorUnitDatamart"] = gensorUnit_datamarts.all_gensor_unit_datamarts()
    datamart_files["mcoTree"] = termTree_datamart.get_tree()
    datamart_files["downloadableFilesDatamart"] = downloadable_files_dm.get_all_downloadable_docs(arguments.rdbversion, arguments.citation)

    datamartsData = ""
    for collection_name, objects in datamart_files.items():
        print("Writing {} json file,".format(collection_name))
        objects_to_json = []
        for item in objects:
            objects_to_json.append(item.copy())
        print("\t total of {} objects: {}".format(collection_name, len(objects_to_json)))
        datamartsData = \
            datamartsData + "\n ### {}: \n {} Total Objects Generated".format(collection_name, len(objects_to_json))
        objects_to_json = {
            collection_name: objects_to_json
        }
        create_json(objects_to_json, collection_name, "lib/data")
    write_summary(datamartsData, started_time)
