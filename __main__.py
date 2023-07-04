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
    srna_datamarts, \
    dttDatamart, \
    regulatory_network_datamart, \
    listPage_dm


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

    args = parser.parse_args()

    return args


def create_json(objects_collections, filename, output):
    filename = os.path.join(output, filename)
    with open("{}.json".format(filename), 'w') as json_file:
        json.dump(objects_collections, json_file, indent=2, sort_keys=True)


def write_summary(datamarts_info):
    with open("DatamartExtractorSummary.md", 'w') as out_file:
        text = "# DatamartsExtractorSummary \n" \
               "Creation date: " + datetime.today().strftime('%Y-%m-%d') + "\n \n" \
                "Creation time: " + datetime.today().strftime('%H:%M:%S') + "\n \n" \
                "RegulonDB Version: 12.0 \n" \
                "\n" \
                "## RegulonDB Datamarts Summary \n" \
                "" + datamarts_info
        out_file.write(text)


if __name__ == '__main__':
    arguments = load_arguments_parser()
    multigenomic_api.connect(arguments.database, arguments.url)

    datamart_files = dict()

    datamart_files["geneDatamart"] = gene_datamarts.all_genes_datamarts()
    datamart_files["operonDatamart"] = operon_datamarts.all_operon_datamarts()
    datamart_files["regulonDatamart"] = regulon_datamart.all_regulon_datamarts()
    datamart_files["sigmulonDatamart"] = sigmulon_datamarts.all_sigmulon_datamarts()
    datamart_files["srnaDatamart"] = srna_datamarts.all_srna_datamarts()
    datamart_files["dnaFeatures"] = dttDatamart.all_dtt_datamarts()
    datamart_files["regulatoryNetworkDatamart"] = regulatory_network_datamart.all_regulatory_network_nodes()
    datamart_files["listPage"] = listPage_dm.get_all_list_page_docs()

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
    write_summary(datamartsData)
