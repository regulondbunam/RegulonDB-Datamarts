import argparse
import json
import os
from datetime import datetime

import multigenomic_api

from src.datamarts.collections import gene_datamarts, operon_datamarts, regulon_datamart, sigmulon_datamarts


def load_arguments_parser():
    parser = argparse.ArgumentParser(description="Extract and Transforms Datamarts from Multigenomic DB.",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter
                                     )
    parser.add_argument(
        "-db", "--database",
        help="Database where the multigenomic data is been stored",
        metavar="multigenomic",
        default="",
        required=False
    )

    parser.add_argument(
        "-u", "--url",
        help="URL to establish a connection between the process and MongoDB",
        metavar="mongodb://user:pass@localhost:27017/regulondbmultigenomic",
        default="",
        required=False
    )

    arguments = parser.parse_args()

    return arguments


def create_json(objects, filename, output):
    filename = os.path.join(output, filename)
    with open("{}.json".format(filename), 'w') as json_file:
        json.dump(objects, json_file, indent=2, sort_keys=True)


def write_summary(datamarts_info):
    with open("DatamartExtractorSummary.md", 'w') as out_file:
        text = "# DatamartsExtractorSummary \n" \
               "Creation date: " + datetime.today().strftime('%Y-%m-%d') + "\n \n" \
                "Creation time: " + datetime.today().strftime('%H:%M:%S') + "\n \n" \
                "RegulonDB Version: 10.8 \n" \
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
    datamartsData = ""
    for collection_name, objects in datamart_files.items():
        print("Writing {} json file,".format(collection_name))
        objects_to_json = []
        for object in objects:
            objects_to_json.append(object.copy())
        print("\t total of {} objects: {}".format(collection_name, len(objects_to_json)))
        datamartsData = datamartsData + "\n ### {}: \n {} Total Objects Generated".format(collection_name, len(objects_to_json))
        objects_to_json = {
            collection_name: objects_to_json
        }
        create_json(objects_to_json, collection_name, "build")
    write_summary(datamartsData)
