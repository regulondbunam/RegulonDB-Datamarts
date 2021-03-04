import multigenomic_api
from src.datamarts.collections import operon_datamarts
import os
import json

multigenomic_api.connect("regulondbmultigenomic","mongodb://release:regulondb_2020@127.0.0.1:27017/regulondbmultigenomic")

datamart_files = dict()

datamart_files["operonDatamart"] = operon_datamarts.all_operon_datamarts()


def create_json(objects, filename, output):
    filename = os.path.join(output, filename)
    with open("{}.json".format(filename), 'w') as json_file:
        json.dump(objects, json_file, indent=2, sort_keys=True)


for collection_name, objects in datamart_files.items():
    print("Writing {} json file,".format(collection_name))
    objects_to_json = []
    for object in objects:
        objects_to_json.append(object.copy())
    print("\t total of {} objects: {}".format(collection_name, len(objects_to_json)))
    objects_to_json = {
        collection_name: objects_to_json
    }
    create_json(objects_to_json, collection_name, "/Users/pablo/Proyectos/RegulonDB/Results/datamarts")