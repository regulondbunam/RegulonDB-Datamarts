import json
import os

ri_phrases_json = "/Users/pablo/Dropbox (UNAM-CCG)/PGC_CO/Proyectos/12.UBI/6RegulonDB-Ana/1.Entrada/ri.json"

with open(ri_phrases_json, 'r') as ri_fp:
    ri_phrases = ri_fp.read()
    ri_objects = json.loads(ri_phrases)
    for ri_object in ri_objects["ri"]:
        for ri_id, properties in ri_object.items():
            sources = properties["sources"][0]
            for property_name, phrases in sources.items():
                if "general.tfBinding" not in property_name:
                    print(ri_id)
