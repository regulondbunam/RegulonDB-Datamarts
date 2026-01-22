# all_pmids_in_objects.py
import pymongo
import re
from datetime import datetime

mongo_client = None
mongo_db = None

def init_mongo(url):
    global mongo_client, mongo_db
    if mongo_client is None:
        mongo_client = pymongo.MongoClient(url)
        mongo_db = mongo_client["regulondbmultigenomic"]
    return mongo_db


def find_doc_by_id(collection_name, doc_id):
    collection = mongo_db[collection_name]
    return collection.find_one({"_id": doc_id})


def extract_citations(doc, obj_name, obj_type):
    items = []
    pub_to_evidences = {}

    for citation in doc.get("citations", []):
        pub_pmid = ""
        ev_code = ""

        if "publications_id" in citation:
            publication = find_doc_by_id("publications", citation["publications_id"])
            if publication and "pmid" in publication:
                pub_pmid = publication["pmid"]

        if "evidences_id" in citation:
            evidence = find_doc_by_id("evidences", citation["evidences_id"])
            if evidence and "code" in evidence:
                ev_code = evidence["code"]

        if pub_pmid:
            pub_to_evidences.setdefault(pub_pmid, [])
            if ev_code:
                pub_to_evidences[pub_pmid].append(ev_code)

    for pub_pmid, ev_codes in pub_to_evidences.items():
        ev_codes_str = ",".join(ev_codes)
        items.append(f"{pub_pmid}\t{obj_type}\t{doc['_id']}\t{obj_name}\t\t{ev_codes_str}")

    return items


def extract_notes(doc, obj_name, obj_type):
    items = []
    citations_pattern = re.compile(r"\[[0-9]+\]")
    pmids = set(re.findall(citations_pattern, doc.get("note", "")))

    for pmid in pmids:
        pmid = pmid.strip("[]")
        items.append(f"{pmid}\t{obj_type}\t{doc['_id']}\t{obj_name}\tnote")

    return items


def process_collection(collection_name, obj_type):
    collection = mongo_db[collection_name]
    items = []
    for doc in collection.find():
        obj_name = doc.get("abbreviatedName", doc.get("name", ""))
        items.extend(extract_citations(doc, obj_name, obj_type))
        items.extend(extract_notes(doc, obj_name, obj_type))

    return items


def process_regulatory_interactions():
    collection = mongo_db["regulatoryInteractions"]
    items = []

    for doc in collection.find():
        regulator = find_doc_by_id("regulators", doc["regulator"]["_id"])
        obj_name = regulator.get("abbreviatedName", regulator.get("name", "")) if regulator else doc["regulator"]["name"]

        items.extend(extract_citations(doc, obj_name, "ri"))
        items.extend(extract_notes(doc, obj_name, "ri"))

        if "regulatorySites_id" in doc:
            site = find_doc_by_id("regulatorySites", doc["regulatorySites_id"])
            if site:
                items.extend(extract_citations(site, obj_name, "reg_site"))
                items.extend(extract_notes(site, obj_name, "reg_site"))

    return items


def generate_pmids_doc(rdb_version, citation):
    item_list = ["pmid\tobject_type\tobject_id\tobject_name\torigin\tev_code"]

    collections = [
        ("promoters", "promoter"),
        ("genes", "gene"),
        ("products", "product"),
        ("regulators", "regulator"),
        ("regulatoryComplexes", "regulatoryComplex"),
        ("regulatoryContinuants", "regulatoryContinuant"),
        ("sigmaFactors", "sigmaFactor"),
        ("terminators", "terminator"),
        ("transcriptionUnits", "tu"),
    ]

    for collection_name, obj_type in collections:
        item_list.extend(process_collection(collection_name, obj_type))

    item_list.extend(process_regulatory_interactions())

    creation_date = datetime.now()

    return {
        "_id": "RDBECOLIDLF00021",
        "fileName": "allObjectPmids_internal",
        "title": "Complete Object PMIDs Set",
        "fileFormat": "rif-version 1",
        "license": "...",
        "citation": citation,
        "contact": {
            "person": "RegulonDB Team",
            "email": "regulondb@ccg.unam.mx",
        },
        "version": "1.0",
        "creationDate": creation_date.strftime("%m-%d-%Y"),
        "content": " \n".join(item_list),
        "rdbVersion": rdb_version,
        "description": "Collection of all objects with their PMIDs.",
        "group": "EVIDENCE",
    }
