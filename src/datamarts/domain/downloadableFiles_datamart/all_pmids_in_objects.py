import pymongo
import re
from datetime import datetime


def find_doc_by_id(collection_name, doc_id):
    """Generic function to find a document by ID."""
    client = pymongo.MongoClient("mongodb://andresloal:15091052@localhost:27017/")
    db = client["regulondbmultigenomic"]
    collection = db[collection_name]
    return collection.find_one({"_id": doc_id})


def extract_citations(doc, obj_name, obj_type):
    """Extract citations from a document."""
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
            if pub_pmid not in pub_to_evidences:
                pub_to_evidences[pub_pmid] = []
            if ev_code:
                pub_to_evidences[pub_pmid].append(ev_code)

    for pub_pmid, ev_codes in pub_to_evidences.items():
        ev_codes_str = ",".join(ev_codes)
        items.append(f"{pub_pmid}\t{obj_type}\t{doc['_id']}\t{obj_name}\t\t{ev_codes_str}")

    return items


def extract_notes(doc, obj_name, obj_type):
    """Extract PMIDs from notes."""
    items = []
    citations_pattern = re.compile(r"\[[0-9]+\]")
    pmids = set(re.findall(citations_pattern, doc.get("note", "")))

    for pmid in pmids:
        pmid = pmid.strip("[]")
        items.append(f"{pmid}\t{obj_type}\t{doc['_id']}\t{obj_name}\tnote")

    return items


def process_collection(collection_name, obj_type):
    """Process a MongoDB collection and extract relevant rows."""
    client = pymongo.MongoClient("mongodb://andresloal:15091052@localhost:27017/")
    db = client["regulondbmultigenomic"]
    collection = db[collection_name]

    items = []
    for doc in collection.find():
        obj_name = doc.get("abbreviatedName", doc.get("name", ""))
        items.extend(extract_citations(doc, obj_name, obj_type))
        items.extend(extract_notes(doc, obj_name, obj_type))

    return items


def process_regulatory_interactions():
    """Process regulatory interactions with special handling for regulatory sites."""
    client = pymongo.MongoClient("mongodb://andresloal:15091052@localhost:27017/")
    db = client["regulondbmultigenomic"]
    collection = db["regulatoryInteractions"]

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
    """Generate a document containing all extracted PMIDs."""
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

    pmids_doc = {
        "_id": "RDBECOLIDLF00021",
        "fileName": "allObjectPmids_internal",
        "title": "Complete Object PMIDs Set",
        "fileFormat": "rif-version 1",
        "license": (
            "# RegulonDB is free for academic/noncommercial use\n"
            "# User is not entitled to change or erase data sets of the RegulonDB\n"
            "# database or to eliminate copyright notices from RegulonDB. Furthermore,\n"
            "# User is not entitled to expand RegulonDB or to integrate RegulonDB partly\n"
            "# or as a whole into other databank systems, without prior written consent\n"
            "# from CCG-UNAM.\n"
            "# Please check the license at https://regulondb.ccg.unam.mx/manual/aboutUs/terms-conditions"
        ),
        "citation": citation,
        "contact": {
            "person": "RegulonDB Team",
            "webPage": None,
            "email": "regulondb@ccg.unam.mx",
        },
        "version": "1.0",
        "creationDate": creation_date.strftime("%m-%d-%Y"),
        "columnsDetails": (
            "# Columns:\n"
            "# (1) PMID\n"
            "# (2) Type of the biological object\n"
            "# (3) Id of the object\n"
            "# (4) Name of the object\n"
            "# (5) Origin\n"
            "# (6) Evidence codes"
        ),
        "content": " \n".join(item_list),
        "rdbVersion": rdb_version,
        "description": "Collection of all objects with their PMIDs.",
        "group": "EVIDENCE",
    }

    return pmids_doc
