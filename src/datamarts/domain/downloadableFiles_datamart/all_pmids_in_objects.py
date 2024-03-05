import pymongo
import csv
import re
import multigenomic_api
from datetime import datetime


def find_reg_by_id(regulator_id):
    client = pymongo.MongoClient("mongodb://andresloal:15091052@localhost:27017/")
    db = client["regulondbmultigenomic"]
    collection = db["regulators"]
    regulator = collection.find_one({"conformations._id": regulator_id})
    return regulator


def get_pubs_of_note(note, items, doc_id, obj_name, type):
    citations_pattern = re.compile("(\[[0-9]+\])")
    pmids_search = re.findall(citations_pattern, note)
    pmids_search = list(set(pmids_search))
    for pmid in pmids_search:
        pmid = pmid[1:-1]
        cadena = f"{doc_id}\t{obj_name}\t{type}\t{pmid}\tnote"
        if cadena not in items:
            items.append(cadena)
    return items


def get_item_rows(collection, type):
    items = []
    for doc in collection:
        obj_name = ""
        try:
            if doc.abbreviated_name:
                obj_name = doc.abbreviated_name
        except:
            obj_name = doc.name
            pass
        for citation in doc.citations:
            if citation.publications_id:
                publication = multigenomic_api.publications.find_by_id(citation.publications_id)
                if publication:
                    if publication.pmid:
                        cadena = f"{doc.id}\t{obj_name}\t{type}\t{publication.pmid}\t"
                        if cadena not in items:
                            items.append(cadena)
        if doc.note:
            items = get_pubs_of_note(doc.note, items, doc.id, obj_name, type)
    return items


def get_ri_rows(collection):
    items = []
    for doc in collection:
        regulator = find_reg_by_id(doc.regulator.id)
        obj_name = ""
        try:
            if doc.abbreviated_name:
                obj_name = doc.abbreviated_name
        except:
            obj_name = doc.name
            pass
        for citation in doc.citations:
            if citation.publications_id:
                publication = multigenomic_api.publications.find_by_id(citation.publications_id)
                if publication:
                    if publication.pmid:
                        cadena = f"{doc.id}\t{obj_name}\tri\t{publication.pmid}\t"
                        if cadena not in items:
                            items.append(cadena)
        if doc.note:
            items = get_pubs_of_note(doc.note, items, doc.id, obj_name, "ri")
        if doc.regulatory_sites_id:
            site = multigenomic_api.regulatory_sites.find_by_id(doc.regulatory_sites_id)
            for citation in site.citations:
                if citation.publications_id:
                    publication = multigenomic_api.publications.find_by_id(citation.publications_id)
                    if publication:
                        if publication.pmid:
                            cadena = f"{site.id}\t{obj_name}\treg_site\t{publication.pmid}\t"
                            if cadena not in items:
                                items.append(cadena)
            if site.note:
                items = get_pubs_of_note(site.note, items, site.id, obj_name, "reg_site")
    return items


def all_pmids_rows(rdb_version, citation):
    item_list = ["object_id\tobject_name\tobject_type\tpmid\torigin"]

    biological_objects = multigenomic_api.promoters.get_all()
    item_list.extend(get_item_rows(biological_objects, "promoter"))

    biological_objects = multigenomic_api.genes.get_all()
    item_list.extend(get_item_rows(biological_objects, "gene"))

    biological_objects = multigenomic_api.products.get_all()
    item_list.extend(get_item_rows(biological_objects, "product"))

    biological_objects = multigenomic_api.regulators.get_all()
    item_list.extend(get_item_rows(biological_objects, "regulator"))

    biological_objects = multigenomic_api.regulatory_continuants.get_all()
    item_list.extend(get_item_rows(biological_objects, "regulatoryContinuant"))

    biological_objects = multigenomic_api.regulatory_interactions.get_all()
    item_list.extend(get_ri_rows(biological_objects))

    biological_objects = multigenomic_api.sigma_factors.get_all()
    item_list.extend(get_item_rows(biological_objects, "sigmaFactor"))

    biological_objects = multigenomic_api.terminators.get_all()
    item_list.extend(get_item_rows(biological_objects, "terminator"))

    biological_objects = multigenomic_api.transcription_units.get_all()
    item_list.extend(get_item_rows(biological_objects, "tu"))

    creation_date = datetime.now()

    pmids_doc = {
            "_id": "RDBECOLIDLF00021",
            "fileName": "allObjectPmids_internal",
            "title": "Complete Object PMIDs Set",
            "fileFormat": "rif-version 1",
            "license": "# RegulonDB is free for academic/noncommercial use\n# User is not entitled to change or erase data sets of the RegulonDB\n# database or to eliminate copyright notices from RegulonDB. Furthermore,\n# User is not entitled to expand RegulonDB or to integrate RegulonDB partly\n# or as a whole into other databank systems, without prior written consent\n# from CCG-UNAM.\n# Please check the license at https://regulondb.ccg.unam.mx/manual/aboutUs/terms-conditions",
            "citation": citation,
            "contact": {
                "person": "RegulonDB Team",
                "webPage": None,
                "email": "regulondb@ccg.unam.mx"
            },
            "version": "1.0",
            "creationDate": f"{creation_date.strftime('%m-%d-%Y')}",
            "columnsDetails": "# Columns:\n# (1) Id of the object\n# (2) Name of the Object\n# (3) Type of the biological object\n# (4) PMID\n# (5) If the pmid was in the note, this indicates the origin",
            "content": " \n".join(item_list),
            "rdbVersion": rdb_version,
            "description": "Collection of all objets with its pmids.",
            "group": "EVIDENCE"
        }
    return pmids_doc
