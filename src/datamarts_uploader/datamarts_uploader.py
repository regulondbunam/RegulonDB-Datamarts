import pymongo
from pymongo import errors
import json
import arguments
import os


# Connecting to MongoDB in localhost and made some test.
def connect_db(url, used_db):
    client = pymongo.MongoClient(url)
    return client[used_db]


def open_json(file):
    with open(file) as file_pointer:
        json_data = json.load(file_pointer)
    return json_data


def create_collection_with_file(file):
    validator = open_json(file)
    collection_name = list(validator.keys())[0]
    commands = validator[collection_name]
    db.drop_collection(collection_name)
    try:
        print("Creating collection with name \"" + collection_name + "\"")
        db.create_collection(collection_name, **commands)
    except pymongo.errors.OperationFailure as e:
        print(e.code)
        print(e.details)
        print("Collection with name \"" + collection_name + "\" already exists")


def insert_docs(doc_file):
    doc = open_json(doc_file)
    for collection_name, collection_documents in doc.items():
        for collection_document in collection_documents:
            try:
                db[collection_name].insert_one(collection_document)
            except pymongo.errors.DuplicateKeyError as dk_error:
                print(f'Document with id {collection_document["_id"]} already exists', dk_error)
            except pymongo.errors.WriteError as inv_document:
                print(f"An error occurs; check the document, {collection_document['_id']}", inv_document)


if __name__ == "__main__":
    args = arguments.load_arguments()
    db = connect_db(args.url, args.database)

    print("principal de create_collection.py")

    for filename in os.listdir(args.schemas):
        f = os.path.join(args.schemas, filename)
        create_collection_with_file(f)

    for filename in os.listdir(args.collection_data):
        f = os.path.join(args.collection_data, filename)
        print(f"inserting docs of {f}")
        insert_docs(f)

else:
    print(__name__)
