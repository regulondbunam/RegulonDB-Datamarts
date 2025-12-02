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
        print("Creating " + collection_name + " collection")
        db.create_collection(collection_name, **commands)
    except pymongo.errors.OperationFailure as e:
        print(e.code)
        print(e.details)
        print(collection_name + " collection already exists")


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


def create_indexes_from_file(index_file):
    with open(index_file, "r", encoding="utf-8") as f:
        data = json.load(f)

    commands = data.get("commands", [])
    print("\n=== Creating Indexes ===")

    for cmd in commands:
        collection = cmd.get("createIndexes")
        if not collection:
            print("Invalid command:", cmd)
            continue

        print(f"➡ Creating indexes for: {collection}")
        try:
            result = db.command(cmd)
            print(result)
        except Exception as e:
            print(f"❌ Error creating index on {collection}: {e}")

    print("✔ Finished creating indexes.\n")


if __name__ == "__main__":
    args = arguments.load_arguments()
    db = connect_db(args.url, args.database)

    print("Starting creation of collection and adding documents")

    for filename in os.listdir(args.schemas):
        f = os.path.join(args.schemas, filename)
        create_collection_with_file(f)

    for filename in os.listdir(args.collection_data):
        if filename.startswith('.'):
            continue
        f = os.path.join(args.collection_data, filename)
        print(f"inserting docs of {f}")
        insert_docs(f)

    if hasattr(args, "indexes_file"):
        create_indexes_from_file(args.indexes_file)
    else:
        print("⚠ No indexes_file argument provided; skipping index creation.")

else:
    print(__name__)
