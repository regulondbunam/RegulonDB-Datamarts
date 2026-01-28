import pymongo
from pymongo import errors
import json
import arguments
import os
import ijson
from decimal import Decimal
from datetime import datetime


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


def normalize_doc(obj):
    if isinstance(obj, dict):
        return {k: normalize_doc(v) for k, v in obj.items()}
    elif isinstance(obj, list):
        return [normalize_doc(v) for v in obj]
    elif isinstance(obj, Decimal):
        return float(obj)
    else:
        return obj


def insert_docs(json_file, col, batch_size=1000):
    batch = []
    count = 0
    with open(json_file, "rb") as f:
        for doc in ijson.items(f, f"{col.name}.item"):
            doc = normalize_doc(doc)
            batch.append(doc)
            count += 1
            if len(batch) >= batch_size:
                flush_batch(col, batch)
                batch.clear()
            if count % 100_000 == 0:
                print(f"{count:,} inserted into {col.name}")
        if batch:
            flush_batch(col, batch)


def flush_batch(col, batch):
    try:
        col.insert_many(batch, ordered=False)
    except pymongo.errors.BulkWriteError as e:
        for err in e.details.get("writeErrors", []):
            if err["code"] != 11000:
                print("Insert error:", err)


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

    print(f"[{datetime.now():%Y-%m-%d %H:%M:%S}] creation of collection and adding documents")

    for filename in os.listdir(args.schemas):
        f = os.path.join(args.schemas, filename)
        create_collection_with_file(f)

    for filename in os.listdir(args.collection_data):
        if filename.startswith('.'):
            continue
        collection_name = os.path.splitext(filename)[0]
        collection = db[collection_name]
        f = os.path.join(args.collection_data, filename)
        print(f"Inserting {f} into collection '{collection_name}'")
        insert_docs(f, collection)

    if hasattr(args, "indexes_file"):
        create_indexes_from_file(args.indexes_file)
    else:
        print("⚠ No indexes_file argument provided; skipping index creation.")

    print(f"[{datetime.now():%Y-%m-%d %H:%M:%S}] finished upload process")

else:
    print(__name__)
