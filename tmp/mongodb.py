from pymongo import MongoClient

cursor = MongoClient("mongodb://localhost:27017/")

multigenomic_db = cursor["multigenomic"]

collection_info = multigenomic_db.command({'listCollections': 1, 'filter': {'name': 'gene_descripted'}})

print(collection_info)