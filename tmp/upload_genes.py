from pymongo import MongoClient
import json

#  client = MongoClient("mongodb://pablo:pablo123!@ds211259.mlab.com:11259/datamarts?retryWrites=false")
client = MongoClient("mongodb://localhost:27017/")
database = client["regulondbdatamarts"]

collection = database["geneCoexpressions"]

#gene_datamart = open("/Users/pablo/Proyectos/RegulonDB/Results/datamarts/genes.json", "r")
gene_datamart = open("/Users/pablo/Dropbox (UNAM-CCG)/PGC_CO/Proyectos/12.UBI/RegulonDB-Datamarts/4.Seguimiento/CoExpression/coexpression_example.json", "r")
genes = json.loads(gene_datamart.read())
genes = genes["geneCoexpressions"]
for gene in genes:
    collection.insert_one(gene)
