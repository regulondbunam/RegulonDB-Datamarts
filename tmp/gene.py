import multigenomic_api
import os
multigenomic_api.initialize_db('multigenomic', 'mongodb://pablo:pablo@127.0.0.1:27017')

print(multigenomic_api.operon.find_by_gene_id("RDBECOLIGN04684").to_json())
genes = multigenomic_api.gene.get_id_by_bnumbers(["b3883", "b1567"])
print(genes["b1567"])

print(os.path.dirname(__file__))
#RDBECOLIGN04635
#RDBECOLIGN04662
#RDBECOLIGN04683