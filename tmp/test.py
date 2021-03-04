from src import multigenomic_api

multigenomic_api.initialize_db('multigenomic', 'mongodb://pablo:pablo@127.0.0.1:27017')

pmid_length = []
for gene in multigenomic_api.publication.get_all():
    if gene.pmid:
        pmid_length.append(len(gene.pmid))
        if len(gene.pmid) < 5:
            print(gene.pmid)

print(set(pmid_length))

with open("ff", 'r') as fp:
    hola = fp.read()