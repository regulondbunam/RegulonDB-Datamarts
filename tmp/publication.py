from src import multigenomic_api

multigenomic_api.initialize_db('multigenomic', 'mongodb://pablo:pablo@127.0.0.1:27017')
multigenomic_api.publication.find_by_pmid('12068631')

'''
3923479
7862644
4139702
7510397
12762060
18332132
24812390
9353919
'''