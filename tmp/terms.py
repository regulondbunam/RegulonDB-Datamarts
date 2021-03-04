from src import multigenomic_api

multigenomic_api.initialize_db('multigenomic', 'mongodb://pablo:pablo@127.0.0.1:27017')

print(multigenomic_api.term.get_associated_members_ids("RDBECOLITR00001", "products"))