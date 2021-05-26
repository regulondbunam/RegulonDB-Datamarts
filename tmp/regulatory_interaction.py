from src import multigenomic_api

multigenomic_api.initialize_db('multigenomic', 'mongodb://pablo:pablo@127.0.0.1:27017')

for regulator_object in multigenomic_api.regulatory_interaction.find_regulators_by_regulated_entity_ids(["RDBECOLIPM02438"]):
    print(regulator_object.to_dict())

print(multigenomic_api.regulatory_interaction.find_by_regulated_entity_ids([]))