from src import multigenomic_api
from src.multigenomic_api.mongodb.models.gene import Gene
from src.multigenomic_api.mongodb.models.product import Product

multigenomic_api.initialize_db('multigenomic', 'mongodb://pablo:pablo@127.0.0.1:27017')

results = Gene.objects().aggregate(*[{
        '$lookup': {
            'from': "product",
            'localField': '_id',
            'foreignField': 'gene_id',
            'as': 'products'
        }},
        #{"$unwind": "$products"},
        {"$lookup": {
            "from": "motif",
            "localField": "products._id",
            "foreignField": "product_id",
            "as": "motifs"
        }},
        #{"$unwind": "$motifs"},
        {"$project": {"_id": 1, "name": 1, "products._id": 1, "products.name": 1, "motifs._id": 1, "motifs.sequence": 1}},
        {"$count": "passing_scores"}
    ])

