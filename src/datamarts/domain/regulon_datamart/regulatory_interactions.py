import multigenomic_api
from src.datamarts.domain.general.biological_base import BiologicalBase

class RegulatoryInteractions(BiologicalBase):
    def __init__(self, products_ids):
        regulatory_interactions = multigenomic_api.regulatory_interactions.find_regulators_by_regulated_entity_ids(products_ids)
        print(len(regulatory_interactions))

    def to_dict(self):
        regulatory_interactions = {

        }