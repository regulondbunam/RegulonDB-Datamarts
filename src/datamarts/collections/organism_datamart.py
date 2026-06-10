import multigenomic_api
from src.datamarts.domain.general.remove_items import remove_empty_items
from src.datamarts.domain.general.print_progress import print_progress


class OrganismDatamarts:

    @property
    def objects(self):
        organism_objects = multigenomic_api.organisms.get_all()
        for n, org_object in enumerate(organism_objects, start=0):
            print_progress(n + 1, len(organism_objects), "Organisms")
            yield OrganismDatamarts.OrganismDatamart(org_object)
        del organism_objects

    class OrganismDatamart:

        def __init__(self, organism):
            print(organism)
            self.id = organism.id
            self.name = organism.name
            self.description = organism.description


        def to_dict(self):
            organism_datamart = {
                "_id": self.id,
                "name": self.name,
                "description": self.description
            }
            return organism_datamart

def all_organism():
    organisms = OrganismDatamarts()
    json_organisms = []
    for organism in organisms.objects:
        organism_dict = remove_empty_items(organism.to_dict().copy())
        json_organisms.append(remove_empty_items(organism_dict))
    return json_organisms