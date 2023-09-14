import json

from src.datamarts.domain.gensorUnit_datamart.gene_ontology import GeneOntology

from src.datamarts.domain.general.remove_items import remove_empty_items


def load_gus():
    with open("./lib/GensorUnits.json", 'r') as json_gus:
        rules = json.load(json_gus)
        return rules


class GensorUnitsDatamarts:

    @property
    def objects(self):
        gensor_units = load_gus()
        for gu in gensor_units:
            print(gu["_id"])
            gu_dm = GensorUnitsDatamarts.GensorUnitDatamart(gu)
            yield gu_dm
        del gensor_units

    class GensorUnitDatamart:
        def __init__(self, gensor_unit):
            self.gu = gensor_unit
            self.gene_ontology = gensor_unit

        @property
        def gene_ontology(self):
            return self._gene_ontology

        @gene_ontology.setter
        def gene_ontology(self, gensor_unit):
            self._gene_ontology = GeneOntology(gensor_unit).to_dict()

        def to_dict(self):
            self.gu["gensorUnit"]["geneOntology"] = self._gene_ontology
            return self.gu


def all_gensor_unit_datamarts():
    gus = GensorUnitsDatamarts()
    json_sigmulon = []
    for gu in gus.objects:
        gu_dict = remove_empty_items(gu.to_dict().copy())
        json_sigmulon.append(remove_empty_items(gu_dict))
    return json_sigmulon
