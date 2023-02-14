import multigenomic_api
from src.datamarts.domain.general.remove_items import remove_empty_items


class ListOperonDM:

    @property
    def objects(self):
        operon_objects = multigenomic_api.operons.get_all()
        for operon in operon_objects:
            print(operon.id)
            operon_datamart = ListOperonDM.ListOperon(operon)
            yield operon_datamart
        del operon_objects

    class ListOperon:

        def __init__(self, operon):
            trans_units = multigenomic_api.transcription_units.find_by_operon_id(operon.id)
            self.id = operon.id
            self.operon = operon
            self.trans_units = trans_units
            self.promoters = trans_units
            self.genes = trans_units

        @property
        def promoters(self):
            return self._promoters

        @promoters.setter
        def promoters(self, trans_units):
            self._promoters = []
            for tu in trans_units:
                if tu.promoters_id:
                    self._promoters.append(tu.promoters_id)
            self._promoters = set(list(self._promoters))

        @property
        def genes(self):
            return self._genes

        @genes.setter
        def genes(self, trans_units):
            self._genes = []
            for tu in trans_units:
                if tu.genes_ids:
                    for gene_id in tu.genes_ids:
                        if gene_id not in self._genes:
                            self._genes.append(gene_id)

        def to_dict(self):
            operon_datamart = {
                "_id": self.id,
                "name": self.operon.name,
                "synonyms": self.operon.synonyms,
                "statistics": {
                    "transcriptionUnits": len(self.trans_units),
                    "promoters": len(self.promoters),
                    "genes": len(self.genes)
                },
                "datamartType": "operon"
            }
            return operon_datamart


def list_all_operon_datamarts():
    list_operons = ListOperonDM()
    json_operons = []
    for operon in list_operons.objects:
        operon_dict = remove_empty_items(operon.to_dict().copy())
        json_operons.append(remove_empty_items(operon_dict))
    return json_operons
