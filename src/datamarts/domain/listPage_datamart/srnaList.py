import multigenomic_api
from src.datamarts.domain.general.remove_items import remove_empty_items
from src.datamarts.domain.srna_datamart.regulatory_interactions import RegulatoryInteractions
from src.datamarts.domain.srna_datamart.summary import Summary


class ListSrnaDM:

    @property
    def objects(self):
        srna_objects = multigenomic_api.products.get_all()
        for srna_obj in srna_objects:
            if srna_obj.type == "small RNA":
                print(srna_obj.id)
                regulon_datamart = ListSrnaDM.ListSrna(srna_obj)
                yield regulon_datamart
        del srna_objects

    class ListSrna:

        def __init__(self, srna):
            self.id = srna.id
            self.srna = srna
            self.genes = srna.genes_id
            self.summary = srna

        @property
        def genes(self):
            return self._genes

        @genes.setter
        def genes(self, genes_id):
            self._genes = []
            gene = multigenomic_api.genes.find_by_id(genes_id)
            self._genes.append(gene.name)

        @property
        def summary(self):
            return self._summary

        @summary.setter
        def summary(self, srna):
            self._summary = {}
            reg_ints = multigenomic_api.regulatory_interactions.find_by_regulator_id(srna.id)
            final_reg_int = []
            for ri in reg_ints:
                final_reg_int.append(RegulatoryInteractions(ri).to_dict())
            summary = Summary(srna, final_reg_int).to_dict()
            self._summary = summary

        def to_dict(self):
            regulon_datamart = {
                "_id": self.id,
                "name": self.srna.name,
                "synonyms": self.srna.synonyms,
                "encodedGenes": self.genes,
                "summary": self.summary,
                "datamartType": "srna"
            }
            return regulon_datamart


def list_all_srna_datamarts():
    list_srnas = ListSrnaDM()
    json_srnas = []
    for regulon in list_srnas.objects:
        regulon_dict = remove_empty_items(regulon.to_dict().copy())
        json_srnas.append(remove_empty_items(regulon_dict))
    return json_srnas
