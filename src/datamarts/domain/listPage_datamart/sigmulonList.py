import multigenomic_api
from src.datamarts.domain.general.remove_items import remove_empty_items
from src.datamarts.domain.sigmulon_datamart.statistics import Statistics


class ListSigmaFactorDM:

    @property
    def objects(self):
        sigma_factor_objects = multigenomic_api.sigma_factors.get_all()
        for sigma_factor in sigma_factor_objects:
            print(sigma_factor.id)
            sigma_factor_datamart = ListSigmaFactorDM.ListSigmaFactor(sigma_factor)
            yield sigma_factor_datamart
        del sigma_factor_objects

    class ListSigmaFactor:

        def __init__(self, sigma_factor):
            self.id = sigma_factor.id
            self.sigma_factor = sigma_factor
            self.statistics = sigma_factor.id
            self.sigmulon_gene_name = sigma_factor.genes_id

        @property
        def statistics(self):
            return self._statistics

        @statistics.setter
        def statistics(self, sigmulon_id):
            self._statistics = Statistics(sigmulon_id).to_dict()

        @property
        def sigmulon_gene_name(self):
            return self._sigmulon_gene_name

        @sigmulon_gene_name.setter
        def sigmulon_gene_name(self, genes_id):
            gene = multigenomic_api.genes.find_by_id(genes_id)
            self._sigmulon_gene_name = gene.name

        def to_dict(self):
            sigma_factor_datamart = {
                "_id": self.id,
                "name": self.sigma_factor.name,
                "synonyms": self.sigma_factor.synonyms,
                "statistics": self._statistics,
                "sigmulonGeneName": self.sigmulon_gene_name,
                "datamartType": "sigmulon"
            }
            return sigma_factor_datamart


def list_all_sigma_factor_datamarts():
    list_sigma_factors = ListSigmaFactorDM()
    json_sigma_factors = []
    for sigma_factor in list_sigma_factors.objects:
        sigma_factor_dict = remove_empty_items(sigma_factor.to_dict().copy())
        json_sigma_factors.append(remove_empty_items(sigma_factor_dict))
    return json_sigma_factors
