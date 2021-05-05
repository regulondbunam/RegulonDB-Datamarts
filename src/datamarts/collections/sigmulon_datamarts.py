import multigenomic_api
from src.datamarts.domain.sigmulon_datamart.sigma_factor import SigmaFactor
from src.datamarts.domain.sigmulon_datamart.transcribed_promoters import TranscribedPromoters
from src.datamarts.domain.sigmulon_datamart.statistics import Statistics
from src.datamarts.domain.general.biological_base import BiologicalBase


class SigmulonDatamarts:

    @property
    def objects(self):
        sigma_factor_objects = multigenomic_api.sigma_factors.get_all()
        for sigma_factor in sigma_factor_objects:
            sigma_factor_datamart = SigmulonDatamarts.SigmulonDatamart(sigma_factor)
            yield sigma_factor_datamart
        del sigma_factor_objects

    class SigmulonDatamart:

        def __init__(self, sigma_factor):
            self.id = sigma_factor.id
            self.sigma_factor = sigma_factor
            self.transcribed_promoters = sigma_factor
            self.statistics = sigma_factor.id

        @property
        def sigma_factor(self):
            return self._sigma_factor

        @sigma_factor.setter
        def sigma_factor(self, sigma_factor):
            self._sigma_factor = SigmaFactor(sigma_factor).to_dict()

        @property
        def transcribed_promoters(self):
            return self._transcribed_promoters

        @transcribed_promoters.setter
        def transcribed_promoters(self, sigma_factor):
            self._transcribed_promoters = []
            trans_promoters = get_transcribed_promoters(sigma_factor.id)
            for promoter_id in trans_promoters:
                promoter = multigenomic_api.promoters.find_by_id(promoter_id)
                promoter = TranscribedPromoters(promoter).to_dict()
                self._transcribed_promoters.append(promoter)

        @property
        def statistics(self):
            return self._statistics

        @statistics.setter
        def statistics(self, sigmulon_id):
            self._statistics = Statistics(sigmulon_id).to_dict()

        def to_dict(self):
            sigmulon_datamart = {
                "_id": self.id,
                "sigmaFactor": self.sigma_factor,
                "transcribedPromoters": self.transcribed_promoters,
                "statistics": self.statistics,
                "allCitations": BiologicalBase.get_all_citations()
            }
            return sigmulon_datamart


def all_sigmulon_datamarts():
    sigmulon = SigmulonDatamarts()
    json_sigmulon = []
    for sigma in sigmulon.objects:
        sigma_dict = remove_none_fields_empty_lists(sigma.to_dict().copy())
        json_sigmulon.append(remove_none_fields_empty_lists(sigma_dict))
    return json_sigmulon


def remove_none_fields_empty_lists(sigma_object):
    if isinstance(sigma_object, dict):
        return {property: remove_none_fields_empty_lists(property_value) for property, property_value in sigma_object.items() if property_value}
    elif isinstance(sigma_object, list):
        if len(sigma_object) != 0:
            return [remove_none_fields_empty_lists(v) for v in sigma_object]
    else:
        return sigma_object


def get_transcribed_promoters(sigma_factor_id):
    trans_promoters = []
    promoters = multigenomic_api.promoters.find_by_sigma_factor_id(sigma_factor_id)
    for promoter in promoters:
        trans_units = multigenomic_api.transcription_units.find_by_promoter_id(promoter.id)
        if len(trans_units) > 0:
            trans_promoters.append(promoter.id)
    return trans_promoters
