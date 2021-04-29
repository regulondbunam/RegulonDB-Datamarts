import multigenomic_api
from src.datamarts.domain.sigmulon_datamart.sigma_factor import SigmaFactor

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
            self.statistics = sigma_factor

        @property
        def sigma_factor(self):
            return self._sigma_factor

        @sigma_factor.setter
        def sigma_factor(self, sigma_factor):
            self._sigma_factor = SigmaFactor(sigma_factor).to_dict()

        def to_dict(self):
            sigmulon_datamart = {
                "_id": self.id,
                "sigmaFactor": self.sigma_factor,
                # "transcribedPromoters": self.transcribed_promoters,
                # "statistics": self.statistics,
                # "citations": self.citations
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
