import multigenomic_api
from random import randint

from src.datamarts.domain.dnaFeatures_datamart.gene_items import GeneDnaFeatures
from src.datamarts.domain.dnaFeatures_datamart.promoter_items import PromoterDnaFeatures
from src.datamarts.domain.dnaFeatures_datamart.terminator_items import TerminatorDNAFeatures
from src.datamarts.domain.dnaFeatures_datamart.reg_int_items import RegIntDnaFeatures


def all_dtt_datamarts():
    dict_colors = asign_colors_to_gene_multifun_type()
    gene_items = GeneDnaFeatures(dict_colors)
    promoter_items = PromoterDnaFeatures(dict_colors)
    terminator_items = TerminatorDNAFeatures()
    reg_int_items = RegIntDnaFeatures(dict_colors)
    json_items = []

    for gene in gene_items.objects:
        gene = remove_none_fields_empty_lists(gene.to_dict())
        json_items.append(gene.copy())
    for promoter in promoter_items.objects:
        promoter = remove_none_fields_empty_lists(promoter.to_dict())
        json_items.append(promoter.copy())
    for terminator in terminator_items.objects:
        terminator = remove_none_fields_empty_lists(terminator.to_dict())
        json_items.append(terminator.copy())
    for reg_int in reg_int_items.objects:
        reg_int = remove_none_fields_empty_lists(reg_int.to_dict())
        json_items.append(reg_int.copy())
    return json_items


def asign_colors_to_gene_multifun_type():
    colors = {}
    used_colors = []
    multifun_term_items = multigenomic_api.terms.get_multifun_terms()
    for multifun_item in multifun_term_items:
        rgb_color = ""
        while True:
            rgb_color = f"{randint(0, 255)}, {randint(0, 255)}, {randint(0, 255)}"
            if rgb_color not in used_colors:
                break
        used_colors.append(rgb_color)
        colors.setdefault(multifun_item.id, rgb_color)
    return colors


def remove_none_fields_empty_lists(gene_object):
    if isinstance(gene_object, dict):
        return {property: remove_none_fields_empty_lists(property_value) for property, property_value in gene_object.items() if property_value}
    elif isinstance(gene_object, list):
        if len(gene_object) != 0:
            return [remove_none_fields_empty_lists(v) for v in gene_object]
    else:
        return gene_object