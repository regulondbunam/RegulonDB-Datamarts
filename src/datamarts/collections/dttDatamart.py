import multigenomic_api
from random import randint

from src.datamarts.domain.dnaFeatures_datamart.gene_items import GeneDnaFeatures
from src.datamarts.domain.dnaFeatures_datamart.promoter_items import PromoterDnaFeatures
from src.datamarts.domain.dnaFeatures_datamart.terminator_items import TerminatorDNAFeatures
from src.datamarts.domain.dnaFeatures_datamart.reg_int_items import RegIntDnaFeatures

from src.datamarts.domain.general.remove_items import remove_empty_items


def all_dtt_datamarts():
    dict_colors = assign_colors_to_gene_multifun_type()
    gene_items = GeneDnaFeatures(dict_colors)
    promoter_items = PromoterDnaFeatures(dict_colors)
    terminator_items = TerminatorDNAFeatures()
    reg_int_items = RegIntDnaFeatures(dict_colors)
    json_items = []

    for gene in gene_items.objects:
        gene = remove_empty_items(gene.to_dict())
        json_items.append(gene.copy())
    for promoter in promoter_items.objects:
        promoter = remove_empty_items(promoter.to_dict())
        json_items.append(promoter.copy())
    for terminator in terminator_items.objects:
        terminator = remove_empty_items(terminator.to_dict())
        json_items.append(terminator.copy())
    for reg_int in reg_int_items.objects:
        reg_int = remove_empty_items(reg_int.to_dict())
        json_items.append(reg_int.copy())
    return json_items


def assign_colors_to_gene_multifun_type():
    colors = {}
    used_colors = []
    multifun_term_items = multigenomic_api.terms.get_multifun_terms()
    for multifun_item in multifun_term_items:
        rgb_color = ""
        while True:
            rgb_color = f"{randint(150, 255)}, {randint(150, 255)}, {randint(150, 255)}"
            if rgb_color not in used_colors:
                break
        used_colors.append(rgb_color)
        colors.setdefault(multifun_item.id, rgb_color)
    return colors
