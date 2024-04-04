import multigenomic_api
from datetime import datetime
import re
from src.datamarts.domain.general.remove_items import remove_empty_items


class TermsDatamart:

    @property
    def objects(self):
        terms_objects = multigenomic_api.terms.get_go_terms()
        for term in terms_objects:
            print(term.id)
            ri_row = TermsDatamart.TermFamilyBranch(term)
            yield ri_row
        del terms_objects

    class TermFamilyBranch:
        def __init__(self, term):
            self.term = term
            self.subclass_of = term.sub_class_of
            self.subclasses = term.super_class_of
            self.ontology_id = term.external_cross_references

        @property
        def subclass_of(self):
            return self._subclass_of

        @subclass_of.setter
        def subclass_of(self, subclass_of):
            self._subclass_of = subclass_of

        @property
        def subclasses(self):
            return self._subclasses

        @subclasses.setter
        def subclasses(self, subclasses):
            self._subclasses = subclasses

        @property
        def ontology_id(self):
            return self._ontology_id

        @ontology_id.setter
        def ontology_id(self, ext_cross_refs):
            self._ontology_id = ""
            for ref in self.term.external_cross_references:
                if re.match(r'GO:[0-9]*',ref.object_id):
                    self._ontology_id = ref.object_id

        def to_row(self):
            term_dict = {
                "_id": self.term.id,
                "name": self.term.name,
                "subclasses": self.subclasses,
                "subclassOf": self.subclass_of,
                "description": self.term.definition.text,
                "ontologyId": self.ontology_id,
                "genes": get_members_list(self.term.members)
            }
            return term_dict


def get_tree():
    terms = TermsDatamart()
    terms_family_tree = []
    for term in terms.objects:
        terms_family_tree.append(remove_empty_items(term.to_row().copy()))
    return terms_family_tree


def get_members_list(members):
    members_list = []
    if members:
        if members.products:
            for prod_id in members.products:
                product = multigenomic_api.products.find_by_id(prod_id)
                gene = multigenomic_api.genes.find_by_id(product.genes_id)
                members_list.append({
                    "_id": gene.id,
                    "name": gene.name,
                    "productName": product.name
                })
    return members_list
