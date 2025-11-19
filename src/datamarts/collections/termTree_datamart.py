import multigenomic_api
from datetime import datetime
import re
from src.datamarts.domain.general.remove_items import remove_empty_items
from src.datamarts.domain.general.print_progress import print_progress


class TermsDatamart:

    @property
    def objects(self):
        terms_objects = multigenomic_api.terms.get_all()
        for n, term in enumerate(terms_objects, start=0):
            print_progress(n + 1, len(terms_objects), "Term Tree")
            yield TermsDatamart.TermFamilyBranch(term)
        del terms_objects

    class TermFamilyBranch:
        def __init__(self, term):
            self.term = term
            self.subclass_of = term.sub_class_of
            self.subclasses = term.super_class_of

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

        def to_row(self):
            description = None
            if self.term.ontologies_id == "RDBONTOLGON00001":
                description = self.term.definition.text
            elif self.term.ontologies_id == "RDBONTOLMCO00001":
                description = self.term.description
            term_dict = {
                "_id": self.term.id,
                "name": self.term.name,
                "subclasses": self.subclasses,
                "subclassOf": self.subclass_of,
                "description": description,
                "ontologyId": self.term.ontologies_id,
                "iri": self.term.iri,
                "oboId": self.term.obo_id,
                "genes": get_members_list(self.term.members),
                "synonyms": self.term.synonyms,
                "hasRelatedSynonyms": self.term.has_related_synonyms,
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
