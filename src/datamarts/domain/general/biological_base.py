import re

import multigenomic_api

from .external_cross_reference import ExternalCrossReference
from .evidence import Evidence
from .publication import Publication


class BiologicalBase(object):

    all_citations = []

    def __init__(self, external_cross_references, citations, note):
        self.external_cross_references = external_cross_references
        self.citations = citations
        self.formatted_note = note

    @property
    def external_cross_references(self):
        return self._external_cross_references

    @external_cross_references.setter
    def external_cross_references(self, external_cross_references):
        self._external_cross_references = []
        for external_cross_reference in external_cross_references:
            ext_cross_ref_id = external_cross_reference.external_cross_references_id
            object_id = external_cross_reference.object_id
            external_cross_reference = multigenomic_api.external_cross_references.find_by_id(ext_cross_ref_id)
            external_cross_reference = ExternalCrossReference(external_cross_reference, object_id)
            # I added this to avoid duplicates
            if external_cross_reference.to_json() not in self._external_cross_references:
                self._external_cross_references.append(external_cross_reference.to_json())

    @property
    def citations(self):
        return self._citations

    @citations.setter
    def citations(self, citations):
        self._citations = []
        for citation in citations:
            citation_object = {}

            if citation.evidences_id is not None:
                evidence = Evidence(multigenomic_api.evidences.find_by_id(citation.evidences_id))
                citation_object["evidence"] = evidence.to_dict()

            if citation.publications_id is not None:
                publication = Publication(multigenomic_api.publications.find_by_id(citation.publications_id))
                citation_object["publication"] = publication.to_dict()

            self._citations.append(citation_object.copy())
            # all_citations is been shared through all the classes that extends from BiologicalBase
            # having said this, we need a class variable for the all_citations process while processing each
            # datamart document
            # I.E. while working with gene datamart, it will store all the citations from gene, products, motifs,
            # terms, then we will be able to obtain them all and use them in the property all_citations from the
            # datamart as all_citations
            if citation_object not in BiologicalBase.all_citations:
                BiologicalBase.all_citations.append(citation_object.copy())

    @property
    def formatted_note(self):
        return self._formatted_note

    @formatted_note.setter
    def formatted_note(self, note):
        if note is not None:
            if "CITS" in note:
                note = BiologicalBase.format_note(note)
        self._formatted_note = note

    @staticmethod
    def format_note(note):
        citations_pattern = re.compile("(\[[0-9]+\])")
        pmids_search = re.findall(citations_pattern, note)
        pmids_search = list(set(pmids_search))
        unmapped_pmid = False
        for pmid in pmids_search:
            pmid = pmid[1:-1]
            try:
                publication = multigenomic_api.publications.find_by_pmid(pmid)
            except:
                publication = None
            if publication is not None:
                note = note.replace(pmid, publication.id)
                BiologicalBase.add_note_publication_to_all_citations(publication)
            else:
                unmapped_pmid = True
        if unmapped_pmid:
            # TODO: Aqui agregaremos mensaje al log
            # TODO: Agregar mensaje al log sobre posible inconsisetncia en nota, verificar nota antes de insercion
            # note = note.replace(pmid, f"{pmid}[ERROR!]")
            pass
        return note

    @classmethod
    def get_all_citations(cls):
        all_citations = cls.all_citations
        cls.all_citations = []
        return all_citations

    @classmethod
    def add_citation_to_all_citations(cls, citation):
        if citation not in cls.all_citations:
            cls.all_citations.append(citation.copy())

    @classmethod
    def add_note_publication_to_all_citations(cls, publication):
        publication_object = Publication(publication)
        citation_object = {
            "publication": publication_object.to_dict().copy()
        }
        cls.add_citation_to_all_citations(citation_object)
