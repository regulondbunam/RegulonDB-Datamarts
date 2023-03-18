from src.datamarts.domain.general.biological_base import BiologicalBase


class Gene(BiologicalBase):

    def __init__(self, gene):
        super().__init__(gene.external_cross_references, gene.citations, gene.note)
        self.gene = gene
        self.fragments = gene.fragments
        self.terms = gene.terms

    def to_dict(self):
        gene = {
            "bnumber": self.gene.bnumber,
            "centisomePosition": self.gene.centisome_position,
            "citations": self.citations,
            "externalCrossReferences": self.external_cross_references,
            "fragments": self.fragments,
            "gcContent": self.gene.gc_content,
            "_id": self.gene.id,
            "leftEndPosition": self.gene.left_end_position,
            'multifunTerms': self.terms,
            "name": self.gene.name,
            "note": self.formatted_note,
            "rightEndPosition": self.gene.right_end_position,
            "sequence": self.gene.sequence,
            "strand": self.gene.strand,
            "synonyms": self.gene.synonyms,
            # "terms": self.terms,
            "type": self.gene.type
        }
        return gene

    @property
    def terms(self):
        return self._terms

    @terms.setter
    def terms(self, terms):
        self._terms = []
        for term in terms:
            term = {
                '_id': term.terms_id,
                'label': term.term_label,
                'name': term.terms_name
            }
            self._terms.append(term.copy())

    @property
    def fragments(self):
        return self._fragments

    @fragments.setter
    def fragments(self, fragments):
        self._fragments = []
        for fragment in fragments:
            fragment = {
                "_id": fragment.id,
                "name": fragment.name,
                "leftEndPosition": fragment.left_end_position,
                "rightEndPosition": fragment.right_end_position,
                "sequence": fragment.sequence,
                "centisomePosition": fragment.centisome_position
            }
            self._fragments.append(fragment.copy())
