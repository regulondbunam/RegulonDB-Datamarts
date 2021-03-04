class Publication:

    def __init__(self, publication, citation=None):
        self.publication = publication
        self.url = self.publication.url
        self.citation = citation

    @property
    def url(self):
        return self._url

    @url.setter
    def url(self, url):
        self._url = url
        if self._url is None and self.publication.pmid is not None:
            self._url = f"https://www.ncbi.nlm.nih.gov/pubmed/{self.publication.pmid}"

    @property
    def citation(self):
        return self._citation

    @citation.setter
    def citation(self, citation):
        self._citation = None
        if citation is None:
            authors = ", ".join(self.publication.authors)
            citation = f"{authors}, {self.publication.year}, {self.publication.title}"
            self._citation = citation

    def to_dict(self):
        publication = {
            "id": self.publication.id,
            "authors": self.publication.authors,
            "citation": self.citation,
            "pmid": self.publication.pmid,
            "title": self.publication.title,
            "url": self.url,
            "year": self.publication.year
        }
        return publication