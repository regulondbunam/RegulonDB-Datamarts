import multigenomic_api
import re


def get_pmids(gen_object):
    pmids = []
    for citation in gen_object.citations:
        if citation.publications_id:
            publication = multigenomic_api.publications.find_by_id(citation.publications_id)
            if publication:
                if publication.pmid:
                    if publication.pmid not in pmids:
                        pmids.append(publication.pmid)
    citations_pattern = re.compile("(\[[0-9]+\])")
    if gen_object.note:
        pmids_search = re.findall(citations_pattern, gen_object.note)
        pmids_search = list(set(pmids_search))
        for pmid in pmids_search:
            pmid = pmid[1:-1]
            if pmid not in pmids:
                pmids.append(pmid)
    return ";".join(pmids)
