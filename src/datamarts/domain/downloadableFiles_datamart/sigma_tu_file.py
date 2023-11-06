import multigenomic_api
from datetime import datetime


class SigmaTU:

    @property
    def objects(self):
        ri_objects = multigenomic_api.regulatory_interactions.get_all()
        for ri_object in ri_objects:
            # print(ri_object.id)
            ri_row = SigmaTU.SigmaTUDatamart(ri_object)
            yield ri_row
        del ri_objects

    class SigmaTUDatamart:
        def __init__(self, ri):
            self.ri = ri
            self.regulated_tu = ri.regulated_entity
            self.sigma = ri.regulated_entity
            self.function = ri.function
            self.ri_evidences = ri.citations

        @property
        def sigma(self):
            return self._sigma

        @sigma.setter
        def sigma(self, regulated_entity):
            self._sigma = ""
            if regulated_entity.type == "transcriptionUnit":
                tu = multigenomic_api.transcription_units.find_by_id(regulated_entity.id)
                if tu.promoters_id:
                    promoter = multigenomic_api.promoters.find_by_id(tu.promoters_id)
                    if promoter.binds_sigma_factor:
                        sigma = multigenomic_api.sigma_factors.find_by_id(promoter.binds_sigma_factor.sigma_factors_id)
                        self._sigma = sigma

        @property
        def function(self):
            return self._function

        @function.setter
        def function(self, function):
            self._function = ""
            if function == "repressor":
                self._function = '-'
            elif function == "activator":
                self._function = '+'

        @property
        def regulated_tu(self):
            return self._regulated_tu

        @regulated_tu.setter
        def regulated_tu(self, regulated_entity):
            self._regulated_tu = None
            if regulated_entity.type == "transcriptionUnit":
                tu = multigenomic_api.transcription_units.find_by_id(regulated_entity.id)
                if tu.promoters_id:
                    promoter = multigenomic_api.promoters.find_by_id(tu.promoters_id)
                    self._regulated_tu = f"{tu.name}[{promoter.name}]"

        @property
        def ri_evidences(self):
            return self._ri_evidences

        @ri_evidences.setter
        def ri_evidences(self, citations):
            self._ri_evidences = []
            for citation in citations:
                if citation.evidences_id:
                    citation_dict = multigenomic_api.evidences.find_by_id(citation.evidences_id)
                    if citation_dict.code not in self._ri_evidences:
                        self._ri_evidences.append(citation_dict.code)
            self._ri_evidences = ",".join(self._ri_evidences)
            self._ri_evidences = f"[{self._ri_evidences}]"

        def to_row(self):
            row = None
            if self.sigma != "":
                if self.regulated_tu:
                    row = f"{self.sigma.abbreviated_name}" \
                          f"\t{self.regulated_tu}" \
                          f"\t{'+'}" \
                          f"\t{self.ri_evidences}" \
                          f"\t{self.ri.confidence_level or '?'}"
            return row


def remove_similar_items(lista):
    resultado = []
    omitir = set()
    for cadena in lista:
        inicio = cadena[:-1]
        if inicio not in omitir:
            resultado.append(cadena)
            omitir.add(inicio)
    return resultado


def remove_repeated_items_by_different_evidences(list):
    resultado = []
    omitir = []
    for i in range(0, len(list) - 1):
        current_item = list[i]
        next_item = list[i + 1]

        cad_inicial_curr = current_item.rsplit('\t')
        cad_inicial_next = next_item.rsplit('\t')
        print(cad_inicial_curr[:-2], cad_inicial_next[:-2])
        if cad_inicial_curr[:-2] == cad_inicial_next[:-2]:
            cad_inicial_curr = "\t".join(cad_inicial_curr[:-2])
            if cad_inicial_curr not in omitir:
                omitir.append(cad_inicial_curr)
                resultado.append(current_item)
        else:
            cad_inicial_curr = "\t".join(cad_inicial_curr[:-2])
            if cad_inicial_curr not in omitir:
                omitir.append(cad_inicial_curr)
                resultado.append(current_item)
    return resultado


def get_all_rows():
    trans_factors = SigmaTU()
    tfs_content = [
        "1)sigmaName\t2)regulatedName\t3)function\t4)evidences\t5)confidenceLevel"]
    for tf in trans_factors.objects:
        if tf.to_row() is not None:
            tfs_content.append(tf.to_row())
    tfs_content = list(set(tfs_content))
    tfs_content = sorted(remove_similar_items(tfs_content))
    tfs_content = remove_repeated_items_by_different_evidences(tfs_content)
    creation_date = datetime.now()
    tfs_doc = {
        "_id": "RDBECOLIDLF00019",
        "fileName": "NetworkSigmaTU",
        "title": "Complete Sigma-TU Network Set",
        "fileFormat": "rif-version 1",
        "license": "RegulonDB is free for academic/noncommercial use\nUser is not entitled to change or erase data sets of the RegulonDB\ndatabase or to eliminate copyright notices from RegulonDB. Furthermore,\nUser is not entitled to expand RegulonDB or to integrate RegulonDB partly\nor as a whole into other databank systems, without prior written consent\nfrom CCG-UNAM.\nPlease check the license at https://regulondb.ccg.unam.mx/manual/aboutUs/terms-conditions",
        "citation": "Salgado H., Gama-Castro S. et al (2023). RegulonDB 12.0: A Comprehensive resource of transcriptional regulation in E. coli K-12",
        "contact": {
            "person": "RegulonDB Team",
            "webPage": None,
            "email": "regulondb@ccg.unam.mx"
        },
        "version": "1.0",
        "creationDate": f"{creation_date.strftime('%m-%d-%Y')}",
        "columnsDetails": "Columns:\n(1) sigmaName. Sigma Name\n(2) regulatedName. TU regulated by the Sigma Factor (regulated TU)\n(3) function. Regulatory Function of the Regulator on the regulated Gene (+ activator, - repressor, -+ dual, ? unknown)\n(4)riEvidences. Evidence that supports the existence of the regulatory interaction\n(5) confidenceLevel. RI confidence level based on its evidence (Values: Confirmed[C], Strong[S], Weak[W], Unknown[?])",
        "content": " \n".join(tfs_content),
        "rdbVersion": "12.0"
    }
    return tfs_doc
