import multigenomic_api
from datetime import datetime


class SigmaTU:

    @property
    def objects(self):
        promoter_objects = multigenomic_api.promoters.get_all()
        for promoter_obj in promoter_objects:
            #print(promoter_obj.id)
            ri_row = SigmaTU.SigmaTUDatamart(promoter_obj)
            yield ri_row
        del promoter_objects

    class SigmaTUDatamart:
        def __init__(self, promoter):
            self.promoter = promoter
            self.sigma = promoter.binds_sigma_factor
            self.sigma_evidences = promoter
            self.regulated_tus = promoter

        @property
        def sigma(self):
            return self._sigma

        @sigma.setter
        def sigma(self, binds_sigma_factor):
            self._sigma = ""
            if binds_sigma_factor:
                self._sigma = multigenomic_api.sigma_factors.find_by_id(binds_sigma_factor.sigma_factors_id)

        @property
        def regulated_tus(self):
            return self._regulated_tus

        @regulated_tus.setter
        def regulated_tus(self, promoter):
            self._regulated_tus = []
            trans_units = multigenomic_api.transcription_units.find_by_promoter_id(promoter.id)
            for tu in trans_units:
                self._regulated_tus.append(f"{tu.name}[{promoter.name}]")

        @property
        def sigma_evidences(self):
            return self._sigma_evidences

        @sigma_evidences.setter
        def sigma_evidences(self, promoter):
            self._sigma_evidences = []
            citations = []
            if promoter.binds_sigma_factor:
                citations.extend(promoter.binds_sigma_factor.citations)
            citations.extend(promoter.citations)
            for citation in citations:
                if citation.evidences_id:
                    citation_dict = multigenomic_api.evidences.find_by_id(citation.evidences_id)
                    if citation_dict.code not in self._sigma_evidences:
                        self._sigma_evidences.append(citation_dict.code)
            self._sigma_evidences = ",".join(self._sigma_evidences)
            self._sigma_evidences = f"[{self._sigma_evidences}]"

        def to_row(self):
            response = []
            if self.sigma != "":
                for tu in self.regulated_tus:
                    row = f"{self.sigma.abbreviated_name}" \
                          f"\t{tu}" \
                          f"\t{'+'}" \
                          f"\t{self.sigma_evidences}" \
                          f"\t{self.promoter.confidence_level or '?'}"
                    response.append(row)
            return response


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
        "1)sigmaName\t2)regulatedName\t3)function\t4)promoterEvidence\t5)confidenceLevel"]
    for tf in trans_factors.objects:
        if tf.to_row() is not None:
            tfs_content.extend(tf.to_row())
    tfs_content = list(set(tfs_content))
    tfs_content = sorted(remove_similar_items(tfs_content))
    tfs_content = remove_repeated_items_by_different_evidences(tfs_content)
    creation_date = datetime.now()
    tfs_doc = {
        "_id": "RDBECOLIDLF00019",
        "fileName": "NetworkSigmaTU",
        "title": "Complete Sigma-TU Network Set",
        "fileFormat": "rif-version 1",
        "license": "# RegulonDB is free for academic/noncommercial use\n# User is not entitled to change or erase data sets of the RegulonDB\n# database or to eliminate copyright notices from RegulonDB. Furthermore,\n# User is not entitled to expand RegulonDB or to integrate RegulonDB partly\n# or as a whole into other databank systems, without prior written consent\n# from CCG-UNAM.\n# Please check the license at https://regulondb.ccg.unam.mx/manual/aboutUs/terms-conditions",
        "citation": "# Heladia Salgado, Socorro Gama-Castro, et al., RegulonDB v12.0: a comprehensive resource of transcriptional regulation in E. coli K-12,\n# Nucleic Acids Research, 2023;, gkad1072, https://doi.org/10.1093/nar/gkad1072",
        "contact": {
            "person": "RegulonDB Team",
            "webPage": None,
            "email": "regulondb@ccg.unam.mx"
        },
        "version": "1.0",
        "creationDate": f"{creation_date.strftime('%m-%d-%Y')}",
        "columnsDetails": "# Columns:\n# (1) sigmaName. Sigma Name\n# (2) regulatedName. TU regulated by the Sigma Factor (regulated TU)\n# (3) function. Regulatory Function of the Sigma on the regulated TU (+ activator, - repressor, -+ dual, ? unknown)\n# (4)promoterEvidences. Evidence that supports the regulation the sigma on the promoter and the promoter evidences\n# (5) confidenceLevel. RI confidence level based on its evidence (Values: Confirmed[C], Strong[S], Weak[W], Unknown[?])",
        "content": " \n".join(tfs_content),
        "rdbVersion": "12.0"
    }
    return tfs_doc