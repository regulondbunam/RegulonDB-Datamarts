import multigenomic_api
from datetime import datetime

from mongoengine.errors import DoesNotExist


class Regulators:

    @property
    def objects(self):
        regulator_objects = get_all_regulators(multigenomic_api.regulatory_interactions.get_all())
        tf_objects = multigenomic_api.transcription_factors.get_all()
        for tf_obj in tf_objects:
            print(tf_obj.id)
            tf_obj["regulator_type"] = "transcriptionFactor"
            regulon_datamart = Regulators.RegulatorDatamart(tf_obj)
            yield regulon_datamart
        for reg_obj in regulator_objects:
            print(reg_obj.id)
            regulator = reg_obj
            if reg_obj.type == "srna":
                regulator = multigenomic_api.products.find_by_id(reg_obj.id)
            if reg_obj.type == "compound":
                regulator = multigenomic_api.regulatory_continuants.find_by_id(reg_obj.id)
            regulator["regulator_type"] = reg_obj.type
            regulon_datamart = Regulators.RegulatorDatamart(regulator)
            yield regulon_datamart
        del regulator_objects
        del tf_objects

    class RegulatorDatamart:
        def __init__(self, regulator):
            active_conformations = get_all_act_conf(regulator)
            self.regulator_name = regulator
            self.regulator = regulator
            self.regulator_synonyms = regulator.synonyms
            self.genes = regulator
            self.active_conf = active_conformations
            self.inactive_conf = regulator
            self.active_conf_syn = active_conformations
            self.inactive_conf_syn = regulator
            self.active_conf_effectors = active_conformations
            self.inactive_conf_effectors = regulator
            self.active_conf_effectors_syn = active_conformations
            self.inactive_conf_effectors_syn = regulator
            self.symmetry = regulator
            self.regulator_evidences = regulator.citations
            self.additive_evidences = regulator.additive_evidences_ids
            self.regulator_conf_pmids = active_conformations

        @property
        def regulator_name(self):
            return self._regulator_name

        @regulator_name.setter
        def regulator_name(self, regulator):
            self._regulator_name = ""
            if regulator['regulator_type'] == "compound":
                self._regulator_name = regulator.name
            else:
                self._regulator_name = regulator.abbreviated_name

        @property
        def regulator_synonyms(self):
            return self._regulator_synonyms

        @regulator_synonyms.setter
        def regulator_synonyms(self, synonyms):
            self._regulator_synonyms = None
            if len(synonyms) > 0:
                self._regulator_synonyms = ""
                for synonym in synonyms:
                    self._regulator_synonyms += f"{synonym};"
            if len(self._regulator_synonyms) > 0:
                self._regulator_synonyms = self._regulator_synonyms[:-1]

        @property
        def genes(self):
            return self._genes

        @genes.setter
        def genes(self, regulator):
            self._genes = ""
            if regulator.regulator_type == "transcriptionFactor":
                for product_id in regulator.products_ids:
                    product = multigenomic_api.products.find_by_id(product_id)
                    gene = multigenomic_api.genes.find_by_id(product.genes_id)
                    self._genes += f"{gene.name};"
            if regulator.regulator_type == "srna":
                gene = multigenomic_api.genes.find_by_id(regulator.genes_id)
                self._genes += f"{gene.name};"
            if len(self._genes) > 0:
                self._genes = self._genes[:-1]

        @property
        def active_conf(self):
            return self._active_conf

        @active_conf.setter
        def active_conf(self, active_confs):
            self._active_conf = ""
            if len(active_confs) > 0:
                for conf in active_confs:
                    conf_obj = {}
                    if conf['type'] == "product":
                        reg = multigenomic_api.products.find_by_id(conf["id"])
                        conf_obj = reg.abbreviated_name or reg.name
                    elif conf['type'] == "regulatoryComplex":
                        reg = multigenomic_api.regulatory_complexes.find_by_id(conf["id"])
                        conf_obj = reg.abbreviated_name or reg.name
                    elif conf['type'] == "regulatoryContinuant":
                        reg = multigenomic_api.regulatory_continuants.find_by_id(conf["id"])
                        conf_obj = reg.name
                    self._active_conf += f"{conf_obj};"
            if len(self._active_conf) > 0:
                self._active_conf = self._active_conf[:-1]

        @property
        def inactive_conf(self):
            return self._inactive_conf

        @inactive_conf.setter
        def inactive_conf(self, regulator):
            self._inactive_conf = ""
            if regulator.regulator_type == "transcriptionFactor":
                for conf in regulator.inactive_conformations:
                    if conf.type == "product":
                        conf = multigenomic_api.products.find_by_id(conf.id)
                    elif conf.type == "regulatoryComplex":
                        conf = multigenomic_api.regulatory_complexes.find_by_id(conf.id)
                    self._inactive_conf += f"{conf.abbreviated_name or conf.name};"
                if len(self._inactive_conf) > 0:
                    self._inactive_conf = self._inactive_conf[:-1]

        @property
        def active_conf_syn(self):
            return self._active_conf_syn

        @active_conf_syn.setter
        def active_conf_syn(self, active_confs):
            self._active_conf_syn = ""
            for conf in active_confs:
                for synonym in conf["synonyms"]:
                    self._active_conf_syn += f"{synonym};"
            if len(self._active_conf_syn) > 0:
                self._active_conf_syn = self._active_conf_syn[:-1]

        @property
        def inactive_conf_syn(self):
            return self._inactive_conf_syn

        @inactive_conf_syn.setter
        def inactive_conf_syn(self, regulator):
            self._inactive_conf_syn = ""
            if regulator.regulator_type == "transcriptionFactor":
                for conf in regulator.inactive_conformations:
                    if conf.type == "product":
                        conf = multigenomic_api.products.find_by_id(conf.id)
                    elif conf.type == "regulatoryComplex":
                        conf = multigenomic_api.regulatory_complexes.find_by_id(conf.id)
                    for synonym in conf.synonyms:
                        self._inactive_conf_syn += f"{synonym};"
                if len(self._inactive_conf_syn) > 0:
                    self._inactive_conf_syn = self._inactive_conf_syn[:-1]

        @property
        def active_conf_effectors(self):
            return self._active_conf_effectors

        @active_conf_effectors.setter
        def active_conf_effectors(self, active_confs):
            self._active_conf_effectors = ""
            for conf in active_confs:
                if conf["type"] == "regulatoryComplex" and conf["effectors"] != "":
                    self._active_conf_effectors += f"{conf['effectors']};"
            if len(self._active_conf_effectors) > 0:
                self._active_conf_effectors = self._active_conf_effectors[:-1]

        @property
        def inactive_conf_effectors(self):
            return self._inactive_conf_effectors

        @inactive_conf_effectors.setter
        def inactive_conf_effectors(self, regulator):
            self._inactive_conf_effectors = ""
            if regulator.regulator_type == "transcriptionFactor":
                for conf in regulator.inactive_conformations:
                    if conf.type == "regulatoryComplex":
                        conf = multigenomic_api.regulatory_complexes.find_by_id(conf.id)
                        for continuant_id in conf.regulatory_continuants_ids:
                            continuant = multigenomic_api.regulatory_continuants.find_by_id(continuant_id)
                            self._inactive_conf_effectors += f"{continuant.name};"
                if len(self._inactive_conf_effectors) > 0:
                    self._inactive_conf_effectors = self._inactive_conf_effectors[:-1]

        @property
        def active_conf_effectors_syn(self):
            return self._active_conf_effectors_syn

        @active_conf_effectors_syn.setter
        def active_conf_effectors_syn(self, active_confs):
            self._active_conf_effectors_syn = ""
            for conf in active_confs:
                if conf["type"] == "regulatoryComplex" and conf["effectorsSynonyms"] != "":
                    self._active_conf_effectors_syn += f"{conf['effectorsSynonyms']}"
            if len(self._active_conf_effectors_syn) > 0:
                self._active_conf_effectors_syn = self._active_conf_effectors_syn[:-1]

        @property
        def inactive_conf_effectors_syn(self):
            return self._inactive_conf_effectors_syn

        @inactive_conf_effectors_syn.setter
        def inactive_conf_effectors_syn(self, regulator):
            self._inactive_conf_effectors_syn = ""
            if regulator.regulator_type == "transcriptionFactor":
                for conf in regulator.inactive_conformations:
                    if conf.type == "regulatoryComplex":
                        conf = multigenomic_api.regulatory_complexes.find_by_id(conf.id)
                        for continuant_id in conf.regulatory_continuants_ids:
                            continuant = multigenomic_api.regulatory_continuants.find_by_id(continuant_id)
                            for synonym in continuant.synonyms:
                                self._inactive_conf_effectors_syn += f"{synonym};"
                if len(self._inactive_conf_effectors_syn) > 0:
                    self._inactive_conf_effectors_syn = self._inactive_conf_effectors_syn[:-1]

        @property
        def symmetry(self):
            return self._symmetry

        @symmetry.setter
        def symmetry(self, regulator):
            self._symmetry = ""
            if regulator.regulator_type == "transcriptionFactor":
                if regulator.symmetry:
                    for symm_ele in regulator.symmetry:
                        self._symmetry += f"{symm_ele};"
                    if len(self._symmetry) > 0:
                        self._symmetry = self._symmetry[:-1]

        @property
        def regulator_evidences(self):
            return self._regulator_evidences

        @regulator_evidences.setter
        def regulator_evidences(self, citations):
            self._regulator_evidences = []
            for citation in citations:
                if citation.evidences_id:
                    citation_dict = multigenomic_api.evidences.find_by_id(citation.evidences_id)
                    citation_item = f"[{citation_dict.code}:{citation_dict.type}]"
                    if citation_item not in self._regulator_evidences:
                        self._regulator_evidences.append(citation_item)
            self._regulator_evidences = "".join(self._regulator_evidences)

        @property
        def additive_evidences(self):
            return self._additive_evidences

        @additive_evidences.setter
        def additive_evidences(self, additive_evs_ids):
            self._additive_evidences = ""
            if len(additive_evs_ids) > 0:
                self._additive_evidences = ""
                for additive_evs_id in additive_evs_ids:
                    additive_evidence_dict = multigenomic_api.additive_evidences.find_by_id(additive_evs_id)
                    self._additive_evidences += f"[{additive_evidence_dict.code}:{additive_evidence_dict.confidence_level}]"

        @property
        def regulator_conf_pmids(self):
            return self._regulator_conf_pmids

        @regulator_conf_pmids.setter
        def regulator_conf_pmids(self, active_conformations):
            self._regulator_conf_pmids = []
            for act_conf in active_conformations:
                if act_conf["citations"]:
                    for citation in act_conf["citations"]:
                        if citation.publications_id:
                            publication = multigenomic_api.publications.find_by_id(citation.publications_id)
                            if f"{publication.pmid}; " not in self._regulator_conf_pmids:
                                self._regulator_conf_pmids.append(f"{publication.pmid}; ")
            if self.regulator.regulator_type == "transcriptionFactor":
                if self.regulator.inactive_conformations:
                    if len(self.regulator.inactive_conformations) > 0:
                        for conf in self.regulator.inactive_conformations:
                            if conf.type == "product":
                                conf = multigenomic_api.products.find_by_id(conf.id)
                            elif conf.type == "regulatoryComplex":
                                conf = multigenomic_api.regulatory_complexes.find_by_id(conf.id)
                            for citation in conf["citations"]:
                                if citation.publications_id:
                                    publication = multigenomic_api.publications.find_by_id(citation.publications_id)
                                    if f"{publication.pmid}; "not in self._regulator_conf_pmids:
                                        self._regulator_conf_pmids.append(f"{publication.pmid}; ")
            self._regulator_conf_pmids = "".join(self._regulator_conf_pmids)[:-2]

        def to_row(self):
            # TODO: add pmids
            return f"{self.regulator.id}" \
                   f"\t{self.regulator_name}" \
                   f"\t{self.regulator_synonyms}" \
                   f"\t{self.regulator.regulator_type}" \
                   f"\t{self.genes}" \
                   f"\t{self.active_conf}" \
                   f"\t{self.inactive_conf}" \
                   f"\t{self.active_conf_syn}" \
                   f"\t{self.inactive_conf_syn}" \
                   f"\t{self.active_conf_effectors}" \
                   f"\t{self.inactive_conf_effectors}" \
                   f"\t{self.active_conf_effectors_syn}" \
                   f"\t{self.inactive_conf_effectors_syn}" \
                   f"\t{self.symmetry}" \
                   f"\t{self.regulator_evidences}" \
                   f"\t{self.additive_evidences}" \
                   f"\t{self.regulator.confidence_level}" \
                   f"\t{self.regulator_conf_pmids}"


def get_all_regulators(reg_ints):
    regulators = []
    for ri in reg_ints:
        if ri.regulator.type == "product":
            product = multigenomic_api.products.find_by_id(ri.regulator.id)
            if product.type == "small RNA":
                regulator = ri.regulator
                regulator["type"] = "srna"
                if ri.regulator not in regulators:
                    regulators.append(ri.regulator)
        if ri.regulator.type == "regulatoryContinuant":
            regulator = ri.regulator
            regulator["type"] = "compound"
            if ri.regulator not in regulators:
                regulators.append(ri.regulator)
    return regulators


def get_all_act_conf(regulator):
    conformations = []
    reg_complex = None
    try:
        reg_complex = multigenomic_api.regulatory_complexes.find_by_name(regulator.name)
    except DoesNotExist:
        pass
    if regulator.regulator_type == "transcriptionFactor":
        if reg_complex is not None:
            reg_ints = multigenomic_api.regulatory_interactions.find_by_regulator_id(
                reg_complex.id)
            if len(reg_ints) > 0:
                effectors = ""
                effectors_syn = ""
                for reg_cont_id in reg_complex.regulatory_continuants_ids:
                    continuant = multigenomic_api.regulatory_continuants.find_by_id(reg_cont_id)
                    effectors += f"{continuant.name};"
                    for synonym in continuant.synonyms:
                        effectors_syn += f"{synonym};"
                conformations.append({"name": regulator.name,
                                      "id": reg_complex.id,
                                      "synonyms": regulator.synonyms,
                                      "type": "regulatoryComplex",
                                      "effectors": effectors,
                                      "effectorsSynonyms": effectors_syn,
                                      "citations": regulator.citations})
        else:
            for product_id in regulator.products_ids:
                reg_ints = multigenomic_api.regulatory_interactions.find_by_regulator_id(product_id)
                if len(reg_ints) > 0:
                    prod = multigenomic_api.products.find_by_id(product_id)
                    conformations.append({"name": prod.name,
                                          "id": prod.id,
                                          "synonyms": prod.synonyms,
                                          "type": "product",
                                          "citations": prod.citations})
            conformation_object = None
            for conformation in regulator.active_conformations:
                if conformation.type == "product":
                    prod = multigenomic_api.products.find_by_id(conformation.id)
                    conformation_object = {"name": prod.name,
                                           "id": prod.id,
                                           "synonyms": prod.synonyms,
                                           "type": "product",
                                           "citations": prod.citations}
                elif conformation.type == "regulatoryComplex":
                    complx = multigenomic_api.regulatory_complexes.find_by_id(conformation.id)
                    effectors = ""
                    effectors_syn = ""
                    for reg_cont_id in complx.regulatory_continuants_ids:
                        continuant = multigenomic_api.regulatory_continuants.find_by_id(reg_cont_id)
                        effectors += continuant.name
                        for synonym in continuant.synonyms:
                            effectors_syn += f"{synonym};"
                    conformation_object = {"name": complx.name,
                                           "id": complx.id,
                                           "synonyms": complx.synonyms,
                                           "type": "regulatoryComplex",
                                           "effectors": effectors,
                                           "effectorsSynonyms": effectors_syn,
                                           "citations": complx.citations}
                conformations.append(conformation_object)
    return conformations


def all_regulators_rows():
    regulators = Regulators()
    regulators_content = ["1)regulatorId\t2)regulatorName\t3)regulatorSynonyms\t4)regulatorType\t5)geneCodingForRegulator\t6)regulatorActiveConformations\t7)regulatorInactiveConformations\t8)regulatorActiveConformationsSynonyms\t9)regulatorInactiveConformationsSynonyms\t10)regulatorActiveConformationsEffector\t11)regulatorInactiveConformationsEffector\t12)regulatorActiveConformationsEffectorSynonyms\t13)regulatorInactiveConformationsEffectorSynonyms\t14)symmetry\t15)regulatorEvidences\t16)additiveEvidences\t17)confidenceLevel\t18)regulatorConformationPMID"]
    for regulator in regulators.objects:
        regulators_content.append(regulator.to_row())
    creation_date = datetime.now()
    regulators_doc = {
        "_id": "RDBECOLIDLF00009",
        "fileName": "RegulatorSet",
        "title": "Complete Regulators Set",
        "fileFormat": "rif-version 1",
        "license": "RegulonDB is free for academic/noncommercial use\t\tUser is not entitled to change or erase data sets of the RegulonDB\tdatabase or to eliminate copyright notices from RegulonDB. Furthermore,\tUser is not entitled to expand RegulonDB or to integrate RegulonDB partly\tor as a whole into other databank systems, without prior written consent\regulatorrom CCG-UNAM.\t\tPlease check the license at http://regulondb.ccg.unam.mx/menu/download/full_version/terms_and_conditions.jsp",
        "citation": "Tierrafría, V. H. et al. (2022). RegulonDB 11.0: Comprehensive high-throughput datasets on transcriptional regulation in Escherichia coli K-12,\tMicrob Genom. 2022 May;8(5). doi: 10.1099/mgen.0.000833. PMID: 35584008. https://doi.org/10.1099/mgen.0.000833",
        "contact": {
            "person": "RegulonDB Team",
            "webPage": "http://regulondb.ccg.unam.mx/menu/about_regulondb/contact_us/index.jsp",
            "email": "regulondb@ccg.unam.mx"
        },
        "version": "",
        "creationDate": f"{creation_date.strftime('%m-%d-%Y')}",
        "columnsDetails": "Columns:\n(1) Regulator identifier assigned by RegulonDB\n(2) Regulator Name\n(3) Regulator Synonyms List\n(4) Regulator Synonyms List\n(5) Gene Coding for the Regulator\n(6) Regulator Active Conformations\n(7) Regulator Inactive Conformations\n(8) Regulator Active Conformations  Synonyms List\n(9) Regulator Inactive Conformations  Synonyms List\n(10) Effector Name related to  Regulator Active Conformations\n(11) Effector Name related to  Regulator Inactive Conformations\n(12) Effector Synonyms List related to Regulator Active Conformations Regulator \n(13) Effector Synonyms List related to Regulator Inactive Conformations Regulator\n(14) Regulator Symmetry\n(15) Evidence that supports the Regulator conformation [Evidence code | Evidence type: C = Confirmed, S = Strong, W = Weak | Evidence name ]\n(16) addEvidence. Additive Evidence [CV(EvidenceCode1/EvidenceCodeN)|Confidence Level]\n(17) confidenceLevel. Confidence level (Values: Confirmed, Strong, Weak)\n(18) Regulator conformation reference identifier (PMID)",
        "content": " \n".join(regulators_content)
    }
    return regulators_doc
