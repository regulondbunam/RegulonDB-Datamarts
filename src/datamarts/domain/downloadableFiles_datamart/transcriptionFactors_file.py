import multigenomic_api
from datetime import datetime

from mongoengine.errors import DoesNotExist


class TranscriptionFactor:

    @property
    def objects(self):
        tf_objects = multigenomic_api.transcription_factors.get_all()
        for tf_object in tf_objects:
            # print(tf_object.id)
            ri_row = TranscriptionFactor.TFDatamart(tf_object)
            yield ri_row
        del tf_objects

    class TFDatamart:
        def __init__(self, tf):
            active_conformations = get_all_act_conf(tf)
            self.tf = tf
            self.tf_synonyms = tf.synonyms
            self.genes = tf.products_ids
            self.active_conf = active_conformations
            self.inactive_conf = tf.inactive_conformations
            self.active_conf_syn = active_conformations
            self.inactive_conf_syn = tf.inactive_conformations
            self.active_conf_effectors = active_conformations
            self.inactive_conf_effectors = tf.inactive_conformations
            self.active_conf_effectors_syn = active_conformations
            self.inactive_conf_effectors_syn = tf.inactive_conformations
            self.symmetry = tf.symmetry
            self.tf_evidences = tf.citations
            self.additive_evidences = tf.additive_evidences_ids
            self.tf_conf_pmids = active_conformations

        @property
        def tf_synonyms(self):
            return self._tf_synonyms

        @tf_synonyms.setter
        def tf_synonyms(self, synonyms):
            self._tf_synonyms = ""
            for synonym in synonyms:
                self._tf_synonyms += f"{synonym};"
            if len(self._tf_synonyms) > 0:
                self._tf_synonyms = self._tf_synonyms[:-1]

        @property
        def genes(self):
            return self._genes

        @genes.setter
        def genes(self, products_ids):
            self._genes = ""
            for product_id in products_ids:
                product = multigenomic_api.products.find_by_id(product_id)
                gene = multigenomic_api.genes.find_by_id(product.genes_id)
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
        def inactive_conf(self, inactive_confs):
            self._inactive_conf = ""
            for conf in inactive_confs:
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
        def inactive_conf_syn(self, inactive_confs):
            self._inactive_conf_syn = ""
            for conf in inactive_confs:
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
        def inactive_conf_effectors(self, inactive_confs):
            self._inactive_conf_effectors = ""
            for conf in inactive_confs:
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
        def inactive_conf_effectors_syn(self, inactive_confs):
            self._inactive_conf_effectors_syn = ""
            for conf in inactive_confs:
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
        def symmetry(self, symmetry):
            self._symmetry = ""
            if symmetry:
                for symm_ele in symmetry:
                    self._symmetry += f"{symm_ele};"
                if len(self._symmetry) > 0:
                    self._symmetry = self._symmetry[:-1]

        @property
        def tf_evidences(self):
            return self._tf_evidences

        @tf_evidences.setter
        def tf_evidences(self, citations):
            self._tf_evidences = []
            for citation in citations:
                if citation.evidences_id:
                    citation_dict = multigenomic_api.evidences.find_by_id(citation.evidences_id)
                    citation_item = f"[{citation_dict.code}:{citation_dict.type}]"
                    if citation_item not in self._tf_evidences:
                        self._tf_evidences.append(citation_item)
            self._tf_evidences = "".join(self._tf_evidences)

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
        def tf_conf_pmids(self):
            return self._tf_conf_pmids

        @tf_conf_pmids.setter
        def tf_conf_pmids(self, active_conformations):
            self._tf_conf_pmids = []
            for act_conf in active_conformations:
                if act_conf["citations"]:
                    for citation in act_conf["citations"]:
                        if citation.publications_id:
                            publication = multigenomic_api.publications.find_by_id(citation.publications_id)
                            if f"{publication.pmid}; " not in self._tf_conf_pmids:
                                self._tf_conf_pmids.append(f"{publication.pmid}; ")
            if self.tf.inactive_conformations:
                if len(self.tf.inactive_conformations) > 0:
                    for conf in self.tf.inactive_conformations:
                        if conf.type == "product":
                            conf = multigenomic_api.products.find_by_id(conf.id)
                        elif conf.type == "regulatoryComplex":
                            conf = multigenomic_api.regulatory_complexes.find_by_id(conf.id)
                        for citation in conf["citations"]:
                            if citation.publications_id:
                                publication = multigenomic_api.publications.find_by_id(citation.publications_id)
                                if f"{publication.pmid}; "not in self._tf_conf_pmids:
                                    self._tf_conf_pmids.append(f"{publication.pmid}; ")
            self._tf_conf_pmids = "".join(self._tf_conf_pmids)[:-2]

        def to_row(self):
            return f"{self.tf.id}" \
                   f"\t{self.tf.abbreviated_name}" \
                   f"\t{self.tf_synonyms}" \
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
                   f"\t{self.tf_evidences}" \
                   f"\t{self.additive_evidences}" \
                   f"\t{self.tf.confidence_level}" \
                   f"\t{self.tf_conf_pmids}"


def get_all_act_conf(tf):
    conformations = []
    reg_complex = None
    try:
        reg_complex = multigenomic_api.regulatory_complexes.find_by_name(tf.name)
    except DoesNotExist:
        pass
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
            conformations.append({"name": tf.name,
                                  "id": reg_complex.id,
                                  "synonyms": tf.synonyms,
                                  "type": "regulatoryComplex",
                                  "effectors": effectors,
                                  "effectorsSynonyms": effectors_syn,
                                  "citations": tf.citations})
    else:
        for product_id in tf.products_ids:
            reg_ints = multigenomic_api.regulatory_interactions.find_by_regulator_id(product_id)
            if len(reg_ints) > 0:
                prod = multigenomic_api.products.find_by_id(product_id)
                conformations.append({"name": prod.name,
                                      "id": prod.id,
                                      "synonyms": prod.synonyms,
                                      "type": "product",
                                      "citations": prod.citations})
        conformation_object = None
        for conformation in tf.active_conformations:
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


def all_tfs_rows():
    trans_factors = TranscriptionFactor()
    tfs_content = ["1)tfId\t2)tfName\t3)tfSynonyms\t4)geneCodingForTF\t5)tfActiveConformations\t6)tfInactiveConformations\t7)tfActiveConformationsSynonyms\t8)tfInactiveConformationsSynonyms\t9)tfActiveConformationsEffector\t10)tfInactiveConformationsEffector\t11)tfActiveConformationsEffectorSynonyms\t12)tfInactiveConformationsEffectorSynonyms\t13)symmetry\t14)tfEvidences\t15)additiveEvidences\t16)confidenceLevel\t17)tfConformationPMID"]
    for tf in trans_factors.objects:
        tfs_content.append(tf.to_row())
    creation_date = datetime.now()
    tfs_doc = {
        "_id": "RDBECOLIDLF00004",
        "fileName": "TFSet",
        "title": "Complete Transcription Factor Set",
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
        "columnsDetails": "Columns:\n(1) Transcription Factor (TF) identifier assigned by RegulonDB\n(2) TF Name\n(3) TF Synonyms List\n(4) Gene Coding for the TF\n(5) TF Active Conformations\n(6) TF Inactive Conformations\n(7) TF Active Conformations  Synonyms List\n(8) TF Inactive Conformations  Synonyms List\n(9) Effector Name related to  TF Active Conformations\n(10) Effector Name related to  TF Inactive Conformations\n(11) Effector Synonyms List related to TF Active Conformations TF \n(12) Effector Synonyms List related to TF Inactive Conformations TF\n(13) TF Symmetry\n(14) Evidence that supports the TF conformation [Evidence code | Evidence type: C = Confirmed, S = Strong, W = Weak | Evidence name ]\n(15) addEvidence. Additive Evidence [CV(EvidenceCode1/EvidenceCodeN)|Confidence Level]\n(16) confidenceLevel. Confidence level (Values: Confirmed, Strong, Weak)\n(17) TF conformation reference identifier (PMID)",
        "content": " \n".join(tfs_content),
        "rdbVersion": "12.0"
    }
    return tfs_doc
