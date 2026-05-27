import multigenomic_api
from datetime import datetime

from mongoengine.errors import DoesNotExist


class TFConformation:

    @property
    def objects(self):
        tf_objects = multigenomic_api.transcription_factors.get_all()
        for tf_object in tf_objects:
            # print(tf_object.id)
            ri_row = TFConformation.TFDatamart(tf_object)
            yield ri_row
        del tf_objects

    class TFDatamart:
        def __init__(self, tf):
            self.tf = tf

        def to_row(self) -> str:
            rows: list[str] = []
            tf_name = self.tf.abbreviated_name
            def build_complex_row(complex_obj) -> str | None:
                products = []
                stoichiometry = 0
                cplx_type = "complex-products"
                if complex_obj.regulatory_continuants_ids:
                    cplx_type = "complex-continuants"
                for prod in complex_obj.products:
                    coefficient = prod.coefficient or ""
                    product = multigenomic_api.products.find_by_id(prod.products_id)
                    products.append(f"{product.abbreviated_name}:{coefficient}")
                    if coefficient:
                        stoichiometry += coefficient
                if not products:
                    return None
                return (
                    f"{tf_name}\t"
                    f"{complex_obj.abbreviated_name or complex_obj.name}\t"
                    f"{cplx_type}\t"
                    f"{';'.join(products)}\t"
                    f"{'' if stoichiometry == 0 else stoichiometry}\t"
                    f"{define_oligomer_class(len(products))}"
                )
            # --- Regulatory complex by TF name ---
            active_reg_complex_ids = {
                act_conf.id
                for act_conf in self.tf.active_conformations
                if act_conf.type == "regulatoryComplex"
            }

            try:
                reg_cplx = multigenomic_api.regulatory_complexes.find_by_name(self.tf.name)
            except DoesNotExist:
                reg_cplx = None
            if reg_cplx and reg_cplx.id not in active_reg_complex_ids:
                reg_ints = multigenomic_api.regulatory_interactions.find_by_regulator_id(reg_cplx.id)
                if len(reg_ints) > 0:
                    row = build_complex_row(reg_cplx)
                    if row:
                        rows.append(row)

            # --- Product by TF name ---
            active_product_ids = {
                act_conf.id
                for act_conf in self.tf.active_conformations
                if act_conf.type == "product"
            }

            try:
                prod = multigenomic_api.products.find_by_name(self.tf.name)
            except DoesNotExist:
                prod = None

            if prod and prod.id not in active_product_ids:
                reg_ints = multigenomic_api.regulatory_interactions.find_by_regulator_id(prod.id)
                if len(reg_ints) > 0:
                    product_name = prod.abbreviated_name or prod.name
                    rows.append(
                        f"{tf_name}\t"
                        f"{product_name}\t"
                        f"product\t"
                        f"{product_name}:\t"
                        f"\t"
                        f"Unknown"
                    )

            # --- Active conformations ---
            for act_conf in self.tf.active_conformations:
                if act_conf.type == "product":
                    product = multigenomic_api.products.find_by_id(act_conf.id)
                    product_name = product.abbreviated_name or product.name
                    rows.append(
                        f"{tf_name}\t"
                        f"{product_name}\t"
                        f"{act_conf.type}\t"
                        f"{product_name}:\t"
                        f"\t"
                        f"Unknown"
                    )
                elif act_conf.type == "regulatoryComplex":
                    regulatory_complex = multigenomic_api.regulatory_complexes.find_by_id(act_conf.id)
                    row = build_complex_row(regulatory_complex)
                    if row:
                        rows.append(row)
            return "\n".join(rows)


def define_oligomer_class(num_products: int) -> str:
    component_class = ""
    if num_products > 0:
        # Homo vs Hetero
        component_class = "Heterotípico " if num_products > 1 else "Homotípico"
    else:
        component_class = "Unknown"
    return f"{component_class}"



def all_tfs_rows(rdb_version, citation):
    trans_factors = TFConformation()
    tfs_content = ["1)tfName\t2)ActiveConformation\t3)ActiveConfType\t4)Components\t5)Stoichiometry\t6)OligomerClass"]
    for tf in trans_factors.objects:
        #print(repr(tf.to_row()))
        tfs_content.append(tf.to_row())
    creation_date = datetime.now()
    tfs_doc = {
        "_id": "RDBECOLIDLF00030",
        "fileName": "TF-ConformationSet",
        "title": "Complete Transcription Factor with Conformations Set",
        "fileFormat": "rif-version 1",
        "license": "# RegulonDB is free for academic/noncommercial use\n# User is not entitled to change or erase data sets of the RegulonDB\n# database or to eliminate copyright notices from RegulonDB. Furthermore,\n# User is not entitled to expand RegulonDB or to integrate RegulonDB partly\n# or as a whole into other databank systems, without prior written consent\n# from CCG-UNAM.\n# Please check the license at https://regulondb.ccg.unam.mx/manual/aboutUs/terms-conditions",
        "citation": citation,
        "contact": {
            "person": "RegulonDB Team",
            "webPage": None,
            "email": "regulondb@ccg.unam.mx"
        },
        "version": "1.0",
        "creationDate": f"{creation_date.strftime('%m-%d-%Y')}",
        "columnsDetails": "# Columns:\n"
                            "# (1) Canonical name of the transcription factor\n"
                            "# (2) Name of the active conformation of the transcription factor (This may include the effector, for example: AraC–arabinose)\n"
                            "# (3) Type of the active conformation of the transcription factor\n"
                            "# (4) List of components that participate in the active conformation, formatted as: component:coefficient (Components are separated by semicolons ';')\n"
                            "# (5) Total number of components participating in the active conformation, calculated as the sum of all coefficients.\n"
                            "# (6) Classification of the active conformation based on its composition",
        "content": "\n".join(sorted(tfs_content)),
        "rdbVersion": rdb_version,
        "description": "Transcription factors and their conformations (subset of RegulatorSet).",
        "group": "REGULATOR"
    }
    return tfs_doc
