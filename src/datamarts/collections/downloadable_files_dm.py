from src.datamarts.domain.downloadableFiles_datamart import regulatoryInteractions_file
from src.datamarts.domain.downloadableFiles_datamart import promoters_file
from src.datamarts.domain.downloadableFiles_datamart import transcriptionUnits_file
from src.datamarts.domain.downloadableFiles_datamart import operons_file
from src.datamarts.domain.downloadableFiles_datamart import transcriptionFactors_file
from src.datamarts.domain.downloadableFiles_datamart import gene_product_file
from src.datamarts.domain.downloadableFiles_datamart import terminators_file
from src.datamarts.domain.downloadableFiles_datamart import gene_sequence_file
from src.datamarts.domain.downloadableFiles_datamart import regulators_file
from src.datamarts.domain.downloadableFiles_datamart import tfgene_file
from src.datamarts.domain.downloadableFiles_datamart import tfgene_file_release4
from src.datamarts.domain.downloadableFiles_datamart import regulatoryInteractions_tf_file
from src.datamarts.domain.downloadableFiles_datamart import object_evidences_file
from src.datamarts.domain.downloadableFiles_datamart import additiveEvidences_file
from src.datamarts.domain.downloadableFiles_datamart import tftu_file
from src.datamarts.domain.downloadableFiles_datamart import tf_tf_file
from src.datamarts.domain.downloadableFiles_datamart import gene_product_ids_file
from src.datamarts.domain.downloadableFiles_datamart import sigma_gene_file
from src.datamarts.domain.downloadableFiles_datamart import sigma_tu_file
from src.datamarts.domain.downloadableFiles_datamart import utr_5_3_sequence_file
from src.datamarts.domain.downloadableFiles_datamart import all_pmids_in_objects
from src.datamarts.domain.downloadableFiles_datamart import tu_prom_operon_tf_file
from src.datamarts.domain.downloadableFiles_datamart import hth_genes_motif
from src.datamarts.domain.downloadableFiles_datamart import network_conformation_gene


from src.datamarts.domain.general.print_progress import print_progress


def get_all_downloadable_docs(rdb_version, citation, url):
    rdb_files = [
        ("RiSet",                        regulatoryInteractions_file.all_ris_rows),
        ("TF-RiSet",                     regulatoryInteractions_tf_file.all_ris_rows),
        ("RegulatorSet",                 regulators_file.all_regulators_rows),
        ("TFSet",                        transcriptionFactors_file.all_tfs_rows),
        ("Regulator-Gene",               tfgene_file.get_all_rows),
        ("Regulator-Gene-internal",      tfgene_file_release4.get_all_rows),
        ("Gene-Prod",                    gene_product_file.all_gene_rows),
        ("GeneSeq",                      gene_sequence_file.all_gene_rows),
        ("Operons",                      operons_file.all_operons_rows),
        ("TUs",                          transcriptionUnits_file.all_tus_rows),
        ("Promoters",                    promoters_file.all_promoters_rows),
        ("Terminators",                  terminators_file.all_terminators_rows),
        ("Evidences",                    object_evidences_file.all_evidences_rows),
        ("Additive Evidences",           additiveEvidences_file.all_evidences_rows),
        ("Regulator-TU",                 tftu_file.get_all_rows),
        ("Regulator-Regulator",          tf_tf_file.get_all_rows),
        ("Gene-Prod IDs",                gene_product_ids_file.all_gene_rows),
        ("Sigma-Gene",                   sigma_gene_file.get_all_rows),
        ("Sigma-TU",                     sigma_tu_file.get_all_rows),
        ("UTR Sequences",                utr_5_3_sequence_file.all_utr_rows),
        ("PMIDS",                        "PMID_HANDLER"),  # caso especial
        ("Tu-Prom-Operon-TF",            tu_prom_operon_tf_file.all_tus_rows),
        ("HTH Genes",                    hth_genes_motif.all_gene_rows),
        ("Conformation Gene",            network_conformation_gene.get_all_rows),
    ]

    downloadable_files_dm = []

    # Inicializar PMIDS mongo antes del conteo
    all_pmids_in_objects.init_mongo(url)

    total_tasks = len(rdb_files)
    current = 0

    # Procesar cada archivo
    for label, handler in rdb_files:
        current += 1
        print_progress(current, total_tasks, "Downloadable Files")

        if handler == "PMID_HANDLER":
            # caso especial PMIDS
            rows = all_pmids_in_objects.generate_pmids_doc(rdb_version, citation)
        else:
            rows = handler(rdb_version, citation)

        downloadable_files_dm.append(rows)

    return downloadable_files_dm
