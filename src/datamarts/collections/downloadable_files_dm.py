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


def get_all_downloadable_docs(rdb_version, citation):
    downloadable_files_dm = []
    # Regulatory Interactions
    print("RiSet")
    reg_ints = regulatoryInteractions_file.all_ris_rows(rdb_version, citation)
    downloadable_files_dm.append(reg_ints)
    # Regulatory Interactions of TFs
    print("TF-RiSet")
    tf_reg_ints = regulatoryInteractions_tf_file.all_ris_rows(rdb_version, citation)
    downloadable_files_dm.append(tf_reg_ints)
    # Regulators
    print("RegulatorSet")
    regulators = regulators_file.all_regulators_rows(rdb_version, citation)
    downloadable_files_dm.append(regulators)
    # TranscriptionFactors
    print("TFSet")
    tfs = transcriptionFactors_file.all_tfs_rows(rdb_version, citation)
    downloadable_files_dm.append(tfs)
    # NetworkRegulatorGene
    print("Regulator-Gene")
    tf_gene = tfgene_file.get_all_rows(rdb_version, citation)
    downloadable_files_dm.append(tf_gene)
    # NetworkRegulatorGene_internal
    print("Regulator-Gene-internal")
    tf_gene_r4 = tfgene_file_release4.get_all_rows(rdb_version, citation)
    downloadable_files_dm.append(tf_gene_r4)
    # Gene-Product
    print("Gene-Prod")
    genes = gene_product_file.all_gene_rows(rdb_version, citation)
    downloadable_files_dm.append(genes)
    # GeneSequence
    print("GeneSeq")
    geneSequences = gene_sequence_file.all_gene_rows(rdb_version, citation)
    downloadable_files_dm.append(geneSequences)
    # Operons
    print("Operons")
    operons = operons_file.all_operons_rows(rdb_version, citation)
    downloadable_files_dm.append(operons)
    # TranscriptionUnits
    print("TUs")
    tus = transcriptionUnits_file.all_tus_rows(rdb_version, citation)
    downloadable_files_dm.append(tus)
    # Promoters
    print("Promoters")
    promoters = promoters_file.all_promoters_rows(rdb_version, citation)
    downloadable_files_dm.append(promoters)
    # Terminators
    print("Terminators")
    terminators = terminators_file.all_terminators_rows(rdb_version, citation)
    downloadable_files_dm.append(terminators)
    # Object Evidences
    print("Evidences")
    obj_evs = object_evidences_file.all_evidences_rows(rdb_version, citation)
    downloadable_files_dm.append(obj_evs)
    # Additive Evidences
    print("Additive Evidences")
    add_evs = additiveEvidences_file.all_evidences_rows(rdb_version, citation)
    downloadable_files_dm.append(add_evs)
    # NetworkRegulatorGene
    print("Regulator-TU")
    tf_tu = tftu_file.get_all_rows(rdb_version, citation)
    downloadable_files_dm.append(tf_tu)
    # NetworkRegulator-Regulator
    print("Regulator-Regulator")
    tf_tf = tf_tf_file.get_all_rows(rdb_version, citation)
    downloadable_files_dm.append(tf_tf)
    # Gene-Product Ids
    print("Gene-Prod")
    gene_prod_ids = gene_product_ids_file.all_gene_rows(rdb_version, citation)
    downloadable_files_dm.append(gene_prod_ids)
    # NetworkSigma-Gene
    print("Sigma-Gene")
    sigma_gene = sigma_gene_file.get_all_rows(rdb_version, citation)
    downloadable_files_dm.append(sigma_gene)
    # NetworkSigma-TU
    print("Sigma-TU")
    sigma_tu = sigma_tu_file.get_all_rows(rdb_version, citation)
    downloadable_files_dm.append(sigma_tu)
    # NetworkSigma-TU
    print("UTR Sequences")
    utr_seq = utr_5_3_sequence_file.all_utr_rows(rdb_version, citation)
    downloadable_files_dm.append(utr_seq)

    return downloadable_files_dm
