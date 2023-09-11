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


def get_all_downloadable_docs():
    downloadable_files_dm = []
    # Regulatory Interactions
    reg_ints = regulatoryInteractions_file.all_ris_rows()
    downloadable_files_dm.append(reg_ints)
    # Promoters
    promoters = promoters_file.all_promoters_rows()
    downloadable_files_dm.append(promoters)
    # TranscriptionUnits
    tus = transcriptionUnits_file.all_tus_rows()
    downloadable_files_dm.append(tus)
    # Operons
    operons = operons_file.all_operons_rows()
    downloadable_files_dm.append(operons)
    # TranscriptionFactors
    tfs = transcriptionFactors_file.all_tfs_rows()
    downloadable_files_dm.append(tfs)
    # Genes
    genes = gene_product_file.all_gene_rows()
    downloadable_files_dm.append(genes)
    # Terminators
    terminators = terminators_file.all_terminators_rows()
    downloadable_files_dm.append(terminators)
    # GeneSequence
    geneSequences = gene_sequence_file.all_gene_rows()
    downloadable_files_dm.append(geneSequences)
    # Regulators
    regulators = regulators_file.all_regulators_rows()
    downloadable_files_dm.append(regulators)
    # TF-Gene
    tf_gene = tfgene_file.get_all_rows()
    downloadable_files_dm.append(tf_gene)

    return downloadable_files_dm
