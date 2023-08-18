from src.datamarts.domain.downloadableFiles_datamart import regulatoryInteractions_file


def get_all_downloadable_docs():
    downloadable_files_dm = []

    reg_ints = regulatoryInteractions_file.all_ris_rows()
    downloadable_files_dm.append(reg_ints)

    return downloadable_files_dm
