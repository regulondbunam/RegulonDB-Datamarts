def load_arguments():
    import argparse

    parser = argparse.ArgumentParser(description="Loads JSON files into the specified database",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     epilog="Either [Multigenomic or Datamarts] and [Directory or File] should be selected through arguments options"
                                     )

    parser.add_argument(
        "-db", "--database",
        help="Database where the multigenomic data is been stored",
        metavar="regulondbdatamarts",
        default="regulondbdatamarts",
        required=False
    )

    parser.add_argument(
        "-u", "--url",
        help="URL to establish a connection between the process and MongoDB",
        metavar="mongodb://user:pass@127.0.0.1:27017/regulondbdatamarts",
        default="mongodb://andresloal:15091052@127.0.0.1:27017/regulondbdatamarts",
        required=False
    )

    parser.add_argument(
        "-cd", "--collection_data",
        type=str,
        help="Input data files for regulondbdatamarts collections",
        metavar="lib/data",
        default="lib/data",
        required=True
    )

    parser.add_argument(
        "-s", "--schemas",
        type=str,
        help="jsonSchemas folder for regulondbdatamarts collections",
        metavar="lib/jsonschemas",
        default="lib/jsonschemas",
        required=True
    )

    arguments = parser.parse_args()
    return arguments

