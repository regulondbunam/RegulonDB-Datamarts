import os
import json
import re
import multigenomic_api
import argparse
from datetime import datetime


def load_arguments_parser():
    parser = argparse.ArgumentParser(description="Extract and Transforms Datamarts from Multigenomic DB.",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter
                                     )
    parser.add_argument(
        "-db", "--database",
        help="Database where the multigenomic data is been stored",
        metavar="multigenomic",
        default="regulondbmultigenomic",
        required=False
    )

    parser.add_argument(
        "-u", "--url",
        help="URL to establish a connection between the process and MongoDB",
        metavar="mongodb://user:pass@localhost:27017/regulondbmultigenomic",
        default="mongodb://andresloal:15091052@localhost:27017/regulondbmultigenomic",
        required=False
    )

    parser.add_argument(
        "-f", "--file",
        help="file that contains Coexpression Matrix",
        metavar="coexpression-spearman-rank.txt",
        default="",
        required=True
    )

    parser.add_argument(
        "-o", "--output",
        help="folder were file gonna be writted",
        metavar="/build",
        default="../build",
        required=True
    )

    args = parser.parse_args()

    return args


def get_gene_list_dict(db_genes):
    gene_list_dict = []
    for gene in db_genes:
        gene_list_dict.append({
            "_id": gene.id,
            "name": gene.name,
            "locusTag": gene.bnumber
        })
    return gene_list_dict


def get_rgb_color(rank):
    if 0 <= rank < 44:
        return "5, 241, 59"
    elif 44 <= rank < 174:
        return "0, 218, 54"
    elif 174 <= rank < 390:
        return "0, 178, 45"
    elif 390 <= rank < 692:
        return "0, 101, 27"
    elif 692 <= rank < 1081:
        return "32, 64, 1"
    elif 1081 <= rank < 1556:
        return "103, 0, 0"
    elif 1556 <= rank < 2118:
        return "140, 0, 0"
    elif 2118 <= rank < 2766:
        return "178, 0, 0"
    elif 2766 <= rank < 3500:
        return "217, 2, 0"
    elif 3500 <= rank:
        return "255, 2, 0"


if __name__ == '__main__':
    arguments = load_arguments_parser()
    multigenomic_api.connect(arguments.database, arguments.url)
    rdb_genes = multigenomic_api.genes.get_all()

    fin = open(arguments.file, 'r')
    print(f"Working with {arguments.file}")
    pattern = re.compile("^b[0-9]{4}$")
    coexp_items = []
    a = []
    line_count = 1

    for line in fin.readlines():
        a.append([x for x in line.split(' ')])

    log = "Starting at: " + str(datetime.now())

    gene_list = get_gene_list_dict(rdb_genes)

    for line in a[1:]:
        count = line_count
        print(a[0][line_count])
        try:
            gene_1 = next(item for item in gene_list if item["locusTag"] == a[0][line_count])
        except:
            print(f"bnumber {a[count][0]} is not in collection")
            line_count += 1
            continue
        for i in line[line_count:]:
            if pattern.match(i):
                continue
            else:
                if count >= len(a):
                    break
                try:
                    gene_2 = next(item for item in gene_list if item["locusTag"] == a[count][0])
                except:
                    print(f"bnumber {a[count][0]} is not in collection")
                    count += 1
                    continue

                genes = [gene_1, gene_2]
                genes = sorted(genes, key=lambda d: d['locusTag'])
                score = float(i)
                item = {
                    "gene": genes,
                    "rank": score,
                    "rgbColor": get_rgb_color(score),
                    "organisms": {}
                }
                count += 1
                coexp_items.append(item)
        line_count += 1

    log = log + f" \nFile finished at: {str(datetime.now())}"
    
    final_json = {
        "geneCoexpressions": coexp_items
    }

    filename = os.path.join("", "coexpression_items")
    with open("{}.json".format(filename), 'w') as json_file:
        json.dump(final_json, json_file, indent=2, sort_keys=True)

    log = log + f" \nJSON File written at: {str(datetime.now())}"
    print(log)
