import os
import re
import pandas as pd
from my_functions.genomic_functions import extract_denovo_info

from my_functions.paths import GENOMES_LIST, DENSE_DIR
clusters_file = "/home/eliott.tempez/Documents/M2_Stage_I2BC/results/18_clustering/good_candidates_clustering.tsv"
out_folder = "/home/eliott.tempez/Documents/M2_Stage_I2BC/results/18_clustering/"






if __name__ == "__main__":
    # Read list of genomes
    with open(GENOMES_LIST, "r") as f:
        genomes = f.readline().split()
    genomes = [re.sub('"', '', g) for g in genomes]

    # Extract denovo info
    denovo_genes = {}
    for genome in genomes:
        local_dict = extract_denovo_info(genome)
        for denovo in local_dict:
            local_dict[denovo]["genome"] = genome
            denovo_genes[f"{denovo}_elongated"] = local_dict[denovo]

    # Read clusters
    clusters = {}
    with open(clusters_file, "r") as f:
        for line in f:
            line_split = line.strip().split("\t")
            cluster_name = f"{line_split[0]}_elongated"
            gene_name = f"{line_split[1]}_elongated"
            if cluster_name not in clusters:
                clusters[cluster_name] = []
            clusters[cluster_name].append(gene_name)

    ####### Orfans ######
    orfans = {"true": [], "false": []}
    for cluster in clusters:
        if len(clusters[cluster]) != 1:
            continue
        gene_name = clusters[cluster][0]
        genome = denovo_genes[gene_name]["genome"]
        # Open denovo match matrix
        trg_file = os.path.join(DENSE_DIR, genome, "TRG_match_matrix.tsv")
        trg_pd = pd.read_csv(trg_file, sep="\t")
        # Get the row corresponding to the orfan de novo gene
        match_row = trg_pd.loc[trg_pd["CDS"] == gene_name].iloc[:, 2:]
        # Get the genomes for which we have a homolog cds
        homologs = list(match_row.loc[:, (match_row == "CDS").any()].columns)
        denovo_genes[gene_name]["homologs"] = homologs
        if len(homologs) == 0:
            orfans["true"].append(gene_name)
        else:
            orfans["false"].append(gene_name)

    print(f"\nOut of the {len(clusters)} clusters, there are {len(orfans["true"]) + len(orfans["false"])} for which the LCA is at the leaf level.")
    print(f"Out of these, {len(orfans['true'])} are orfans, and {len(orfans['false'])} have homologs in other genomes.\n")

    # Get info on each homolog
    false_orfans_info = {}
    for gene in orfans["false"]:
        false_orfans_info[gene] = {}
        for hom_sp in denovo_genes[gene]["homologs"]:
            # Get the homolog gene name
            blastp_file = os.path.join(DENSE_DIR, denovo_genes[gene]["genome"], "blast_out",
                                       f"TRG_multielongated_blastp_{hom_sp}_CDS_elongated.out")
            with open(blastp_file, "r") as f:
                for line in f:
                    if gene in line and not line.startswith("#"):
                        line_split = line.strip().split("\t")
                        homolog_gene_name = line_split[1]
                        eval = float(line_split[10])
                        identity = float(line_split[2])
                        qcov = float(line_split[13])
                        break

            # Look if the homolog gene is a denovo gene
            is_denovo = homolog_gene_name in denovo_genes
            # Look if the homolog gene is a TRG
            trg_file = os.path.join(DENSE_DIR, hom_sp, "TRG_match_matrix.tsv")
            trg_pd = pd.read_csv(trg_file, sep="\t")
            is_trg = homolog_gene_name in trg_pd["CDS"].values
            trg_line = ""
            if is_trg:
                trg_line = " ".join(list(trg_pd.loc[trg_pd["CDS"] == homolog_gene_name].values)[0][1:])
            
            # Add to dict
            false_orfans_info[gene][hom_sp] = {
                "homolog_gene_name": homolog_gene_name,
                "eval": eval,
                "identity": identity,
                "qcov": qcov,
                "is_denovo": is_denovo,
                "is_trg": is_trg,
                "trg_line": trg_line
            }


    # Print the false orfans info to file
    with open(os.path.join(out_folder, "false_orfans.tsv"), "w") as f:
        f.write("#gene\thomolog gene name\tis denovo\tis trg\tevalue\tidentity\tqcov\n")
        for gene in orfans["false"]:
            f.write(f">{gene}|{denovo_genes[gene]["genome"]}\n")
            for homolog in denovo_genes[gene]["homologs"]:
                name = false_orfans_info[gene][homolog]["homolog_gene_name"]
                is_denovo = false_orfans_info[gene][homolog]["is_denovo"]
                is_trg = false_orfans_info[gene][homolog]["is_trg"]
                trg_line = false_orfans_info[gene][homolog]["trg_line"]
                eval = false_orfans_info[gene][homolog]["eval"]
                identity = false_orfans_info[gene][homolog]["identity"]
                qcov = false_orfans_info[gene][homolog]["qcov"]
                f.write(f"{homolog}\t{name}\t{is_denovo}\t{is_trg}\t{eval}\t{identity}\t{qcov}\n")
                if is_trg:
                    f.write(f"{trg_line}\n")
            f.write("\n")




