"""
This script takes all the clusters that have more than one de novo in them.
For these clusters, it takes into account only the clusters that have at least one genome 
under their lca without a de novo gene.

For each of these genomes without de novo, it does a blast search against neighbours
to see if there are traces of the gene family in its genome. If it finds matches, it
does an integrity search and only takes into account the matches that have a
conserved integrity, ie. that could be coding genes.
"""

import tempfile
import os
import subprocess
from Bio import SeqIO
import pandas as pd
from my_functions.genomic_functions import get_sequence_from_loci


from my_functions.paths import FA_DIR, CDS_DIR
INPUT_FILE = "/home/eliott.tempez/Documents/M2_Stage_I2BC/results/18_clustering/genomes_in_clusters.tsv"
DE_NOVO_FASTA_FILE = "/home/eliott.tempez/Documents/M2_Stage_I2BC/results/18_clustering/good_candidates.fasta"
INTEGRITY_THRESHOLD = 70




def blast(query_fasta_file, subject_fasta_file, blast_type):
    # Create temporary output file
    with tempfile.NamedTemporaryFile(mode='w+', delete=False, suffix='.tsv') as temp_output:
        output_file = temp_output.name
    # Run the blast command
    out_command = "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen qcovhsp sframe"
    command = [
        blast_type,
        "-query", query_fasta_file,
        "-subject", subject_fasta_file,
        "-outfmt", out_command,
        "-evalue", "1e-3"
    ]
    subprocess.run(command, check=True)

    # Read the output
    # Check the file isn't blank
    if os.path.getsize(output_file) == 0:
        return None
    blast_results = pd.read_csv(output_file, sep="\t", header=None)
    blast_results.columns = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", 
                             "qstart", "qend", "sstart", "send", "evalue", "bitscore", 
                             "qlen", "qcovhsp", "sframe"]
    # Remove the temporary output file
    os.remove(output_file)
    print(blast_results)
    return blast_results



def only_keep_integral_matches(genome, blast_result, check_stops=True):
    if blast_result is None:
        print("NONE")
        return None
    
    # Check the qcov for each match
    blast_result = blast_result[blast_result["qcovhsp"] >= INTEGRITY_THRESHOLD]
    if blast_result.empty:
        print("EMPTY")
        return None
    
    print(check_stops)
    # If non-coding, check for stops
    if check_stops:
        # Extract nc sequence for each match
        for index, row in blast_result.iterrows():
            strand = "+" if int(row["sframe"]) > 0 else "-"
            seq = get_sequence_from_loci(
                genome, row["qstart"] -1, row["qend"], strand
            )
            print(seq)
            # Translate the sequence
            translated_seq = seq.translate()
            print(translated_seq)
            # Extract the length of the longest orf
            orfs = translated_seq.split("*")
            longest_orf_length = max(len(orf) for orf in orfs)
            print(longest_orf_length)
            # Calculate the longest orf qcov
            longest_orf_qcov = (longest_orf_length / row["qlen"]) * 100
            print(longest_orf_qcov)
            # If the longest orf qcov is less than the integrity threshold, remove the match
            if longest_orf_qcov < INTEGRITY_THRESHOLD:
                blast_result = blast_result.drop(index)
        if blast_result.empty:
            return None
    
    print(blast_result)
    return blast_result






if __name__ == "__main__":
    # Read the list of clusters and the genomes contained in it
    genomes_in_clusters = {}
    with open(INPUT_FILE, "r") as f:
        for line in f:
            line_split = line.strip().split("\t")
            cluster = line_split[0]
            genome = line_split[1]
            de_novo = line_split[2]
            if cluster not in genomes_in_clusters:
                genomes_in_clusters[cluster] = {"with_de_novo": [], "no_de_novo": []}
            # Add whether the genome has a de_novo or not
            if de_novo == "FALSE":
                genomes_in_clusters[cluster]["no_de_novo"].append(genome)
            else:
                genomes_in_clusters[cluster]["with_de_novo"].append((genome, de_novo))


    # Read the de novo sequences
    de_novo_sequences = {}
    for record in SeqIO.parse(DE_NOVO_FASTA_FILE, "fasta"):
        de_novo_sequences[record.name] = str(record.seq)


    # Keep only the clusters with a coverage < 100%
    for cluster in genomes_in_clusters:
        if genomes_in_clusters[cluster]["no_de_novo"] == []:
            continue
        print(cluster)

        ## Blast de novo sequences against genomes with no de novo sequences
        # Create a temporary fasta file for all queries
        queries_names = [de_novo for _, de_novo in genomes_in_clusters[cluster]["with_de_novo"]]
        with tempfile.NamedTemporaryFile(mode='w+', delete=False, suffix='.fa') as queries_fasta:
            for query_name in queries_names:
                queries_fasta.write(f">{query_name}\n{de_novo_sequences[query_name]}\n")
        queries_fasta_file = queries_fasta.name

        # For each genome in the cluster with no de novo, do the blast
        for no_denovo_genome in genomes_in_clusters[cluster]["no_de_novo"]:
            print(no_denovo_genome)
            no_denovo_fasta_file_genome = f"{FA_DIR}{no_denovo_genome}.fa"
            no_denovo_fasta_file_cds = f"{CDS_DIR}{no_denovo_genome}_CDS.faa"

            # Start with blastp
            print("blastp")
            result = blast(queries_fasta_file, no_denovo_fasta_file_cds, "blastp")
            # Do integrity analysis
            result = only_keep_integral_matches(no_denovo_genome, result, check_stops=False)

            # Then blastn if no result with blastp
            if result:
                continue
            print("\ntblastn")
            result = blast(queries_fasta_file, no_denovo_fasta_file_genome, "tblastn")
            # Do integrity analysis
            result = only_keep_integral_matches(no_denovo_genome, result, check_stops=True)



            print()

        # Remove the temporary fasta file
        os.remove(queries_fasta_file)
