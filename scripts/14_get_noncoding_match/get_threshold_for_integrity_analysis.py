import re
import random
import os
import tempfile
import subprocess
from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd
import matplotlib.pyplot as plt

from my_functions.paths import GENOMES_LIST, CDS_DIR, FA_DIR
IORFS_FILE = "/home/eliott.tempez/Documents/M2_Stage_I2BC/results/17_compare_denovo_sequences/iorfs.txt"
NEIGHBOUR_GENOMES_FILE = "/home/eliott.tempez/Documents/M2_Stage_I2BC/results/2_initial_data_analysis_figures/neighbour_genomes.csv"
OUT_FOLDER = "/home/eliott.tempez/Documents/M2_Stage_I2BC/results/14_get_noncoding_match/"


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
        "-out", output_file,
        "-evalue", "1e-3"
    ]
    subprocess.run(command, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

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
    return blast_results


def get_qcovs_from_best_match(blast_result):
    if blast_result is None:
        return []
    
    # Sort by sseqid then evalue then qcovhsp
    blast_result = blast_result.sort_values(by=["qseqid", "evalue", "qcovhsp"], ascending=[True, True, False])
    # Keep only the best match for each sseqid
    blast_result = blast_result.drop_duplicates(subset=["qseqid"], keep="first")
    qcovs = blast_result["qcovhsp"].tolist()
    return qcovs




if __name__ == "__main__":
    # Read list of genomes
    with open(GENOMES_LIST, "r") as f:
        genomes = f.readline().split()
    genomes = [re.sub('"', '', g) for g in genomes]

    # Read closest genome for each genome
    neighbours = {}
    with open(NEIGHBOUR_GENOMES_FILE, "r") as f:
        for line in f:
            line_split = line.strip().split("\t")
            genome = line_split[0]
            neighbour = line_split[1]
            neighbours[genome] = neighbour

    # Set results
    qcovs_iorfs = []
    qcovs_cdss = []

    # For each genome, extract 100 random CDSs and iORFs
    for genome in genomes:
        # Extract 100 random iORFs
        with open(IORFS_FILE, "r") as f:
            for line in f:
                if line.startswith(f">{genome}"):
                    iorfs_str = f.readline().strip()
                    break
        iorfs = iorfs_str.split()
        iorfs_sample = random.sample(iorfs, 100)
        # Translate sequences
        iorfs_sample = [str(Seq(iorf).translate()) for iorf in iorfs_sample]
        # Create temp fasta file for iORFs
        iorfs_temp_file = tempfile.NamedTemporaryFile(delete=False, suffix=".fasta")
        i = 0
        with open(iorfs_temp_file.name, "w") as temp_f:
            for iorf in iorfs_sample:
                i += 1
                temp_f.write(f">iorf_{i}\n{iorf}\n")
                
        # Extract 100 random CDSs
        fasta_file = f"{CDS_DIR}{genome}_CDS.faa"
        all_cdss = list(SeqIO.parse(fasta_file, "fasta"))
        cdss_sample = [str(cds.seq) for cds in random.sample(all_cdss, 100)]
        # Create temp fasta file for CDSs
        cdss_temp_file = tempfile.NamedTemporaryFile(delete=False, suffix=".fasta")
        i = 0
        with open(cdss_temp_file.name, "w") as temp_f:
            for cds in cdss_sample:
                i += 1
                temp_f.write(f">cds_{i}\n{cds}\n")

        ## Blast iorfs and cdss against the closest genome
        neighbour_cds_fasta = f"{CDS_DIR}{neighbours[genome]}_CDS.faa"
        neighbour_genome_fasta = f"{FA_DIR}{neighbours[genome]}.fa"

        # Blast iORFs
        result_iorf = blast(iorfs_temp_file.name, neighbour_genome_fasta, "tblastn")
        # Keep only best match
        qcovs_iorfs += get_qcovs_from_best_match(result_iorf)

        # Blast CDSs
        result_cds = blast(cdss_temp_file.name, neighbour_cds_fasta, "blastp")
        # Keep only best match
        qcovs_cdss += get_qcovs_from_best_match(result_cds)

        # Remove temporary files
        os.remove(iorfs_temp_file.name)
        os.remove(cdss_temp_file.name)


    # Plot the results
    # Get bins in common
    bins = [i for i in range(0, 101, 2)]
    # Create histogram
    plt.hist(qcovs_iorfs, bins=bins, color='purple', alpha=0.7, label='iORFs')
    plt.hist(qcovs_cdss, bins=bins, color='orange', alpha=0.7, label='CDSs')
    # Add labels and legend
    plt.xlabel('Query Coverage (%)')
    plt.ylabel('Frequency')
    plt.title('Distribution of Query Coverage for iORFs and CDSs')
    plt.legend()
    # Save to file
    plt.savefig(os.path.join(OUT_FOLDER, "qcovs_distribution.png"))


