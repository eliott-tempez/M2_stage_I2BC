import os
import glob
import re
import random
import subprocess
import tempfile
import argparse
import pandas as pd
import Bio.SeqIO as SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import six_frame_translations
import concurrent.futures


OUTPUT_DIR = "/home/eliott.tempez/Documents/M2_Stage_I2BC/results/calculate_conservation/"
GENOMES_LIST = "/home/eliott.tempez/Documents/M2_Stage_I2BC/scripts/genera_archaea/genomes_list.txt"
GENERA_DIR = "/home/eliott.tempez/Documents/archaea_data/genera/out/"
DENSE_DIR = "/home/eliott.tempez/Documents/archaea_data/dense/"
DATA_DIR = "/home/eliott.tempez/Documents/archaea_data/complete_122/"
TRG_RANK = 7.0
NCPUS = 12



############################## FUNCTIONS ##############################
def initialise_conservation_df(genomes):
    """Initialise the dataframe with empty rows and columns"""
    col = ["n_denovo", "n_trg", "n_cds", "n_intergenic", "ssearch_f0", "ssearch_f1", "ssearch_f2", "ssearch_f0_comp", "ssearch_f1_comp", "ssearch_f2_comp"]
    conservation_df = pd.DataFrame(columns=col, index=genomes)
    return conservation_df


def get_seqs_from_gene_names(focal_sp, gene_names, extension_type):
    """Get the sequences of the genes from the gene names"""
    fasta_file = os.path.join(DATA_DIR, "CDS/" + focal_sp + "_CDS." + extension_type)
    fasta = SeqIO.parse(fasta_file, "fasta")
    seqs = {}
    for record in fasta:
        if record.name in gene_names:
            seqs[record.name] = record.seq
    return seqs


def extract_de_novo(focal_species):
    """Extract the name of all de novo genes of the focal species"""
    # Make sure the file exists and read it
    de_novo_file = os.path.join(DENSE_DIR, focal_species, "denovogenes.tsv")
    if os.path.exists(de_novo_file):
        # Get the number of lines in file
        with open(de_novo_file, "r") as f:
            de_novo_lines = f.readlines()
            if len(de_novo_lines) == 1:
                return []
        de_novo_df = pd.read_csv(de_novo_file, sep="\t", header=0)
    else:
        raise FileNotFoundError(f"No file matching pattern {de_novo_file}")
    # Extract the gene names
    de_novo_genes = de_novo_df["CDS"].tolist()
    return de_novo_genes


def extract_focal_TRGs(focal_species, trg_threshold):
    """Extract at most 1000 TRGs of the focal species"""
    # Make sure the file exists and read it
    trg_pattern = os.path.join(GENERA_DIR, focal_species, "*gene_ages.tsv")
    trg_files = glob.glob(trg_pattern)
    if trg_files:
        trg_file = trg_files[0]
        trg_df = pd.read_csv(trg_file, sep="\t", header=0)
    else:
        raise FileNotFoundError(f"No file matching pattern {trg_pattern}")
    # Extract the TRGs
    focal_TRGs = trg_df[trg_df["rank"] >= trg_threshold]["#gene"].tolist()
    # Keep 1000 at most
    if len(focal_TRGs) > 1000:
        focal_TRGs = random.sample(focal_TRGs, 1000)
    return focal_TRGs


def extract_focal_CDSs(focal_species):
    """Extract at most 1000 CDSs of the focal species"""
    # Extract all CDSs >=302 nt of the focal species
    CDS_fasta_file = os.path.join(DATA_DIR, "CDS/" + focal_species + "_CDS.fna")
    CDS_fasta = SeqIO.parse(CDS_fasta_file, "fasta")
    focal_CDSs = {}
    for record in CDS_fasta:
        if len(record.seq) >= 302:
            focal_CDSs[record.name] = record.seq
    # Keep 1000 at most
    focal_CDSs_names = list(focal_CDSs.keys())
    focal_CDSs = {k: focal_CDSs[k] for k in focal_CDSs_names}
    print("Here are the CDSs for genome a:")
    for k, v in focal_CDSs.items():
        print(f"{k}: {v} - length: {len(v)}")
    print("\n")
    return focal_CDSs


def extract_intergenic(species):
        """Extract all intergenic sequences"""
        # Make sure the file exists and read it
        gbk_pattern = os.path.join(DATA_DIR, "annotation", species + ".gb*")
        gbk_files = [f for f in glob.glob(gbk_pattern) if not f.endswith(".fai")]
        if gbk_files:
            gbk_file = gbk_files[0]
            gbk_content = list(SeqIO.parse(gbk_file, "genbank"))
        else:
            raise FileNotFoundError(f"No file matching pattern {gbk_pattern}")

        # Extract data for each contig
        intergenic_segments = {}
        for record in gbk_content:
            contig_length = len(record.seq)
            intergenic_segments[record.name] = []
            # Get all nucl positions in CDSs from the genbank file
            coding_loci = []
            for feature in record.features:
                if feature.type == "CDS":
                    # Make sure we don't include the bugged CDS that is as long as the genome
                    if feature.location.end - feature.location.start < contig_length:
                        coding_loci += list(range(feature.location.start, feature.location.end))
            # Get all the integers that are not in the list of coding loci
            full_set = set(range(min(coding_loci), max(coding_loci) + 1))
            noncoding_positions = list(full_set - set(coding_loci))
            # Convert to consecutive ranges
            noncoding_positions.sort()
            start = noncoding_positions[0]
            for i in range(1, len(noncoding_positions)):
                if noncoding_positions[i] != noncoding_positions[i - 1] + 1:
                    intergenic_segments[record.name].append((start, noncoding_positions[i - 1] + 1))
                    start = noncoding_positions[i]
            intergenic_segments[record.name].append((start, noncoding_positions[-1] + 1))

        # Extract the corresponding sequences from fasta file
        fa_file = os.path.join(DATA_DIR, "fasta_renamed/" + species + ".fa")
        fa_content = list(SeqIO.parse(fa_file, "fasta"))
        intergenic_dict = {}
        i = 0
        for contig in fa_content:
            if contig.name in intergenic_segments:
                for j in range(len(intergenic_segments[contig.name])):
                    i += 1
                    start, end = intergenic_segments[contig.name][j]
                    intergenic_dict[f"{contig.name}_{i}"] = contig.seq[start:end]
        return intergenic_dict


def cut_chunks(seqs, length):
    """Cut the CDSs in lengths-long chunks randomly"""
    for k, v in seqs.items():
        len_over = len(v) - length
        if len_over == 0:
            len_over = 1
        # Chose chunk randomly in-frame
        possible_starts = list(range(0, len_over, 3))
        start = random.choice(possible_starts)
        seqs[k] = v[start:(start + length)]
    return seqs


def get_db(species, colname, db):
    """Get the database fasta file for the species depending on the type of sequences"""
    if colname in ["n_denovo", "n_trg", "n_cds"]:
        db_fasta_file = os.path.join(DATA_DIR, "CDS/" + species + "_CDS.faa")

    elif colname in ["n_f1", "n_f2", "ssearch_f0", "ssearch_f1", "ssearch_f2", "ssearch_f0_comp", "ssearch_f1_comp", "ssearch_f2_comp"]:
        # Get the dna sequences for the matching genes
        gene_names = list(db.values())
        gene_dict_nt = get_seqs_from_gene_names(species, gene_names, "fna")
        print("Here are the sequences for the homologs:")
        for k, v in gene_dict_nt.items():
            print(f"{k}: {v} - length: {len(v)}")
        print("\n")
        if "comp" in colname:
            gene_dict_nt = {k: v.reverse_complement() for k, v in gene_dict_nt.items()}
        # Translate in all 6 frames
        if "0" in colname:
            gene_dict = {k: v[:-3].translate(table=11) for k, v in gene_dict_nt.items()}
        elif "1" in colname:
            gene_dict = {k: v[1:-2].translate(table=11) for k, v in gene_dict_nt.items()}
        elif "2" in colname:
            gene_dict = {k: v[2:-1].translate(table=11) for k, v in gene_dict_nt.items()}
        # Write the sequences to a temporary file
        with tempfile.NamedTemporaryFile(mode='w', delete=False) as db_file:
            for gene, seq in gene_dict.items():
                db_file.write(f">{gene}\n{seq}\n")
            db_fasta_file = db_file.name

    elif colname == "n_intergenic":
        # Extract intergenic sequences
        intergenic_dict = extract_intergenic(species)
        # Write the sequences to a temporary file
        with tempfile.NamedTemporaryFile(mode='w', delete=False) as db_file:
            for gene, seq in intergenic_dict.items():
                db_file.write(f">{gene}\n{seq}\n")
            db_fasta_file = db_file.name

    return db_fasta_file


def run_blast(query_sequences, db_faa_file, blast_type, keep_homologs, db):
    """Run BLAST for the query sequences against the database and return the number of matches"""
    # Create a temporary file for the query sequences
    with tempfile.NamedTemporaryFile(mode='w', delete=False) as query_file:
        for trg_id, trg_seq in query_sequences.items():
            query_file.write(f">{trg_id}\n{trg_seq}\n")
        query_file_path = query_file.name
    # Create a temporary file for the BLAST output
    with tempfile.NamedTemporaryFile(mode='w', delete=False) as output_file:
        output_file_path = output_file.name

    print("Here is the content of the fasta file for species b against which those CDSs are run:")
    with open(db_faa_file, "r") as f:
        content = f.read()
        print(content)
    print("\n")

    # Run the BLAST
    try:
        output = subprocess.run([blast_type, "-query", query_file_path, "-subject", db_faa_file, "-out", output_file_path, "-outfmt", "6 qseqid sseqid qlen evalue qcovs", "-evalue", "1e-3", "-num_threads", str(NCPUS)], capture_output=True, check=True)
    except subprocess.CalledProcessError as e:
        print(f"BLAST command failed with return code {e.returncode}: {e.stderr.decode()}")
        raise

    # Check the file isn't blank
    with open(output_file_path, "r") as f:
        content = f.read().strip()
    if content == "":
        os.remove(query_file_path)
        os.remove(output_file_path)
        return 0, None

    # Read the output file
    result = pd.read_csv(output_file_path, sep="\t", header=None)
    result.columns = ["qseqid", "sseqid", "qlen", "evalue", "qcov"]
    # If we have a match couple already known
    if db:
        # Keep only the matches in the right CDS
        result = result[result.apply(lambda row: db.get(row ["qseqid"]) == row ["sseqid"], axis=1)]
    # Keep only the best hit for each query
    result = result.sort_values("evalue").drop_duplicates("qseqid")
    print("Here is the result of the BLAST:")
    print(result)
    # Keep homolog sequences
    if keep_homologs:
        homologs = {row["qseqid"]: row["sseqid"] for i, row in result.iterrows()}
        return len(result), homologs
    # Remove the temporary files
    os.remove(query_file_path)
    os.remove(output_file_path)
    # Return number of matches
    return len(result), None


def process_conservation_for_species(species, focal_sp, conservation_df, colname, query_sequences, db=None):
    """Process the conservation for the focal species and add to the dataframe"""
    # Keep homolog sequences for the CDSs
    homologs, keep_homologs = None, False
    if colname == "n_cds":
        keep_homologs = True

    if species == focal_sp:
        conservation_df.loc[species, colname] = len(query_sequences)
        return conservation_df, homologs
    
    # For all the other species, run BLAST for all the sequences
    # Get the database
    db_sp = db[species] if db else None
    db_fasta_file = get_db(species, colname, db_sp)
    # Get the type of BLAST to run
    if colname in ["n_denovo", "n_trg", "n_cds"]:
        blast_type = "blastp"
    elif colname == "n_intergenic":
        blast_type = "tblastx"

    # run blast
    nb_matches, homologs = run_blast(query_sequences, db_fasta_file, blast_type, keep_homologs, db_sp)
    print(f"#Species {species}: {nb_matches} matches\n")
    conservation_df.loc[species, colname] = nb_matches
    
    # Delete db file if it was created
    if colname in ["n_f1", "n_f2", "n_intergenic"]:
        os.remove(db_fasta_file)

    return conservation_df, homologs 


def process_conservation_parallel(focal_sp, conservation_df, colname, query_sequences, db=None):
    """Process the conservation for all species in parallel"""
    homologs, keep_homologs = None, False
    if colname == "n_cds":
        keep_homologs = True
        homologs = {}

    with concurrent.futures.ThreadPoolExecutor(max_workers=NCPUS) as executor:
        futures = {executor.submit(process_conservation_for_species, species, focal_sp, conservation_df, colname, query_sequences, db): species for species in conservation_df.index}
        for future in concurrent.futures.as_completed(futures):
            species = futures[future]
            try:
                conservation_df, hom = future.result()
                if keep_homologs:
                    homologs[species] = hom
            except Exception as exc:
                print(f"Species {species} generated an exception: {exc}")

    return conservation_df, homologs


def ssearch_for_species(species, focal_sp, conservation_df, colname, query_sequences, homologs):
    """Run ssearch for each cds for the focal species"""
    if species == focal_sp:
        print(f"#Species {species}: {len(query_sequences)} matches")
        conservation_df.loc[species, colname] = len(query_sequences)
        return conservation_df

    nb_matches = 0
    db = get_db(species, colname, homologs)
    for query_gene in homologs:
        db_gene = homologs[query_gene]
        print(f"Query gene: {query_gene} - Homolog: {db_gene}")
        print(f"Query sequence: {query_sequences[query_gene]}")

        # Write temporary files for the query and the homolog
        with tempfile.NamedTemporaryFile(mode='w', delete=False) as query_file:
            query_file.write(f">{query_gene}\n{str(query_sequences[query_gene])}\n")
            query_file_path = query_file.name
        fasta = SeqIO.parse(db, "fasta")
        for record in fasta:
            if record.name == db_gene:
                db_seq = str(record.seq)
                pass
        with tempfile.NamedTemporaryFile(mode='w', delete=False) as db_file:
            db_file.write(f">{db_gene}\n{str(db_seq)}\n")
            db_file_path = db_file.name
        # Write temporary file for the output
        with tempfile.NamedTemporaryFile(mode='w', delete=False) as output_file:
            output_file_path = output_file.name
        print(f"Homolog sequence: {str(db_seq)}")

        # Run ssearch
        try:
            output = subprocess.run(["ssearch36", "-m", "8", "-E", "1e-3", query_file_path, db_file_path, "-O", output_file_path], capture_output=True, check=True)
        except subprocess.CalledProcessError as e:
            print(f"Water command failed with return code {e.returncode}: {e.stderr.decode()}")
            raise
        
        print("Here is the content of the output file:")
        with open(output_file_path, "r") as f:
            content = f.read()
            print(content)
        print("\n")

        # Check the file isn't blank
        with open(output_file_path, "r") as f:
            content = f.read().strip()
        if content == "":
            os.remove(query_file_path)
            os.remove(db_file_path)
            os.remove(output_file_path)
            continue


        # If it isn't, then we have a match
        nb_matches += 1
        # Remove the temporary files
        os.remove(query_file_path)
        os.remove(db_file_path)
        os.remove(output_file_path)

    print(f"#Species {species}: {nb_matches} matches")
    conservation_df.loc[species, colname] = nb_matches
    return conservation_df


def run_ssearch_parallel(focal_sp, conservation_df, colname, query_sequences, homologs):
    """Run ssearch in parallel"""
    with concurrent.futures.ThreadPoolExecutor(max_workers=NCPUS) as executor:
        futures = {executor.submit(ssearch_for_species, species, focal_sp, conservation_df, colname, query_sequences, homologs[species]): species for species in conservation_df.index}
        for future in concurrent.futures.as_completed(futures):
            species = futures[future]
            try:
                conservation_df = future.result()
            except Exception as exc:
                print(f"Species {species} generated an exception: {exc}")
    return conservation_df





############################## MAIN ##############################
if __name__ == "__main__":
    # Parse command line arguments for focal species
    parser = argparse.ArgumentParser(description='Calculate TRG conservation.')
    parser.add_argument('--focal_species', type=str, required=True, help='Focal species identifier')
    args = parser.parse_args()
    FOCAL_SPECIES = args.focal_species

    # Read the list of genomes
    """with open(GENOMES_LIST, "r") as f:
        genomes = f.readline().split()
    genomes = [re.sub('"', '', g) for g in genomes]"""
    genomes = ["a", "b"]
    # Initiate the dataframe
    conservation_df = initialise_conservation_df(genomes)
    print("\n")

    print("This is a test run to check the code works properly.")
    print("The code is run on 2 species only: a and b. The focal species will be a. The lines starting with # are the original output lines fot the program.\n")


    ########################################
    ############### De novo ################
    ########################################
    """print("Calculating de novo conservation...")
    # Extract the name of all de novo genes of the focal species
    de_novo_genes = extract_de_novo(FOCAL_SPECIES)
    if de_novo_genes == []:
        conservation_df["n_denovo"] = 0
        print("No de novo genes found for the focal species")
    else:
        # Get the protein sequences for each de novo gene
        de_novo_seqs = get_seqs_from_gene_names(FOCAL_SPECIES, de_novo_genes, "faa")
        print(f"Extracted {len(de_novo_seqs)} de novo genes for the focal species {FOCAL_SPECIES}\n")
        # Calculate the de novo conservation and add to dataframe
        conservation_df, _ = process_conservation_parallel(FOCAL_SPECIES, conservation_df, "n_denovo", de_novo_seqs)
    print("\n\n")"""



    ########################################
    ################# CDSs #################
    ########################################
    print("#Calculating CDS conservation...")
    # Extract 1000 CDss of the focal species
    focal_CDSs_nt = extract_focal_CDSs(FOCAL_SPECIES)
    print(f"#Extracted {len(focal_CDSs_nt)} CDSs for the focal species {FOCAL_SPECIES}\n")
    # Translate the CDSs
    focal_CDSs = {k: v.translate(table=11) for k, v in focal_CDSs_nt.items()}
    print(f"Here are the translated CDSs:")
    for k, v in focal_CDSs.items():
        print(f"{k}: {v} - length: {len(v)}")
    print("\n")
    # Calculate the CDS conservation and add to dataframe + keep homolog sequences
    conservation_df, homologs = process_conservation_parallel(FOCAL_SPECIES, conservation_df, "n_cds", focal_CDSs)
    print("Here is the dictionnary of homologs:")
    print(homologs)
    print("\n\n")



    ########################################
    ################# TRGs #################
    ########################################
    """print("Calculating TRG conservation...")
    # Extract the name of 1000 (at most) TRGs of the focal species
    focal_TRG_genes = extract_focal_TRGs(FOCAL_SPECIES, TRG_RANK)
    # Get the protein sequences for each TRG
    focal_TRGs = get_seqs_from_gene_names(FOCAL_SPECIES, focal_TRG_genes, "faa")
    print(f"Extracted {len(focal_TRGs)} TRGs for the focal species {FOCAL_SPECIES}\n")
    # Calculate the TRG conservation and add to dataframe
    conservation_df, _ = process_conservation_parallel(FOCAL_SPECIES, conservation_df, "n_trg", focal_TRGs)
    print("\n\n")"""



    ########################################
    ############## Intergenic ##############
    ########################################
    """print("Calculating intergenic conservation...")
    # Extract 1000 intergenic 100 sequences of the focal species
    focal_intergenic_nt = extract_intergenic(FOCAL_SPECIES)
    # Keep only the intergenic sequences of at least 100 nt and cut chunks
    focal_intergenic = {k: v for k, v in focal_intergenic_nt.items() if len(v) >= 100}
    focal_intergenic = cut_chunks(focal_intergenic, 100)
    # Keep 1000 at most
    focal_intergenic_names = list(focal_intergenic.keys())
    if len(focal_intergenic) > 1000:
        focal_intergenic_names = random.sample(focal_intergenic_names, 1000)
    focal_intergenic = {k: focal_intergenic[k] for k in focal_intergenic_names}
    print(f"Extracted {len(focal_intergenic)} intergenic sequences for the focal species {FOCAL_SPECIES}\n")
    # Calculate the intergenic conservation and add to dataframe
    conservation_df, _ = process_conservation_parallel(FOCAL_SPECIES, conservation_df, "n_intergenic", focal_intergenic)
    print("\n\n")"""



    ########################################
    ################ Ssearch ###############
    ########################################
    # Cut 300-long chunks randomly for the CDSs
    print("Now that we have the homolog couples, we can evaluate the conservation for all 6 frames of the CDSs.")
    focal_CDSs_nt = cut_chunks(focal_CDSs_nt, 300)
    print("Here is a random 300-long chunk for each CDS in a, in the +0 frame:")
    for k, v in focal_CDSs_nt.items():
        print(f"{k}: {v} - length: {len(v)}")
    print("\n")
    focal_CDSs_nt_comp = {k: v.reverse_complement() for k, v in focal_CDSs_nt.items()}
    print("Here is the reverse complements for the same chunks:")
    for k, v in focal_CDSs_nt_comp.items():
        print(f"{k}: {v} - length: {len(v)}")
    print("\n")
    # Translate in the 6 frames
    print("Here are all 6 frame translations for the CDSs:")
    focal_CDSs_f0 = {k: v[:-3].translate(table=11) for k, v in focal_CDSs_nt.items()}
    focal_CDSs_f1 = {k: v[1:-2].translate(table=11) for k, v in focal_CDSs_nt.items()}
    focal_CDSs_f2 = {k: v[2:-1].translate(table=11) for k, v in focal_CDSs_nt.items()}
    focal_CDSs_f0_comp = {k: v[:-3].translate(table=11) for k, v in focal_CDSs_nt_comp.items()}
    focal_CDSs_f1_comp = {k: v[1:-2].translate(table=11) for k, v in focal_CDSs_nt_comp.items()}
    focal_CDSs_f2_comp = {k: v[2:-1].translate(table=11) for k, v in focal_CDSs_nt_comp.items()}
    print("Frame 0:")
    for k, v in focal_CDSs_f0.items():
        print(f"{k}: {v} - length: {len(v)}")
    print("Frame 1:")
    for k, v in focal_CDSs_f1.items():
        print(f"{k}: {v} - length: {len(v)}")
    print("Frame 2:")
    for k, v in focal_CDSs_f2.items():
        print(f"{k}: {v} - length: {len(v)}")
    print("Frame -0:")
    for k, v in focal_CDSs_f0_comp.items():
        print(f"{k}: {v} - length: {len(v)}")
    print("Frame -1:")
    for k, v in focal_CDSs_f1_comp.items():
        print(f"{k}: {v} - length: {len(v)}")
    print("Frame -2:")
    for k, v in focal_CDSs_f2_comp.items():
        print(f"{k}: {v} - length: {len(v)}")
    print("\n")

    # Run ssearch for the 6 frames
    print(f"#Running ssearch for all 6 frames of the {len(focal_CDSs_nt)} CDSs...\n")
    print("#Frame 0...")
    conservation_df = run_ssearch_parallel(FOCAL_SPECIES, conservation_df, "ssearch_f0", focal_CDSs_f0, homologs)
    print("\n#Frame 1...")
    conservation_df = run_ssearch_parallel(FOCAL_SPECIES, conservation_df, "ssearch_f1", focal_CDSs_f1, homologs)
    print("\n#Frame 2...")
    conservation_df = run_ssearch_parallel(FOCAL_SPECIES, conservation_df, "ssearch_f2", focal_CDSs_f2, homologs)
    print("\n#Frame -0...")
    conservation_df = run_ssearch_parallel(FOCAL_SPECIES, conservation_df, "ssearch_f0_comp", focal_CDSs_f0_comp, homologs)
    print("\n#Frame -1...")
    conservation_df = run_ssearch_parallel(FOCAL_SPECIES, conservation_df, "ssearch_f1_comp", focal_CDSs_f1_comp, homologs)
    print("\n#Frame -2...")
    conservation_df = run_ssearch_parallel(FOCAL_SPECIES, conservation_df, "ssearch_f2_comp", focal_CDSs_f2_comp, homologs)
    print("\n\n")





    # Old code to calculate +1/+2 conservation with blastp
    """########################################
    ############ CDSs / +1 / +2 ############
    ########################################
    print("Calculating CDS conservation...")
    # Extract 1000 CDss of the focal species
    focal_CDSs_nt = extract_focal_CDSs(FOCAL_SPECIES)
    print(f"Extracted {len(focal_CDSs_nt)} CDSs for the focal species {FOCAL_SPECIES}\n")
    # Cut 302-long chunks randomly
    focal_CDSs_nt = cut_chunks(focal_CDSs_nt, 302)
    # Translate to get 100 aa-long protein sequences
    focal_CDSs = {k: v[:300].translate(table=11) for k, v in focal_CDSs_nt.items()}
    # Calculate the CDS conservation and add to dataframe + keep homolog sequences
    conservation_df, homologs = process_conservation_parallel(FOCAL_SPECIES, conservation_df, "n_cds", focal_CDSs)
    print("\n\n")


    # Frames 1 and 2
    # Translate alt frames to get 100 aa-long protein sequences
    focal_CDS_frame_1 = {k: v[1:301].translate(table=11) for k, v in focal_CDSs_nt.items()}
    focal_CDS_frame_2 = {k: v[2:302].translate(table=11) for k, v in focal_CDSs_nt.items()}
    print("Calculating conservation for +1 frame...")
    # Calculate the conservation for +1 frame and add to dataframe
    conservation_df, _ = process_conservation_parallel(FOCAL_SPECIES, conservation_df, "n_f1", focal_CDS_frame_1, homologs)
    print("\n\n")
    print("Calculating conservation for +2 frame...")
    # Calculate the conservation for +2 frame and add to dataframe
    conservation_df, _ = process_conservation_parallel(FOCAL_SPECIES, conservation_df, "n_f2", focal_CDS_frame_2, homologs)
    print("\n\n")"""



    # Write result to file
    """conservation_df.to_csv(os.path.join(OUTPUT_DIR, f"conservation_df_{FOCAL_SPECIES}.tsv"), sep="\t")"""

