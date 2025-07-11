import os
import re
import tempfile
import subprocess
import numpy as np
import pandas as pd
import psa
from my_functions.genomic_functions import extract_denovo_info, get_sequence_from_loci


from my_functions.paths import GENOMES_LIST
GOOD_CANDIDATES_FILE = "/home/eliott.tempez/Documents/M2_Stage_I2BC/results/14_get_noncoding_match/good_candidates.txt"
OUTGROUP_NUMBER = 2


def get_extended_matched_seq(genome, contig, start, end, strand, missing_cter, missing_nter):
    """
    Get the sequence (nucleotides) of the match, extended on each relevant side.
    """
    # Extend sequence on both sides
    if strand == "+":
        add_right = (missing_cter * 3) * 2 if missing_cter > 4 else 0
        add_left = (missing_nter * 3) * 2 if missing_nter > 4 else 0
    elif strand == "-":
        add_left = (missing_cter * 3) * 2 if missing_cter > 4 else 0
        add_right = (missing_nter * 3) * 2 if missing_nter > 4 else 0
    match_start = start - add_left
    match_end = end + add_right
    # Extract extended sequence
    seq = get_sequence_from_loci(genome, contig, match_start, match_end, strand)
    # Get the limits of the actual match
    limit_start = add_left if strand == "+" else add_right
    limit_end = (len(seq) - add_right) if strand == "+" else (len(seq) - add_left)
    return seq, limit_start, limit_end


def in_frame(pos, start_or_end):
    """Get the closest number that is a 3 multiple"""
    if pos % 3 == 0:
        return pos
    is_start = start_or_end == "start"
    if is_start:
        if (pos - 1) % 3 == 0:
            return pos - 1
        return pos - 2
    if (pos + 1) % 3 == 0:
        return pos + 1
    return pos + 2


def smith_waterman(ref_seq, subject_seq):
    aln_dict = {}
    # Remove stops
    ref_seq = re.sub(r"\*", "", str(ref_seq))
    subject_seq = re.sub(r"\*", "", str(subject_seq))
    # Align
    aln = psa.water(moltype = "prot", qseq = ref_seq, sseq = subject_seq)
    # Store alignment in dict
    aln_dict["pval"] = aln.pvalue()
    aln_dict["psim"] = aln.psimilarity
    aln_dict["length"] = aln.length
    aln_dict["qstart"], aln_dict["qend"] = aln.qstart - 1, aln.qend
    aln_dict["sstart"], aln_dict["send"] = aln.sstart - 1, aln.send
    aln_dict["raw"] = aln.raw
    return aln_dict


def is_significant(alignment):
    if not alignment:
        return False
    aln_len = alignment["length"]
    aln_pval = alignment["pval"]
    if aln_len >= 5 and aln_pval <= 1e-2:
        return True
    return False


def get_absolute_match(aln, start_pos_query, start_pos_subject, is_already_nu):
    # For the qstart
    relative_qstart = aln["qstart"]
    absolute_qstart = relative_qstart + start_pos_query
    aln["qstart"] = absolute_qstart
    # For the sstart
    frame = aln["frame"]
    relative_sstart = aln["sstart"]
    if is_already_nu:
        absolute_sstart = relative_sstart + start_pos_subject
    else:
        absolute_sstart = relative_sstart * 3 + start_pos_subject + frame
    aln["sstart"] = absolute_sstart
    # For the qend
    relative_qend = aln["qend"]
    absolute_qend = relative_qend + start_pos_query
    aln["qend"] = absolute_qend
    # For the send
    relative_send = aln["send"]
    if is_already_nu:
        absolute_send = relative_send + start_pos_subject
    else:
        absolute_send = relative_send * 3 + start_pos_subject + frame
    aln["send"] = absolute_send
    return aln


def get_segments_from_set(indexes, threshold):
    """
    Take a set of continuous numbers and extract the segments of consecutive numbers.
    """
    continuous_segments = []
    unmatched_indexes = sorted(list(indexes))
    if len(indexes) < threshold:
        return []
    segment_start = unmatched_indexes[0]
    for i in range(len(unmatched_indexes) - 1):
        if unmatched_indexes[i+1] - unmatched_indexes[i] > 1:
            segment_end = unmatched_indexes[i]
            if segment_end - segment_start >= threshold - 1:
                continuous_segments.append((segment_start, segment_end))
            segment_start = unmatched_indexes[i+1]
    if unmatched_indexes[i+1] - unmatched_indexes[i] == 1:
        if unmatched_indexes[i+1] - segment_start >= threshold - 1:
            continuous_segments.append((segment_start, unmatched_indexes[i+1]))
    return continuous_segments


def get_uncovered_segments(matches, query_len, subject_len):
    """
    Get all the segments > 4 aa in the query and subject sequences that haven't found a match.
    """
    covered_segments_query = []
    covered_segments_subject = []
    for match in matches:
        qstart = match["qstart"]
        qend = match["qend"]
        covered_segments_query += list(range(qstart, qend))
        sstart = match["sstart"]
        send = match["send"]
        covered_segments_subject += list(range(sstart, send))
    uncovered_pos_query = set(range(query_len)) - set(covered_segments_query)
    uncovered_pos_subject = set(range(subject_len)) - set(covered_segments_subject)
    uncovered_segments_query = get_segments_from_set(uncovered_pos_query, 5)
    uncovered_segments_subject = get_segments_from_set(uncovered_pos_subject, 15)
    # Get the start and end points
    qstarts, qends, sstarts, sends = [], [], [], []
    for segment in uncovered_segments_query:
        qstarts.append(segment[0])
        qends.append(segment[1])
    for segment in uncovered_segments_subject:
        sstarts.append(segment[0])
        sends.append(segment[1])
    return qstarts, qends, sstarts, sends


def recursively_align(query_seq, subject_seq_nu, start_pos_query_l, end_pos_query_l, start_pos_subject_l, end_pos_subject_l, matches):
    best_aln = None
    # For all segment combinations
    best_pval = 1
    best_psim = 0
    best_length = 0
    for i in range(len(start_pos_query_l)):
        start_pos_query = start_pos_query_l[i]
        end_pos_query = end_pos_query_l[i]
        for j in range(len(start_pos_subject_l)):
            start_pos_subject = in_frame(start_pos_subject_l[j], "start")
            end_pos_subject = in_frame(end_pos_subject_l[j], "end")
            # Get the cut sequences
            query_seq_segment = query_seq[start_pos_query:end_pos_query]
            subject_seq_nu_segment = subject_seq_nu[start_pos_subject:end_pos_subject]
            # Align in all 3 frames and keep the best one
            for frame in [0, 1, 2]:
                subject_seq_aa = subject_seq_nu_segment[frame:(frame-3)].translate(table=11)
                aln = smith_waterman(query_seq_segment, subject_seq_aa)
                pval = aln["pval"]
                psim = aln["psim"]
                if pval < best_pval:
                    best_pval = pval
                    best_psim = psim
                    best_aln = aln
                    best_aln.update({"frame": frame})
                    best_start_pos_query = start_pos_query
                    best_start_pos_subject = start_pos_subject
                # If we have the same pval, chose the match with the highest similarity
                elif pval == best_pval:
                    if psim > best_psim:
                        best_psim = psim
                        best_aln = aln
                        best_aln.update({"frame": frame})
                        best_start_pos_query = start_pos_query
                        best_start_pos_subject = start_pos_subject
                    # If we have the same similarity, chose the match with the highest length
                    elif psim == best_psim:
                        if aln["length"] > best_length:
                            best_length = aln["length"]
                            best_aln = aln
                            best_aln.update({"frame": frame})
                            best_start_pos_query = start_pos_query
                            best_start_pos_subject = start_pos_subject
    # Check the best alignment is qualitative
    if is_significant(best_aln):
        # Replace the relative positions by absolute ones
        best_aln = get_absolute_match(best_aln, best_start_pos_query, best_start_pos_subject, False)
        # Save the alignment
        matches.append(best_aln)
        # Get the segments that aren't covered by a match
        qstarts, qends, sstarts, sends = get_uncovered_segments(matches, len(query_seq), len(subject_seq_nu))
        recursively_align(query_seq, subject_seq_nu, qstarts, qends, sstarts, sends, matches)


def order_matches(matches):
    qstarts = np.array([m["qstart"] for m in matches])
    sorted_index = np.argsort(qstarts)
    sorted_matches = []
    for i in sorted_index:
        sorted_matches.append(matches[i])
    return sorted_matches


def get_significant_blastp(query, subject, frame):
    hits = []
    # Create temp files
    with tempfile.NamedTemporaryFile(mode='w', delete=False) as query_file:
        query_file.write(f">query\n{query}\n")
        query_file_path = query_file.name
    with tempfile.NamedTemporaryFile(mode='w', delete=False) as subject_file:
        subject_file.write(f">subject\n{subject}\n")
        subject_file_path = subject_file.name
    with tempfile.NamedTemporaryFile(mode='w', delete=False) as output_file:
        output_file_path = output_file.name
    with tempfile.NamedTemporaryFile(mode='w', delete=False) as output_file_raw:
        output_file_path_raw = output_file_raw.name
    # Run Blast
    out_format = "6 qseqid sseqid evalue qcovhsp qstart qend sstart send length"
    subprocess.run(["blastp", "-query", query_file_path, "-subject", subject_file_path, "-out", output_file_path, "-outfmt", out_format, "-evalue", "1e-3"], capture_output=True, check=True)
    subprocess.run(["blastp", "-query", query_file_path, "-subject", subject_file_path, "-out", output_file_path_raw, "-outfmt", "0", "-evalue", "1e-3"], capture_output=True, check=True)
    raw = subprocess.run(f"grep -A 2 'Query  ' {output_file_path_raw}", capture_output=True, text=True, shell=True)
    # Check the output isn't blank
    with open(output_file_path, "r") as f:
        content = f.read().strip()
    if content == "":
        os.remove(query_file_path)
        os.remove(subject_file_path)
        os.remove(output_file_path)
        os.remove(output_file_path_raw)
        return []
    # Read the output
    result = pd.read_csv(output_file_path, sep="\t", header=None)
    result.columns = ["qseqid", "sseqid", "evalue", "qcov", "qstart", "qend", "sstart", "send", "length"]
    for row in result.iterrows():
        qstart = int(row[1]["qstart"]) - 1
        qend = int(row[1]["qend"])
        sstart = (int(row[1]["sstart"]) - 1) * 3
        send = int(row[1]["send"]) * 3
        length = int(row[1]["length"])
        eval = int(row[1]["evalue"])
        if length >= 5:
            hits.append({"qstart": qstart, "qend": qend, "sstart": sstart, "send": send, "length": length, "eval": eval, "frame": frame, "raw": raw.stdout})
    # Delete temp files
    os.remove(query_file_path)
    os.remove(subject_file_path)
    os.remove(output_file_path)
    os.remove(output_file_path_raw)
    return hits


def look_for_frameshifts(denovo_seq, denovo_start, denovo_end, extended_match_seq, extended_start, extended_end, use_blast=False):
    if use_blast:
        matches = []
        for frame in [0, 1, 2]:
            extended_match_seq_aa = extended_match_seq[frame:(frame-3)].translate(table=11)
            matches += get_significant_blastp(denovo_seq, extended_match_seq_aa, frame)
        return matches

    extend_left = extended_start != 0
    extend_right = extended_end != len(extended_match_seq)
    matches_left, matches_right = [], []
    # Start with the left
    if extend_left:
        left_denovo = denovo_seq[:denovo_start]
        extended_left_seq = extended_match_seq[:extended_start]
        # Get the matches recursively
        recursively_align(left_denovo, extended_left_seq, [0], [len(left_denovo)], [0], [len(extended_left_seq)], matches_left)
    # Do the right
    if extend_right:
        right_denovo = denovo_seq[denovo_end:]
        extended_right_seq = extended_match_seq[extended_end:]
        # Get the matches recursively
        recursively_align(right_denovo, extended_right_seq, [0], [len(right_denovo)], [0], [len(extended_right_seq)], matches_right)
        matches_right = [get_absolute_match(r, denovo_end, extended_end, True) for r in matches_right]
    return matches_left + matches_right


def grep_after(text, pattern, num_lines=2):
    """Mimic grep -A 2"""
    lines = text.split('\n')
    result = []
    i = 0
    while i < len(lines):
        if pattern in lines[i]:
            result.extend(lines[i:i + num_lines + 1])
            result.append('--') 
            i += num_lines 
        else:
            i += 1
    return '\n'.join(result)


def print_results(denovo, all_matches_recursive, all_matches_blast, qcov_rec, qcov_blast, real_scale=False):
    # Get the frame covers in string form
    min_start_blast = int(min([dic["sstart"] for dic in all_matches_blast]) / 3)
    max_end_blast = int(max([dic["send"] for dic in all_matches_blast]) / 3)
    min_start_recursive = int(min([dic["sstart"] for dic in all_matches_recursive]) / 3)
    max_end_recursive = int(max([dic["send"] for dic in all_matches_recursive]) / 3)
    min_start = int(min(min_start_blast, min_start_recursive))
    max_end = int(max(max_end_blast, max_end_recursive))

    # Print result for recursive algorithm
    strings = {0: [" "] * (min_start_recursive - min_start) + ["-"] * (max_end_recursive - min_start_recursive) + [" "] * (max_end - max_end_recursive),
               1: [" "] * (min_start_recursive - min_start) + ["-"] * (max_end_recursive - min_start_recursive) + [" "] * (max_end - max_end_recursive),
               2: [" "] * (min_start_recursive - min_start) + ["-"] * (max_end_recursive - min_start_recursive) + [" "] * (max_end - max_end_recursive)}
    for match in all_matches_recursive:
        frame = match["frame"]
        str_start = int((match["sstart"] / 3 - min_start))
        str_end = int((match["send"] / 3 - min_start))
        for i in range(str_start, str_end):
            strings[frame][i] = "*"
    print(f"{''.join(strings[0])}\t{denovo}\n")
    print(f"{''.join(strings[1])}\tqcov = {qcov_rec}%\n")
    print(f"{''.join(strings[2])}\tsmith-waterman\n\n")

    for match in all_matches_recursive:
        print(grep_after(match["raw"], "query  "))
    print("\n\n\n\n")


    # Print result for blast algorithm
    strings = {0: [" "] * (min_start_blast - min_start) + ["-"] * (max_end_blast - min_start_blast) + [" "] * (max_end - max_end_blast),
               1: [" "] * (min_start_blast - min_start) + ["-"] * (max_end_blast - min_start_blast) + [" "] * (max_end - max_end_blast),
               2: [" "] * (min_start_blast - min_start) + ["-"] * (max_end_blast - min_start_blast) + [" "] * (max_end - max_end_blast)}
    for match in all_matches_blast:
        frame = match["frame"]
        str_start = int((match["sstart"] / 3 - min_start))
        str_end = int((match["send"] / 3 - min_start))
        for i in range(str_start, str_end):
            strings[frame][i] = "*"
    print(f"{''.join(strings[0])}\t{denovo}\n")
    print(f"{''.join(strings[1])}\tqcov = {qcov_blast}%\n")
    print(f"{''.join(strings[2])}\tblast\n\n")

    for match in all_matches_blast:
        print(match["raw"])


def read_good_candidates():
    good_candidates = []
    with open(GOOD_CANDIDATES_FILE, "r") as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith("#"):
                good_candidates.append(line)
    return good_candidates


def look_for_frameshifts_within(denovo_seq, qstart, qend, extended_match_seq, sstart, send):
    matches = []
    recursively_align(denovo_seq, extended_match_seq, [qstart], [qend], [sstart], [send], matches)
    return(matches)








if __name__ == "__main__":
    with open(GENOMES_LIST, "r") as f:
        genomes = f.readline().split()
    genomes = [re.sub('"', '', g) for g in genomes]

    # Get the denovo info
    denovo_dict = {}
    for genome in genomes:
        new_denovo = extract_denovo_info(genome, OUTGROUP_NUMBER)

        # For each denovo gene
        for denovo in new_denovo:
            """#--------------------------------------------------
            if "HPMEPLIM_01176_gene_mRNA" not in denovo:
                continue
            #--------------------------------------------------"""
            new_denovo[denovo]["genome"] = genome

            # Get the match sequence
            gene_len = len(new_denovo[denovo]["sequence"])
            loci = new_denovo[denovo]["loci"]
            contig = loci[0]
            start = loci[1]
            end = loci[2]
            strand = loci[3]
            outgroup = new_denovo[denovo]["ancestor_sp"]
            # Look if the Cter has been entirely matched
            missing_cter = gene_len - new_denovo[denovo]["qend"]
            # Look if the Nter has been entirely matched
            missing_nter = new_denovo[denovo]["qstart"]
            # Get the extended match sequence in outgroup
            match_seq_extended, extended_start, extended_end = get_extended_matched_seq(outgroup, contig, start, end, strand, missing_cter, missing_nter)
            match_seq = get_sequence_from_loci(outgroup, contig, start, end, strand)
            # Add the info to the dict
            new_denovo[denovo]["extended_match_seq"] = match_seq_extended
            new_denovo[denovo]["extended_start"] = extended_start
            new_denovo[denovo]["extended_end"] = extended_end
            new_denovo[denovo]["match_seq"] = match_seq


        # Add to global dict
        denovo_dict.update(new_denovo)


    """#--------------------------------------------------------------
    # Keep only gene of interest
    dict_interest = {}
    dict_interest["HPMEPLIM_01176_gene_mRNA"] = denovo_dict["HPMEPLIM_01176_gene_mRNA"]
    denovo_dict = dict_interest
    #--------------------------------------------------------------"""
    print(f"OUTGROUP {OUTGROUP_NUMBER}\n\n\n")
    
    # Get the list of good denovo candidates
    good_candidates = read_good_candidates()

    # Iterate over the "bad" denovo candidates
    for denovo in denovo_dict:
        if denovo in good_candidates:
            continue

        # Look for frameshifts within the match
        denovo_seq = denovo_dict[denovo]["sequence"]
        qstart = denovo_dict[denovo]["qstart"]
        qend = denovo_dict[denovo]["qend"]
        extended_match_seq = denovo_dict[denovo]["extended_match_seq"]
        match_seq = denovo_dict[denovo]["match_seq"]

        frameshifts = look_for_frameshifts_within(denovo_seq, qstart, qend, match_seq, 0, len(match_seq))

        frames = [aln["frame"] for aln in frameshifts]
        # For all frameshifts / altframes
        if len(frameshifts) > 1 or set(frames) != {0}:
            max_qcov = 0
            # Print the results
            print(f"{denovo}\n\n")
            for aln in frameshifts:
                # Print the alignment info
                print(f"Query: {denovo_seq} (length {len(denovo_seq)})")
                print("Alignment infos:")
                print([f"{key}: {value}" for key, value in aln.items() if key != "raw\n"])
                # Get the maximum qcov
                qstart = aln["qstart"]
                qend = aln["qend"]
                query_len = len(denovo_seq)
                qcov_aln = (qend - qstart) / query_len * 100
                if qcov_aln > max_qcov:
                    max_qcov = qcov_aln
                print()
                # Print the alignment raw
                print(aln["raw"])
            print(f"\n\nThe maximum qcov is {max_qcov}%\n")
            print("\n\n\n\n")
            print("______________________________________________________________________________\n\n\n")
        










    """for denovo in denovo_dict:
        extended_match_seq = denovo_dict[denovo]["extended_match_seq"]
        extended_start = denovo_dict[denovo]["extended_start"]
        extended_end = denovo_dict[denovo]["extended_end"]
        denovo_seq = denovo_dict[denovo]["sequence"]
        denovo_start = denovo_dict[denovo]["qstart"]
        denovo_end = denovo_dict[denovo]["qend"]

        # Look for frameshift on both sides
        frameshifts_recursive = look_for_frameshifts(denovo_seq, denovo_start, denovo_end, extended_match_seq, extended_start, extended_end, use_blast=False)
        all_matches_blast = order_matches(look_for_frameshifts(denovo_seq, denovo_start, denovo_end, extended_match_seq, extended_start, extended_end, use_blast=True))

        # Get list of all matches for all frames
        origin_match = {"qstart": denovo_start, "qend": denovo_end, "sstart": extended_start, "send": extended_end, "frame": 0, "raw": ""}
        all_matches_recursive = order_matches(frameshifts_recursive + [origin_match])

        # Don't continue if no match
        if len(all_matches_recursive) == 1 and len(all_matches_blast) == 1:
            continue

        # Get the total qcov
        bases_covered_recursive = []
        for dic in all_matches_recursive:
            bases_covered_recursive += list(range(dic["qstart"], dic["qend"]))
        total_qcov_rec = round((len(set(bases_covered_recursive)) / len(denovo_seq) * 100), 1)
        bases_covered_blast = []
        for dic in all_matches_blast:
            bases_covered_blast += list(range(dic["qstart"], dic["qend"]))
        total_qcov_blast = round((len(set(bases_covered_blast)) / len(denovo_seq) * 100), 1)

        # Print the results
        print(f"{denovo}\n")
        print(f"With water & recursion: found {len(all_matches_recursive)-1} significant alignments\n")
        print(f"With tblastn: found {len(all_matches_blast)-1} significant alignments\n\n\n")
        print_results(denovo, all_matches_recursive, all_matches_blast, total_qcov_rec, total_qcov_blast, real_scale=True)
        print(f"\n\n{"_" * 300}\n\n")"""