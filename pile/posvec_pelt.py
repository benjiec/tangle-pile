import argparse
import re
import sys
import numpy as np
import ruptures as rpt


def parse_mpileup_alignment_tracks(filepath, target_tx):
    """Parses the raw text file to extract feature tracks and explicitly

    injects missing coordinates to preserve spatial scale across structural gaps.
    """
    ins_regex = re.compile(r"\+[0-9]+[ACGTNacgtn]+")
    del_regex = re.compile(r"-[0-9]+[ACGTNacgtn]+")
    snp_regex = re.compile(r"[ACGTNacgtn]")

    positions = []
    features = []
    last_pos = None

    """
    print(
        f"[*] Extracting positions and detecting gaps for {target_tx}...",
        file=sys.stderr,
    )
    """

    with open(filepath, "r") as f:
        for line in f:
            cols = line.strip().split("\t")
            if len(cols) < 6 or cols[0] != target_tx:
                continue

            pos = int(cols[1])
            depth = int(cols[3])
            pileup_str = cols[4]

            # Gap detection and explicit zero injection
            if last_pos is not None and (pos - last_pos) > 1:
                for missing_pos in range(last_pos + 1, pos):
                    positions.append(missing_pos)
                    features.append([0.0, 0.0, 0.0, 0.0])

            last_pos = pos

            if depth == 0 or pileup_str == "*":
                positions.append(pos)
                features.append([0.0, 0.0, 0.0, 0.0])
                continue

            # Clean read start/end markers
            clean_str = re.sub(r"\^.", "", pileup_str).replace("$", "")

            ins_count = len(ins_regex.findall(clean_str))
            del_count = len(del_regex.findall(clean_str))

            clean_str = ins_regex.sub("", clean_str)
            clean_str = del_regex.sub("", clean_str)
            snp_count = len(snp_regex.findall(clean_str))

            clip_count = pileup_str.count("^") + pileup_str.count("$")

            positions.append(pos)
            features.append(
                [
                    float(depth),
                    snp_count / depth,
                    (ins_count + del_count) / depth,
                    clip_count / depth,
                ]
            )

    if not features:
        raise ValueError(f"Transcript '{target_tx}' not found in file.")

    return np.array(positions), np.array(features)


def standardize_signal_variance(feature_matrix):
    """Standardizes variance across multi-scale channels using log and Z-score transformations."""

    """
    print("[*] Normalizing multi-channel signal tracks...", file=sys.stderr)
    """

    norm_matrix = np.copy(feature_matrix)

    # Log-transform depth channel
    norm_matrix[:, 0] = np.log2(norm_matrix[:, 0] + 1)

    # Z-score normalization across all 4 channels
    means = np.mean(norm_matrix, axis=0)
    stds = np.std(norm_matrix, axis=0) + 1e-6
    return (norm_matrix - means) / stds


def detect_signal_changepoints(normalized_matrix, positions, penalty, min_size):
    """Runs PELT change-point detection to identify segment array indexes."""

    """
    print(
        f"[*] Computing PELT segmentation (Penalty={penalty}, Min_Size={min_size})...",
        file=sys.stderr,
    )
    """

    algo = rpt.KernelCPD(kernel="rbf", min_size=min_size).fit(normalized_matrix)
    return algo.predict(pen=penalty)


def compile_segment_summaries(
    feature_matrix, raw_positions, result_indexes, global_mean_depth
):
    """Aggregates data within windows, calculates text classification tags,

    and generates the structured CIGAR-like summary records.
    """

    """
    print("[*] Generating categorical segment summaries...", file=sys.stderr)
    """

    segment_records = []
    start_idx = 0

    for i, idx in enumerate(result_indexes):
        idx_adjusted = min(idx, len(raw_positions) - 1)

        start_pos = raw_positions[start_idx]
        end_pos = raw_positions[idx_adjusted]
        length = idx_adjusted - start_idx + 1

        # Slice features for this segment window
        seg_data = feature_matrix[start_idx : idx_adjusted + 1]
        m_dp = np.mean(seg_data[:, 0])
        m_snp = np.mean(seg_data[:, 1])
        m_indel = np.mean(seg_data[:, 2])
        m_clip = np.mean(seg_data[:, 3])

        # Classification Engine logic
        tags = []
        if m_dp == 0:
            tags.append("Z")
        else:
            if m_dp >= 1.5 * global_mean_depth:
                tags.append("H")
            elif m_dp <= 0.5 * global_mean_depth:
                tags.append("L")

            if m_snp >= 0.15:
                tags.append("V")
            if m_indel >= 0.10:
                tags.append("I")
            if m_clip >= 0.30:
                tags.append("C")

        tag_str = "".join(tags) if tags else "M"

        record = {
            "id": i + 1,
            "start": start_pos,
            "end": end_pos,
            "len": length,
            "tag": tag_str,
            "mean_depth": m_dp,
            "mean_snp": m_snp,
            "mean_indel": m_indel,
            "mean_clip": m_clip,
        }
        segment_records.append(record)
        start_idx = idx

    return segment_records


def generate_summary_string(transcript_id, raw_positions, global_mean_depth, segments):
    seg_summaries = [
        f"{seg['start']}-{seg['end']}:{seg['tag']}" for seg in segments
    ]
    return ",".join(seg_summaries)


def print_report(transcript_id, raw_positions, global_mean_depth, segments):
    """Console print formatter for the compiled segmentation report."""
    print(f"Transcript ID: {transcript_id}")
    print(f"Total Reference Span: {raw_positions[0]} - {raw_positions[-1]} bp")
    print(f"Global Base Mean Depth: {global_mean_depth:.2f}")
    print("-" * 105)
    print(
        f"{'Segment':<12}{'Start':<8}{'End':<8}{'Len':<6}{'Tag':<10}{'Mean_Dp':<9}{'Mean_SNP':<10}{'Mean_Indel':<12}{'Mean_Clip':<10}"
    )
    print("-" * 105)
    for seg in segments:
        print(
            f"Seg_{seg['id']:<4}\t"
            f"{seg['start']:<8}"
            f"{seg['end']:<8}"
            f"{seg['len']:<6}"
            f"{seg['tag']:<10}"
            f"{seg['mean_depth']:<9.1f}"
            f"{seg['mean_snp']:<10.3f}"
            f"{seg['mean_indel']:<12.3f}"
            f"{seg['mean_clip']:<10.3f}"
        )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Modular Multi-channel Transcript Segmentation"
    )
    parser.add_argument("-i", "--input", required=True, help="Mpileup file path")
    parser.add_argument("-t", "--transcript", required=True, help="Transcript ID")
    parser.add_argument("-p", "--penalty", type=float, default=5.0)
    parser.add_argument("-m", "--min_size", type=int, default=5)

    args = parser.parse_args()

    try:
        # Pipeline execution blocks
        raw_positions, feature_matrix = parse_mpileup_alignment_tracks(
            args.input, args.transcript
        )

        normalized_matrix = standardize_signal_variance(feature_matrix)

        result_indexes = detect_signal_changepoints(
            normalized_matrix, raw_positions, args.penalty, args.min_size
        )

        global_mean_depth = np.mean(feature_matrix[:, 0])
        segments = compile_segment_summaries(
            feature_matrix, raw_positions, result_indexes, global_mean_depth
        )

        # print_report(args.transcript, raw_positions, global_mean_depth, segments)
        print("%s\t%s" % (args.transcript, generate_summary_string(args.transcript, raw_positions, global_mean_depth, segments)))

    except Exception as e:
        print(f"[!] Execution failure: {e}", file=sys.stderr)
        sys.exit(1)
