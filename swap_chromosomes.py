#!/usr/bin/env python3
import argparse
import csv
import os
import re
import sys
import math
from pathlib import Path
from random import Random
from typing import Dict, List, Tuple, Optional, Union

# Try to enable headless plotting
try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.ticker as mticker
except Exception:
    plt = None

# ---------- Utilities ----------

# Replace these in your script
SAMPLE_TAG_RE = re.compile(r"(1604-\d+(?:-\d+){1,3}-5x[A-Za-z0-9]+)")

def extract_sample_tag(filename: str) -> str:
    """
    Return a tag that ALWAYS starts with '1604-...'.
    Works with names that have leading prefixes (e.g., 'A64-1604-...') and
    normalizes underscores to hyphens.
    """
    name = Path(filename).stem.replace("_", "-")
    m = SAMPLE_TAG_RE.search(name)
    if m:
        return m.group(1)
    i = name.find("1604-")
    if i != -1:
        return name[i:]
    return name  # last resort

def list_tsv_files(d: Path) -> List[Path]:
    return sorted([p for p in d.iterdir() if p.is_file() and p.suffix.lower() == ".tsv"])

def detect_header_third_is_chr(row: List[str]) -> bool:
    return len(row) >= 3 and row[2].strip().lower() == "chr"

def read_group_by_chr(tsv_path: Path) -> Tuple[Optional[List[str]], Dict[str, List[List[str]]], List[str]]:
    """
    Read a TSV and group rows by the value in the 3rd column (1-based).
    Returns (header, groups, chr_order).
    groups contains ONLY data rows (header removed if detected).
    """
    groups: Dict[str, List[List[str]]] = {}
    chr_order: List[str] = []
    header: Optional[List[str]] = None
    with tsv_path.open(newline="") as f:
        reader = csv.reader(f, delimiter="\t")
        first = True
        for row in reader:
            if not row:
                continue
            if first:
                first = False
                if detect_header_third_is_chr(row):
                    header = row
                    continue
            if len(row) < 3:
                continue
            chr_val = row[2].strip()
            if chr_val not in groups:
                groups[chr_val] = []
                chr_order.append(chr_val)
            groups[chr_val].append(row)
    return header, groups, chr_order

def natural_chr_key(c: str) -> Tuple[int, str]:
    """
    Sort key: chr1..chr22 (1..22) get numeric ordering; all others (e.g., chrX, chrY, chrM)
    sort after those by name.
    """
    c_low = c.lower()
    if c_low.startswith("chr"):
        tail = c_low[3:]
        if tail.isdigit():
            return (int(tail), "")
    return (9999, c_low)

def ensure_outdir(p: Path):
    p.mkdir(parents=True, exist_ok=True)

def load_chr_sizes(path: Path) -> Dict[str, int]:
    """
    Read a whitespace/tab-delimited file with at least two columns:
      <chrom> <length>
    Returns dict chrom -> int(length). Ignores blanks and lines starting with '#'.
    """
    sizes: Dict[str, int] = {}
    try:
        with path.open("r") as fh:
            for line in fh:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                parts = line.split()
                if len(parts) < 2:
                    continue
                chrom, length_str = parts[0], parts[1]
                try:
                    sizes[chrom] = int(length_str)
                except ValueError:
                    continue
    except FileNotFoundError:
        print(f"[WARN] Chrom sizes file not found: {path} — swapped length will be 0.", file=sys.stderr)
    return sizes

def write_swapped(dipl_path: Path,
                  trip_path: Path,
                  selected_chroms: List[str],
                  out_path: Path) -> Tuple[int, List[str], List[str], int, int, int, int]:
    """
    Perform the swap and write out_path.

    Returns:
      written_rows: int (data rows written to output)
      used_selected: List[str] (selected chromosomes present in BOTH dip & trip)
      missing_selected_in_trip: List[str] (selected chromosomes missing from TRIP)
      dip_total_rows: int (total data rows in dip file)
      trip_total_rows: int (total data rows in trip file)
      dip_rows_to_swap: int (sum of dip rows for used_selected chrs)
      trip_rows_to_swap: int (sum of trip rows for used_selected chrs)
    """
    d_header, d_groups, d_order = read_group_by_chr(dipl_path)
    t_header, t_groups, t_order = read_group_by_chr(trip_path)

    selected = list(selected_chroms)  # keep order
    missing_in_trip = [c for c in selected if c not in t_groups]
    used = [c for c in selected if (c in t_groups) and (c in d_groups)]

    dip_total_rows = sum(len(v) for v in d_groups.values())
    trip_total_rows = sum(len(v) for v in t_groups.values())
    dip_rows_to_swap = sum(len(d_groups[c]) for c in used)
    trip_rows_to_swap = sum(len(t_groups[c]) for c in used)

    with out_path.open("w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        if d_header is not None:
            writer.writerow(d_header)
        written = 0
        for c in d_order:
            use_trip = c in used
            rows = t_groups[c] if use_trip else d_groups.get(c, [])
            for r in rows:
                writer.writerow(r)
                written += 1

    return written, used, missing_in_trip, dip_total_rows, trip_total_rows, dip_rows_to_swap, trip_rows_to_swap

def pick_chroms_for_pair(rng: Random,
                         override_chroms: List[str],
                         n_chroms: int,
                         dipl_path: Path,
                         trip_path: Path) -> List[str]:
    """
    If override_chroms provided, use them (skip those not present).
    Else, sample n_chroms from intersection of chromosomes in dip & trip.
    """
    _, d_groups, _ = read_group_by_chr(dipl_path)
    _, t_groups, _ = read_group_by_chr(trip_path)
    d_chrs = set(d_groups.keys())
    t_chrs = set(t_groups.keys())
    available = sorted(d_chrs & t_chrs, key=natural_chr_key)
    if override_chroms:
        return [c for c in override_chroms if c in available]
    k = min(n_chroms, len(available))
    if k <= 0:
        return []
    choices = available[:]
    rng.shuffle(choices)
    return sorted(choices[:k], key=natural_chr_key)

# ---------- Histogram helpers ----------

def percentile(sorted_vals: List[float], p: float) -> float:
    if not sorted_vals:
        return 0.0
    if len(sorted_vals) == 1:
        return float(sorted_vals[0])
    k = (len(sorted_vals) - 1) * (p / 100.0)
    f = math.floor(k); c = math.ceil(k)
    if f == c:
        return float(sorted_vals[int(k)])
    d0 = float(sorted_vals[f]) * (c - k)
    d1 = float(sorted_vals[c]) * (k - f)
    return d0 + d1

def compute_bins(values: List[int], strategy: Union[str, int, None]) -> Tuple[List[float], List[int]]:
    """
    Compute histogram bins and counts in base pairs.
    Returns (bin_edges_bp, counts) where len(edges)=len(counts)+1.
    """
    n = len(values)
    if n == 0:
        return [0.0, 1.0], [0]
    vmin = float(min(values)); vmax = float(max(values))
    if vmin == vmax:
        return [vmin - 0.5, vmax + 0.5], [n]
    # number of bins
    if isinstance(strategy, int):
        nbins = max(1, strategy)
    else:
        s = (strategy or "fd").lower()
        if s in ("fd", "freedman-diaconis"):
            sv = sorted(values); q1 = percentile(sv, 25); q3 = percentile(sv, 75); iqr = max(0.0, q3 - q1)
            if iqr <= 0.0:
                nbins = max(1, int(math.ceil(math.log2(n) + 1)))
            else:
                bin_width = 2.0 * iqr * (n ** (-1.0/3.0))
                nbins = max(1, int(math.ceil((vmax - vmin) / bin_width)))
        elif s in ("sturges", "s"):
            nbins = max(1, int(math.ceil(math.log2(n) + 1)))
        elif s in ("auto",):
            if n <= 2:
                nbins = 1
            else:
                sv = sorted(values); q1 = percentile(sv, 25); q3 = percentile(sv, 75); iqr = max(0.0, q3 - q1)
                if iqr <= 0.0:
                    nbins = max(1, int(math.ceil(math.log2(n) + 1)))
                else:
                    bin_width = 2.0 * iqr * (n ** (-1.0/3.0))
                    nbins = max(1, int(math.ceil((vmax - vmin) / bin_width)))
        else:
            nbins = max(1, int(math.ceil(math.log2(n) + 1)))
    edges = [vmin + (vmax - vmin) * i / nbins for i in range(nbins + 1)]
    counts = [0] * nbins
    for v in values:
        idx = nbins - 1 if v == vmax else int((v - vmin) / (vmax - vmin) * nbins)
        counts[min(max(idx, 0), nbins - 1)] += 1
    return edges, counts

def write_hist_tsv(edges: List[float], counts: List[int], out_tsv: Path):
    with out_tsv.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["bin_left_bp", "bin_right_bp", "count"])
        for i, c in enumerate(counts):
            w.writerow([int(edges[i]), int(edges[i+1]), c])

def plot_histogram(edges_bp: List[float], counts: List[int], out_png: Path, title: str):
    """
    Plot histogram with x-axis in megabases (scientific notation).
    """
    if plt is None:
        print("[WARN] matplotlib not available; skipping PNG plot.", file=sys.stderr)
        return
    edges_mb = [e / 1e6 for e in edges_bp]
    widths_mb = [edges_mb[i+1] - edges_mb[i] for i in range(len(counts))]
    centers_mb = [edges_mb[i] for i in range(len(counts))]
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.bar(centers_mb, counts, width=widths_mb, align="edge")
    ax.set_xlabel("Total Trisomy Size (Mb)")
    ax.set_ylabel("count")
    ax.set_title(title)
    sf = mticker.ScalarFormatter(useMathText=True)
    sf.set_scientific(True)
    sf.set_powerlimits((0, 0))  # always sci notation
    ax.xaxis.set_major_formatter(sf)
    fig.tight_layout()
    fig.savefig(out_png, dpi=150)
    plt.close(fig)

def read_swapped_bp_from_log(log_path: Path) -> List[int]:
    """
    Read swapped_chr_total_bp column from swap_runs.log (tab-delimited).
    Returns list of ints. Skips header automatically.
    """
    values: List[int] = []
    with log_path.open("r", newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader, None)
        idx = None
        if header:
            for i, col in enumerate(header):
                if col.strip() == "swapped_chr_total_bp":
                    idx = i
                    break
        for row in reader:
            if not row:
                continue
            try:
                values.append(int(row[idx if idx is not None else 4]))
            except (ValueError, IndexError, TypeError):
                continue
    return values

# ---------- Main workflow ----------

def main():
    parser = argparse.ArgumentParser(
        description="Swap selected chromosome rows (3rd column == 'chr') "
                    "from TRIP tsvs into DIP tsvs. Outputs TSVs, a run log, and an optional histogram."
    )
    parser.add_argument("--dipl-dir", default="/tbusa/scratch/kdistor/sequencing_Data/simul_trip_downsampled/dipl",
                        type=Path, help="Directory containing diploid TSV files")
    parser.add_argument("--trip-dir", default="/tbusa/scratch/kdistor/sequencing_Data/simul_trip_downsampled/trip",
                        type=Path, help="Directory containing triploid TSV files")
    parser.add_argument("--out-dir", default="swapped_tsvs",
                        type=Path, help="Directory to write outputs and run log")

    # Chromosome selection
    parser.add_argument("--chroms", nargs="+",
                        help="Explicit list of chromosomes to swap (e.g. --chroms chr1 chr17 chrX). "
                             "If omitted, --n-chroms random chromosomes will be chosen per pair.")
    parser.add_argument("--n-chroms", type=int, default=3,
                        help="When --chroms is not provided, choose this many chromosomes per pair (default: 3)")

    # Pairing control
    parser.add_argument("--trip-uses", type=int, default=2,
                        help="Number of times each triploid file will be used (default: 2)")

    # Seeds
    parser.add_argument("--pair-seed", type=int, default=100,
                        help="Random seed for pairing trip->dip (default: 100)")
    parser.add_argument("--chrom-seed", type=int, default=100,
                        help="Random seed for random chromosome selection (default: 100)")

    # Chromosome sizes file (for totaling swapped bp)
    parser.add_argument("--chr-sizes",
        type=Path,
        default=Path("/tbusa/scratch/kdistor/sequencing_Data/simul_trip_downsampled/hg38_chr_sizes_1_22_XY"),
        help="Two-column chrom sizes file (e.g., 'chr1<TAB>248956422').")

    # Histogram/plot controls
    parser.add_argument("--hist-bins",
        default="fd",
        help="Histogram binning: integer number of bins or one of {'fd','sturges','auto'}. Default: fd.")
    parser.add_argument("--no-plot",
        action="store_true",
        help="Disable creation of swapped bp histogram plot/tsv.")

    args = parser.parse_args()

    # Prep
    ensure_outdir(args.out_dir)
    chr_sizes = load_chr_sizes(args.chr_sizes)
    if not chr_sizes:
        print(f"[WARN] No chromosome sizes loaded from {args.chr_sizes}. Sums will be 0.", file=sys.stderr)

    log_path = args.out_dir / "swap_runs.log"
    with log_path.open("w", newline="") as logf:
        log_writer = csv.writer(logf, delimiter="\t")
        log_writer.writerow([
            "run_number",
            "triploid_file",
            "diploid_file",
            "chromosomes_swapped",
            "swapped_chr_total_bp",
            "dip_rows_to_swap",
            "trip_rows_to_swap",
            "dip_total_rows",
            "trip_total_rows",
            "swapped_file_rows",
            "output_file"
        ])

        # Gather files
        dipl_files = list_tsv_files(args.dipl_dir)
        trip_files = list_tsv_files(args.trip_dir)

        if not trip_files:
            print(f"[ERROR] No TSVs found in trip dir: {args.trip_dir}", file=sys.stderr)
            sys.exit(1)
        if not dipl_files:
            print(f"[ERROR] No TSVs found in dipl dir: {args.dipl_dir}", file=sys.stderr)
            sys.exit(1)

        print(f"[INFO] Found {len(trip_files)} triploid files and {len(dipl_files)} diploid files.")
        print(f"[INFO] Each triploid file will be used {args.trip_uses} time(s).")

        # Deterministic pairing setup
        pair_rng = Random(args.pair_seed)
        chrom_rng_base = args.chrom_seed

        # Pre-shuffle diploid list once; for each trip we pick a random starting offset
        dipls_shuffled = dipl_files[:]
        pair_rng.shuffle(dipls_shuffled)

        run_counter = 0

        for _, trip_path in enumerate(trip_files):
            start_offset = pair_rng.randrange(len(dipls_shuffled))
            for k in range(args.trip_uses):
                dip_path = dipls_shuffled[(start_offset + k) % len(dipls_shuffled)]
                run_counter += 1

                # Choose chromosomes
                chrom_rng = Random(chrom_rng_base + run_counter)
                override = args.chroms if args.chroms else []
                chosen = pick_chroms_for_pair(chrom_rng, override, args.n_chroms, dip_path, trip_path)
                chosen_sorted = sorted(chosen, key=natural_chr_key)

                # Build output name
                dip_tag = extract_sample_tag(dip_path.name)
                trip_tag = extract_sample_tag(trip_path.name)
                chrom_suffix = "_".join(chosen_sorted) if chosen_sorted else "nochr"
                out_name = f"{dip_tag}-{trip_tag}_{chrom_suffix}_run{run_counter:03d}.tsv"
                out_path = args.out_dir / out_name

                print(f"[RUN {run_counter:03d}] Trip: {trip_path.name}")
                print(f"[RUN {run_counter:03d}] Dip : {dip_path.name}")
                print(f"[RUN {run_counter:03d}] Chr : {', '.join(chosen_sorted) if chosen_sorted else '(none)'}")
                print(f"[RUN {run_counter:03d}] Out : {out_name}")

                # Do the swap + counts
                written, used_chroms, missing_in_trip, \
                    dip_total_rows, trip_total_rows, \
                    dip_rows_to_swap, trip_rows_to_swap = write_swapped(dip_path, trip_path, chosen_sorted, out_path)

                if missing_in_trip:
                    print(f"[WARN] Missing in trip ({trip_path.name}): {', '.join(missing_in_trip)} "
                          f"-> kept diploid rows for those.", file=sys.stderr)

                # Sum only the chromosomes that were actually swapped (present in both dip & trip)
                swapped_bp = sum(chr_sizes.get(c, 0) for c in used_chroms)
                if any(c not in chr_sizes for c in used_chroms):
                    unknown = [c for c in used_chroms if c not in chr_sizes]
                    print(f"[WARN] No size for {', '.join(unknown)} in {args.chr_sizes}; counted as 0.", file=sys.stderr)

                print(f"[RUN {run_counter:03d}] Swapped bp: {swapped_bp}")
                print(f"[RUN {run_counter:03d}] Rows -> dip_to_swap={dip_rows_to_swap}, trip_to_swap={trip_rows_to_swap}, "
                      f"dip_total={dip_total_rows}, trip_total={trip_total_rows}, swapped_file={written}")

                # Log
                log_writer.writerow([
                    f"{run_counter:03d}",
                    str(trip_path),
                    str(dip_path),
                    ",".join(chosen_sorted),
                    str(swapped_bp),
                    str(dip_rows_to_swap),
                    str(trip_rows_to_swap),
                    str(dip_total_rows),
                    str(trip_total_rows),
                    str(written),
                    str(out_path)
                ])

    # ---------- Build histogram from swap_runs.log ----------
    if not args.no_plot:
        try:
            swapped_values = read_swapped_bp_from_log(log_path)
            try:
                bins_arg: Union[int, str, None] = int(args.hist_bins)
            except (TypeError, ValueError):
                bins_arg = args.hist_bins
            edges_bp, counts = compute_bins(swapped_values, bins_arg)

            hist_tsv = args.out_dir / "swapped_bp_hist.tsv"
            write_hist_tsv(edges_bp, counts, hist_tsv)

            hist_png = args.out_dir / "swapped_bp_hist.png"
            plot_histogram(edges_bp, counts, hist_png, "Counts vs Total Trisomy Size (Mb)")

            print(f"[PLOT] Histogram table: {hist_tsv}")
            if plt is not None:
                print(f"[PLOT] Histogram image: {hist_png}")
            else:
                print(f"[PLOT] Skipped image (matplotlib not available).", file=sys.stderr)
        except Exception as e:
            print(f"[WARN] Failed to build histogram: {e}", file=sys.stderr)

    print(f"[DONE] Wrote outputs to: {args.out_dir}")
    print(f"[DONE] Run log saved to: {log_path}")

if __name__ == "__main__":
    main()
