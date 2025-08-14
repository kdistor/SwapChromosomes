# swap\_chromosomes.py — README

Swap selected chromosome rows (3rd column == `chr`) from **TRIP** TSVs into **DIP** TSVs. The script writes new swapped TSVs, a detailed run log, and an optional histogram summarizing total trisomy sizes.

---

## Highlights

- **Deterministic pairing & selection** via seeds
- **Smart sample tags** in output filenames (always start with `1604-...`)
- **Header-aware** TSV parsing (keeps header if 3rd column equals `chr`)
- **Chrom-size totals** per run (`swapped_chr_total_bp`)
- **Row accounting** per run (rows to swap and totals for DIP/TRIP/output)
- **Histogram** of counts vs **Total Trisomy Size (Mb)** (PNG + TSV), x‑axis in Mb with scientific notation

---

## Requirements

- Python 3.8+
- Standard library only for core swap & logging
- Optional (for plotting): `matplotlib`
  - If not installed, the script still runs; it just skips the PNG plot

---

## Input Assumptions

- Inputs are **tab-delimited TSV** files.
- The **3rd column** (1‑based) contains chromosome names (e.g., `chr1`).
- A header row is detected **iff** the 3rd column literally equals `chr` (case-insensitive). That header is preserved in outputs.

### Chromosome Sizes File

Used to sum `swapped_chr_total_bp` for each run.

- Default: `/tbusa/scratch/kdistor/sequencing_Data/simul_trip_downsampled/hg38_chr_sizes_1_22_XY`
- Format: Two columns, whitespace/tab delimited: `<chrom> <length>` (e.g., `chr1 248956422`)
- Lines beginning with `#` or empty lines are ignored

---

## Determinism & Selection

- **Pairing**: Each TRIP file is paired with DIP files deterministically using `--pair-seed` (pre‑shuffle + per‑TRIP random start offset).
- **Chromosome choice**: If `--chroms` is not provided, random selection of `--n-chroms` from the **intersection** of chromosomes present in both files, using `--chrom-seed` + run number for reproducibility.
- **Sorting**: `chr1..chr22` are ordered numerically; everything else (e.g., `chrX`, `chrY`, `chrM`) sorts lexicographically after those.

---

## Outputs

- **Swapped TSVs** in `--out-dir` (default `swapped_tsvs`)
- **Run log**: `swap_runs.log` (tab‑delimited) in `--out-dir`
- **Histogram outputs** (unless `--no-plot`):
  - `swapped_bp_hist.tsv` (binned counts; edges in **bp**)
  - `swapped_bp_hist.png` (x‑axis **Total Trisomy Size (Mb)** in scientific notation; title **Counts vs Total Trisomy Size (Mb)**)

### Output Filename Pattern

```
<dip_tag>-<trip_tag>_<chrlist>_runNNN.tsv
```

- `dip_tag` / `trip_tag`: extracted from filenames; always start at the first `1604-...`; underscores normalized to hyphens
- `chrlist`: selected chromosomes in DIP/TRIP‑aware order (or `nochr`)

### Run Log Schema (`swap_runs.log`)

Columns (tab‑delimited):

1. `run_number` — zero‑padded index per run
2. `triploid_file` — path to TRIP TSV
3. `diploid_file` — path to DIP TSV
4. `chromosomes_swapped` — comma‑joined list requested/used for the run (filtered by availability)
5. `swapped_chr_total_bp` — sum of chromosome sizes for chromosomes **actually** swapped (present in both DIP & TRIP)
6. `dip_rows_to_swap` — number of DIP data rows replaced (across used chromosomes)
7. `trip_rows_to_swap` — number of TRIP data rows inserted (across used chromosomes)
8. `dip_total_rows` — total DIP data rows (header excluded)
9. `trip_total_rows` — total TRIP data rows (header excluded)
10. `swapped_file_rows` — total data rows written to the output (header excluded)
11. `output_file` — path to the newly written TSV

> **Note:** If a requested chromosome isn’t present in TRIP (or DIP), it’s not counted in `used` and contributes 0 to `swapped_chr_total_bp` and to the `*_rows_to_swap` numbers.

---

## Usage

```
python3 swap_chromosomes.py \
  --dipl-dir /tbusa/scratch/kdistor/sequencing_Data/simul_trip_downsampled/dipl \
  --trip-dir /tbusa/scratch/kdistor/sequencing_Data/simul_trip_downsampled/trip \
  --out-dir swapped_tsvs \
  [--chroms chr1 chr7 chrX] \
  [--n-chroms 3] \
  [--trip-uses 2] \
  [--pair-seed 100] \
  [--chrom-seed 100] \
  [--chr-sizes /path/to/chrom_sizes] \
  [--hist-bins fd|sturges|auto|<int>] \
  [--no-plot]
```

### Common Examples

**1) Minimal (defaults; 3 random chromosomes; each TRIP used twice)**

```
python3 swap_chromosomes.py \
  --dipl-dir /tbusa/scratch/kdistor/sequencing_Data/simul_trip_downsampled/dipl \
  --trip-dir /tbusa/scratch/kdistor/sequencing_Data/simul_trip_downsampled/trip \
  --out-dir swapped_tsvs
```

**2) Specific chromosomes; 3 uses per TRIP; Sturges bins**

```
python3 swap_chromosomes.py \
  --chroms chr1 chr7 chrX \
  --trip-uses 3 \
  --hist-bins sturges \
  --out-dir swapped_tsvs_chr_specific
```

**3) 5 random chromosomes; custom seeds; explicit chrom sizes**

```
python3 swap_chromosomes.py \
  --n-chroms 5 \
  --pair-seed 123 --chrom-seed 456 \
  --chr-sizes /tbusa/scratch/kdistor/sequencing_Data/simul_trip_downsampled/hg38_chr_sizes_1_22_XY \
  --out-dir swapped_tsvs_5rand
```

> If you `chmod +x` the script, you can run `./swap_chromosomes.py` instead of `python3 swap_chromosomes.py`.

---

## Notes & Tips

- **Performance:** The script reads each TSV once per relevant decision; for very large files, consider running on local SSD or using compressed TSV + preprocessing if I/O becomes the bottleneck.
- **Headers:** Only preserved if the first non‑empty row has `chr` in the 3rd column (case‑insensitive). All row counts in the log exclude headers.
- **Missing chrom sizes:** Any chromosome absent from the sizes file contributes 0 to `swapped_chr_total_bp` (a warning is printed).
- **Binning:** `--hist-bins fd` uses Freedman–Diaconis (robust to outliers). `sturges` is a common alternative for smaller sample sizes. You can also pass an integer count of bins.

---

## Troubleshooting

- ``** for **``**:** The run still proceeds; size totals will be 0 and a warning prints. Provide the correct path or disable size summing by ignoring that column.
- **No PNG produced:** Likely `matplotlib` isn’t installed or there’s no display. The script uses a headless backend; install `matplotlib` or use only the TSV.
- **Empty histogram / one bin:** Can happen with very few runs or identical totals. Try a different binning strategy or run more pairs.
- **Unexpected filenames (tags):** The script extracts tags beginning at the first `1604-` and normalizes underscores to hyphens. Adjust `SAMPLE_TAG_RE` if your naming changes.

---

## Changelog (relevant excerpts)

- **Aug 14, 2025**
  - Added chromosome size summation (`swapped_chr_total_bp`)
  - Added histogram export (`swapped_bp_hist.tsv`, `swapped_bp_hist.png`)
  - X‑axis label set to **Total Trisomy Size (Mb)**; title to **Counts vs Total Trisomy Size (Mb)**; scientific notation
  - Logged row counts: `dip_rows_to_swap`, `trip_rows_to_swap`, `dip_total_rows`, `trip_total_rows`, `swapped_file_rows`
  - Sorting no longer special‑cases `chrX`/`chrY` numerically (only chr1–22 have numeric priority)

---

## License

Internal/research use. Add a license here if distributing.

