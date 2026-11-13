# fCpG_Mitotic_Clock

Code for **“Early Colorectal Cancer Growth Is Often Slow and Prolonged.”**

## Overview

This repository contains Python programs used to generate candidate colorectal cancer (CRC) growth histories and simulate rapidly fluctuating CpG (fCpG) methylation patterns for approximate Bayesian computation (ABC).

The programs form a sequential workflow. Program1 generates growth-parameter files. Program2 uses those files to construct cell-number histories with spatial gland splitting and terminal Phase 5 options. Program3 uses the cell-number tables to simulate methylation and calculate summary statistics.

**The original programs used in the study are preserved. Important details about a placeholder in Program1 and output metadata in Program3 are documented below.**

## Files

### Program1.py

Generates candidate growth-parameter CSV files. Division outcomes are represented as `q` (one surviving daughter), `r` (two daughters), and `s` (no surviving daughters). Growth-parameter files have columns `split`, `q`, `r`, `s`, and `growN`.

Current settings:

- Generates 100 accepted candidate files in `variants_glandvar_master/`.
- Uses a gland-splitting threshold of 80 cells and 18 splitting events.
- Phase 1 uses `q=0, r=1, s=0` until the first split.
- The numbers of split intervals assigned to Phases 2–4 are sampled from 0–18, retaining combinations that sum to 18. **These numbers are split intervals, not numbers of cell divisions.**
- Each nonempty phase is assigned subtype A or B with equal probability.
- Subtype A samples net growth (`growN`) from `1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 15, 20, 30, 50, 70, 100`. Counts of division outcomes are then generated subject to the net-growth constraint and current population size.
- Subtype B samples `r` from `0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10, 0.12, 0.15, 0.20, 0.30, 0.50`, with `q=1-r` and `s=0` in the parameter file.
- Each list is sampled using Python's `random.choice`.

**Placeholder note:** `compute_gland_variance_from_table()` does not calculate methylation variance. It returns a random value sampled uniformly from 0.045–0.095 and rounded to six decimal places. Candidates are retained when this value is between 0.048 and 0.080, inclusive. This is preliminary random screening, not methylation-based ABC acceptance. The value is included in filenames as `gvar` and must not be interpreted as a biological measurement. Actual methylation and its variances are calculated in Program3. The placeholder also consumes random numbers and therefore forms part of the original stochastic workflow.

Output: Candidate parameter CSV files beginning with `var`, in `variants_glandvar_master/`. Temporary cell-number tables are generated during screening and deleted afterward.

### Program2.py

Reads Program1's candidate parameter files and generates cell-number histories. Each table contains:

`division, desired_population, q, r, s, split`

The script regenerates the cell-number history from the candidate parameters; **it does not reuse the temporary history generated during Program1 screening.** Consequently, stochastic outcome counts and total division numbers can differ between the two realizations.

Current settings:

- Input directory: `variants_glandvar_master/`.
- Output directory: `cellnumber_tables_master_phase5/`.
- Splitting threshold: 80 cells; number of splitting events: 18.
- The first split is labeled `AB`, separating the two tumor sides.
- Three within-side branching configurations are generated: early (`E`: split events 3, 4, 5), middle (`M`: 8, 9, 10), and late (`L`: 12, 13, 14).
- Terminal Phase 5A durations of 10, 50, 100, and 200 additional steps are enabled.
- In Phase 5A, `r=s=round(0.05 × population)` and `q=population-r-s`, maintaining the population size with turnover.
- The no-Phase-5 option and Phase 5B options are disabled in the current main block. If enabled, Phase 5B uses `q=population`, `r=s=0`.

With 100 input candidates, these settings generate **1,200 cell-number tables**: 100 candidates × 4 Phase 5 durations × 3 branching configurations.

Output filenames begin with `cellnumber_split_` and include the candidate filename, Phase 5 setting, and `E`, `M`, or `L` suffix.

### Program3.py

Reads the cell-number tables and `lookup5_table.csv`, simulates methylation on two alleles per cell, and reconstructs eight sampled glands: `A1–A4` and `B1–B4`.

Current settings:

- Input directory: `cellnumber_tables_master_phase5/`.
- Lookup table: `lookup5_table.csv` in the working directory.
- Number of simulation runs per input table: `num_runs=1`.
- All `.csv` files in the input directory are processed in sorted filename order.

Outputs, written in the working directory:

- `sumlineage_arrays.csv`: per-site mean methylation values and fCpG variance for each gland.
- `summary_stats.csv`: candidate metadata, eight gland variances, variance of the across-gland mean methylation profile, pairwise Pearson correlations, and mean variance/correlation summaries for the two sides and combined glands. It also includes early-lineage barcode counts and maximum barcode proportions.

The seven aggregate methylation summaries are `Ave_Gland_Var_A`, `Ave_Gland_Var_B`, `Ave_Gland_Var_All`, `Overall_Variance`, `Ave_Pearson_A`, `Ave_Pearson_B`, and `Ave_Pearson_Cross`.

**Output metadata note:** The original filename parser reads Phase B rates as integers. A filename rate such as `0.02` or `0.50` is therefore reported as `0` in the corresponding `Phase2B_r`, `Phase3B_r`, or `Phase4B_r` column. Do not use these parsed columns to recover the rates; use the decimal values retained in the full filename in the `Var` column or the original Program1 parameter files. Simulation outcome counts are read from the cell-number table, not from these parsed rate columns.

**The `Divisions` and `GVar` metadata fields are also parsed from the original candidate filename.** `GVar` is the Program1 placeholder value. `Divisions` describes the temporary Program1 realization, not necessarily the regenerated Program2 history or its added Phase 5 steps. Read the cell-number table to determine its actual division indices.

**Barcode implementation:** The barcode columns use the first four entries of each lineage array and count the 16 possible sequences of `1` and `2`. Columns ending in `maxprop_4div` report the maximum such barcode count divided by the number of cells in that gland. These are implementation details of this program; the script writes the quantities but does not apply a clonality-acceptance threshold.

**Both output CSVs are opened in append mode.** Use a fresh working directory for each independent batch, or move previous outputs elsewhere before rerunning. Otherwise, results and headers from separate runs can accumulate in the same files. Avoid running multiple instances in the same working directory: the programs also use fixed temporary filenames.

### lookup5_table.csv

Required by Program3. The first row is a header. Each subsequent row represents a simulated CpG site; columns 2 and 3 contain the starting binary methylation states of the two alleles, and columns 4 and 5 contain the `0→1` and `1→0` transition probabilities. Column 1 is not used by the reader.

The number of simulated sites is determined by the number of data rows in this file, not by a fixed site count in Program3. Use the study's original lookup table; it is not generated by Programs1–3.

## 1. System requirements

- Python 3.
- NumPy.
- SciPy (`scipy.stats.pearsonr`).
- Standard Python libraries used include `csv`, `os`, `random`, `itertools`, and `re`.

### Tested software and operating system

The small demonstration was successfully executed on 13 November 2026 with:

| Component | Tested version |
|---|---|
| Operating system | Debian GNU/Linux 13 (trixie), x86-64; Linux kernel 6.18.57, glibc 2.41 |
| Python | 3.12.14 |
| NumPy | 2.3.5 |
| SciPy | 1.17.0 |

These are versions verified for the newly supplied demonstration, not a claim about the original manuscript-analysis environment. The exact original versions are not recorded in the scripts. Windows and macOS have not been tested for this package.

### Hardware

No GPU, accelerator, cluster, or other non-standard hardware is required by the code. A standard CPU computer capable of running Python, NumPy, and SciPy is sufficient for the small demo. Large simulation collections can require substantially more time and storage.

## 2. Installation guide

1. Download the code package, extract it, and open a terminal in the extracted directory.
2. Use Python 3.12 to match the tested Python release family.
3. Install the tested dependency versions:

```bash
python -m pip install -r requirements-demo.txt
```

`requirements-demo.txt` specifies NumPy 2.3.5 and SciPy 1.17.0. All other imported modules are part of Python's standard library.

**Typical installation time:** allow approximately 1–5 minutes on a desktop with Python already installed, a working internet connection, and prebuilt package wheels. This is a planning estimate, not a measured fresh-install benchmark; the verification environment already had the dependencies installed. Installing Python itself or building packages from source may take longer.

To record the environment used on another computer:

```bash
python --version
python -m pip show numpy scipy
```

## 3. Demo

The package includes a small input dataset, `lookup5_table.csv`, containing 100 simulated CpG sites, and a separate `demo.py` runner. The runner imports the original Programs1–3 without changing their source files. It requests one candidate, one Phase 5A duration (10 steps), and the three E/M/L branching configurations. It sets both random-number seeds to 1729.

### Instructions

With `Program1.py`, `Program2.py`, `Program3.py`, `lookup5_table.csv`, and `demo.py` in the same extracted directory, run:

```bash
python demo.py
```

The runner creates a new `demo_output/` directory. If that directory already exists, the runner stops rather than append to or overwrite prior results. Move or rename the previous directory before rerunning.

### Expected output

| Path under demo_output/ | Expected contents |
|---|---|
| variants_glandvar_master/ | 1 candidate parameter CSV |
| cellnumber_tables_master_phase5/ | 3 cell-number CSVs, one per E/M/L configuration, each with 10 Phase 5A steps |
| summary_stats.csv | Header plus 3 simulation-summary rows |
| sumlineage_arrays.csv | Header plus 24 gland rows: 8 glands × 3 simulations |
| lookup5_table.csv | Copy of the demonstration input |
| demo_run_info.json | Seeds, software environment, elapsed times, and output row counts |

Example outputs from the verified run are supplied in `expected_outputs/`. All seven aggregate methylation summaries were checked for finite values. Floating-point values and stochastic trajectories may differ across software versions; the expected file structure and row counts provide basic execution checks.

### Expected runtime

The demonstration's simulation stages completed in **6.3 seconds** in the verification environment on an Intel Xeon Platinum 8573C CPU, with numerical-library thread counts limited to one. This excludes dependency installation and initial Python/package imports. Program3 accounted for approximately 6.29 seconds; Programs1 and 2 each took less than 0.01 seconds.

This was a CPU-only server/container measurement, not a desktop benchmark. Allow a few minutes on a typical desktop; actual time depends on hardware and software. The measured timing and environment are recorded in `expected_outputs/demo_run_info.json`. The full default workflow processes 1,200 tables and will take much longer than this demonstration.

## 4. Instructions for use

### Full default simulation workflow

Place `Program1.py`, `Program2.py`, `Program3.py`, and the intended `lookup5_table.csv` in a fresh working directory. Run:

```bash
python Program1.py
python Program2.py
python Program3.py
```

The scripts have no command-line argument interface. Their settings are defined in the functions and the `if __name__ == "__main__":` blocks. Paths are relative to the directory from which Python is launched, not automatically to the directory containing the scripts.

Program1 creates `variants_glandvar_master/`. Program2 reads it and creates `cellnumber_tables_master_phase5/`. Program3 reads those tables and writes its two output CSVs in the working directory.

### Running with different simulation inputs

Use a separate working copy when changing settings, preserving the original scripts. Candidate numbers and growth settings are selected in Program1. Branching and Phase 5 options are selected in Program2. The number of repeated methylation simulations and input directory are set in Program3.

For a different set of simulated CpGs, provide a lookup CSV in the documented five-column format, with a header, binary initial states in columns 2–3, and transition probabilities in [0,1] in columns 4–5. Program3 reads the filename `lookup5_table.csv` from the working directory.

Experimental gland methylation data are not a direct input to Programs1–3. These programs generate simulations for subsequent comparison with experimental summaries, as explained below.

## ABC comparison and acceptance

These three scripts generate simulations and summary statistics. **They do not load the experimental CRC methylation summaries or implement the manuscript's final ABC acceptance step.**

In the manuscript analysis, the seven designated methylation summary statistics were each required to fall within ±10% of their corresponding experimental values. Any additional gland-clonality filter and the downstream parameter summaries must be applied separately. The Program1 placeholder screening is not this acceptance step.

## Reproducibility

The programs are stochastic and do not set random-number seeds internally. Program1 and Program2 use Python's `random` module; Program3 uses both `random` and NumPy's random-number generator. Reproducing a particular run requires recording the seeds, software versions, input files, parameters, and input-file processing order. Seeding only one random-number generator is insufficient for the full workflow.

Program2 processes candidate files using `os.listdir` without sorting, so file order can affect the assignment of random draws to candidates. Program3 sorts its input CSV filenames.

The manuscript analysis used a large collection of candidate histories generated across simulation batches. The current default of 100 candidates is a batch setting; running these defaults once does not regenerate the manuscript's entire simulation collection.

## Data

Experimental methylation and sequencing data should be obtained from the sources described in the manuscript and its Data Availability statement. Controlled-access or patient-level data are not included in this repository.

Previously published single-cell colorectal tumor methylation data are available under GEO accession GSE97693. The accession for the tumor-gland EPIC methylation data will be added when available.

## Citation

Ryser MD and Shibata D. *Early Colorectal Cancer Growth Is Often Slow and Prolonged.* Manuscript submitted for publication. The final publication citation and DOI will be added when available.

## License

A software license has not yet been specified.
