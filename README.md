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

## Requirements

- Python 3.
- NumPy.
- SciPy (`scipy.stats.pearsonr`).
- Standard Python libraries used include `csv`, `os`, `random`, `itertools`, and `re`.

Install the external dependencies with:

```bash
python -m pip install numpy scipy
```

Exact Python and package versions used for the original analysis are not specified in these scripts. Record the environment used for any new run:

```bash
python --version
python -m pip show numpy scipy
```

## Running the workflow

Place `Program1.py`, `Program2.py`, `Program3.py`, and `lookup5_table.csv` in one working directory. The commands below assume the scripts have those filenames. Run them sequentially from that directory:

```bash
python Program1.py
python Program2.py
python Program3.py
```

The scripts have no command-line argument interface. Their settings are defined in the functions and the `if __name__ == "__main__":` blocks. Paths are relative to the directory from which Python is launched, not automatically to the directory containing the scripts.

Program1 creates `variants_glandvar_master/`. Program2 reads it and creates `cellnumber_tables_master_phase5/`. Program3 reads those tables and writes its two output CSVs in the working directory.

### Smaller demonstration run

To demonstrate the workflow without generating the default 1,200 tables, make a separate working copy of the scripts and lookup table, preserving the original programs.

1. In the working copy of Program1, change the final main-block call from `num_variants=100` to `num_variants=1`.
2. In the working copy of Program2, leave `use_A10=True` and set `use_A50`, `use_A100`, and `use_A200` to `False`. Keep the other Phase 5 switches `False`.
3. Run the three commands above in the new working directory.

This configuration requests one candidate and generates three cell-number tables, one for each branching configuration, each with 10 terminal Phase 5A steps. Program3 attempts one methylation simulation per table. Summary rows are written only when all eight gland profiles are available. Numerical results vary between runs.

This is a reduced configuration derived from the script settings, **not a separately benchmarked example**. Runtime and peak memory have not been measured here; both depend on the sampled history, lookup-table size, software environment, and hardware. A full end-to-end execution check requires the original `lookup5_table.csv`.

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

A license has not yet been specified. A LICENSE file and the corresponding license designation should be added once selected by the authors.
