"""Run a small demonstration without editing Program1.py, Program2.py or Program3.py.
Usage: python demo.py
Requires these three scripts and lookup5_table.csv beside this file.
"""
from pathlib import Path
import csv
import importlib.util
import json
import os
import platform
import random
import shutil
import time

import numpy as np
import scipy


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main():
    root = Path(__file__).resolve().parent
    required = [root / name for name in ['Program1.py', 'Program2.py', 'Program3.py', 'lookup5_table.csv']]
    for path in required:
        if not path.is_file():
            raise FileNotFoundError(path)
    output = root / 'demo_output'
    output.mkdir(exist_ok=False)  # Protect append-mode results from accidental reruns.
    shutil.copyfile(root / 'lookup5_table.csv', output / 'lookup5_table.csv')
    p1, p2, p3 = [load_module('demo_program' + str(i), required[i-1]) for i in (1, 2, 3)]
    random.seed(1729)
    np.random.seed(1729)
    timings = {}
    start = time.perf_counter()
    old_cwd = Path.cwd()
    try:
        os.chdir(output)
        t = time.perf_counter()
        p1.generate_variant_qrs_flip(num_variants=1, output_dir='variants_glandvar_master', max_attempts=5000000)
        timings['Program1_seconds'] = time.perf_counter() - t
        candidates = list(Path('variants_glandvar_master').glob('var*.csv'))
        if len(candidates) != 1:
            raise RuntimeError('Expected one candidate, found ' + str(len(candidates)))
        t = time.perf_counter()
        p2.batch_generate_cellnumber_tables(input_dirs=['variants_glandvar_master'], output_dir='cellnumber_tables_master_phase5', phase5_options=[('A', 10)])
        timings['Program2_seconds'] = time.perf_counter() - t
        t = time.perf_counter()
        p3.main()
        timings['Program3_seconds'] = time.perf_counter() - t
        rows = {}
        for filename in ['summary_stats.csv', 'sumlineage_arrays.csv']:
            with open(filename, newline='') as f:
                rows[filename] = len(list(csv.DictReader(f)))
        if rows != {'summary_stats.csv': 3, 'sumlineage_arrays.csv': 24}:
            raise RuntimeError('Unexpected output row counts: ' + str(rows))
        info = {'purpose':'New verification demo, not the original manuscript environment',
                'seed':1729, 'python':platform.python_version(), 'numpy':np.__version__,
                'scipy':scipy.__version__, 'platform':platform.platform(),
                'timings_seconds':timings, 'elapsed_seconds':time.perf_counter()-start,
                'row_counts_excluding_header':rows,
                'candidate_files':len(candidates),
                'cellnumber_tables':len(list(Path('cellnumber_tables_master_phase5').glob('*.csv')))}
        Path('demo_run_info.json').write_text(json.dumps(info, indent=2))
        print(json.dumps(info, indent=2))
    finally:
        os.chdir(old_cwd)


if __name__ == '__main__':
    main()
