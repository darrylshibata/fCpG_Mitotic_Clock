import csv
import random
import numpy as np
from scipy.stats import pearsonr
from itertools import combinations, product
import os
import re

def parse_filename_detailed(filename):
    """
    Parses complex metadata from the filename.
    Handles p2A/p2B logic, divisions, gvar, and split character.
    """
    params = {
        'Var': filename, # Per your request: full title of the CSV
        'Phase2A_rows': 0, 'Phase2A_growN': 0, 'Phase2B_rows': 0, 'Phase2B_r': 0,
        'Phase3A_rows': 0, 'Phase3A_growN': 0, 'Phase3B_rows': 0, 'Phase3B_r': 0,
        'Phase4A_rows': 0, 'Phase4A_growN': 0, 'Phase4B_rows': 0, 'Phase4B_r': 0,
        'Divisions': 0, 'GVar': 0.0, 'Split': ''
    }

    # 1. Parse Phases (p2, p3, p4)
    for p_num in ['2', '3', '4']:
        a_match = re.search(f'p{p_num}A(\d+)_(\d+)', filename)
        if a_match:
            params[f'Phase{p_num}A_rows'] = int(a_match.group(1))
            params[f'Phase{p_num}A_growN'] = int(a_match.group(2))
            continue
        
        b_match = re.search(f'p{p_num}B(\d+)_(\d+)', filename)
        if b_match:
            params[f'Phase{p_num}B_rows'] = int(b_match.group(1))
            params[f'Phase{p_num}B_r'] = int(b_match.group(2))
            continue

    div_match = re.search(r'div(\d+)', filename)
    if div_match:
        params['Divisions'] = int(div_match.group(1))

    gvar_match = re.search(r'gvar([\d.]+)', filename)
    if gvar_match:
        params['GVar'] = float(gvar_match.group(1))

    split_match = re.search(r'_([EML])\.csv$', filename)
    if split_match:
        params['Split'] = split_match.group(1)

    return params

def read_cellnumber_table(cellnumber_table_file):
    cellnumber_table = []
    with open(cellnumber_table_file, 'r') as file:
        reader = csv.DictReader(file)
        for row in reader:
            cellnumber_table.append((
                int(row['division']), int(row['desired_population']),
                int(row['q']), int(row['r']), int(row['s']), row['split']
            ))
    return cellnumber_table

def read_lookup_table(lookup_table_file):
    array1, array2, p0_1, p1_0 = [], [], [], []
    with open(lookup_table_file, 'r') as file:
        reader = csv.reader(file)
        next(reader)
        for row in reader:
            array1.append(int(row[1])); array2.append(int(row[2]))
            p0_1.append(float(row[3])); p1_0.append(float(row[4]))
    return (array1, array2, p0_1, p1_0)

def shuffle_qrs_pool(q, r, s):
    qrs_pool = ['q'] * q + ['r'] * r + ['s'] * s
    random.shuffle(qrs_pool)
    return qrs_pool

def cell_division(mother_cells, qrs_pool, cellnumber_table, division_index, saved_discards, current_split):
    daughter_cells = []
    division_number = cellnumber_table[division_index][0]
    for mother_cell in mother_cells:
        lineage_array, array1, array2, p0_1, p1_0 = mother_cell
        marker = qrs_pool.pop(0)
        if marker in ['q', 'r']:
            branches = [1] if marker == 'q' else [1, 2]
            for branch in branches:
                daughter_lineage = np.copy(lineage_array)
                daughter_lineage[division_number] = branch
                d_a1, d_a2 = np.copy(array1), np.copy(array2)

                if current_split == '0':
                    for arr in [d_a1, d_a2]:
                        rand = np.random.rand(len(arr))
                        arr[(arr == 0) & (rand < p0_1)] = 2
                        arr[(arr == 1) & (rand < p1_0)] = 0
                        arr[arr == 2] = 1
                daughter_cells.append((daughter_lineage, d_a1, d_a2, p0_1, p1_0))
        elif marker == 's' and current_split != '0':
            if current_split not in saved_discards: saved_discards[current_split] = []
            saved_discards[current_split].append(mother_cell)
    return daughter_cells

def calculate_and_write_average_bit_values(cells, run, gland, gland_averages, gland_variances, csv_filename):
    if not cells: return
    num_bits, num_cells = len(cells[0][1]), len(cells)
    combined_sums = np.zeros(num_bits)
    for cell in cells: combined_sums += (cell[1] + cell[2])
    averages = combined_sums / (num_cells * 2)
    variance = np.var(averages)
    gland_averages[gland], gland_variances[gland] = averages, variance

    with open('sumlineage_arrays.csv', 'a', newline='') as csvfile:
        writer = csv.writer(csvfile)
        if run == 1 and gland == 'A1' and csv_filename == sorted(os.listdir('cellnumber_tables_master_phase5'))[0]:
            writer.writerow(['CSV_File', 'Run_Gland'] + [f'bit{i+1}' for i in range(num_bits)] + ['Variance'])
        writer.writerow([csv_filename, f'Run{run}_{gland}'] + list(averages) + [variance])

def simulate_from(start_index, initial_mothers, cellnumber_table, run, gland, last_run_lineage_arrays, is_last_run, gland_averages, gland_variances, csv_filename):
    mother_cells = initial_mothers[:]
    if len(mother_cells) > cellnumber_table[start_index][1]:
        random.shuffle(mother_cells)
        mother_cells = mother_cells[:cellnumber_table[start_index][1]]
    saved_discards = {}
    for i in range(start_index, len(cellnumber_table)):
        div, desired, q, r, s, split = cellnumber_table[i]
        total = q + r + s
        if total != len(mother_cells):
            factor = len(mother_cells) / total if total > 0 else 0
            q, r = round(q * factor), round(r * factor)
            s = len(mother_cells) - q - r
        mother_cells = cell_division(mother_cells, shuffle_qrs_pool(q, r, s), cellnumber_table, i, saved_discards, split)
    calculate_and_write_average_bit_values(mother_cells, run, gland, gland_averages, gland_variances, csv_filename)
    if is_last_run: last_run_lineage_arrays[gland] = mother_cells
    return saved_discards

def write_summary_stats(run, gland_averages, gland_variances, csv_filename, first_file, params, last_run_lineage_arrays):
    a_glands = ['A1', 'A2', 'A3', 'A4']
    b_glands = ['B1', 'B2', 'B3', 'B4']
    all_glands = a_glands + b_glands
    if len(gland_averages) != 8:
        return

    num_bits = len(next(iter(gland_averages.values())))
    per_bit_avgs = [np.mean([gland_averages[g][i] for g in all_glands]) for i in range(num_bits)]
    overall_var = np.var(per_bit_avgs)

    p_a = [pearsonr(gland_averages[g1], gland_averages[g2])[0] for g1, g2 in combinations(a_glands, 2)]
    p_b = [pearsonr(gland_averages[g1], gland_averages[g2])[0] for g1, g2 in combinations(b_glands, 2)]
    p_cross = [pearsonr(gland_averages[ga], gland_averages[gb])[0] for ga, gb in product(a_glands, b_glands)]

    meta_keys = [
        'Var', 'Phase2A_rows', 'Phase2A_growN', 'Phase2B_rows', 'Phase2B_r',
        'Phase3A_rows', 'Phase3A_growN', 'Phase3B_rows', 'Phase3B_r',
        'Phase4A_rows', 'Phase4A_growN', 'Phase4B_rows', 'Phase4B_r',
        'Divisions', 'GVar', 'Split'
    ]
    
    header = meta_keys + [f'Gland_Var_{g}' for g in all_glands] + ['Overall_Variance']
    header += [f'Pearson_A_{p[0]}-{p[1]}' for p in combinations(a_glands, 2)]
    header += [f'Pearson_B_{p[0]}-{p[1]}' for p in combinations(b_glands, 2)]
    header += [f'Pearson_Cross_{p[0]}-{p[1]}' for p in product(a_glands, b_glands)]
    header += ['Ave_Gland_Var_A', 'Ave_Gland_Var_B', 'Ave_Gland_Var_All', 'Ave_Pearson_A', 'Ave_Pearson_B', 'Ave_Pearson_Cross']

    # === NEW: Prepare all 16 possible 4-bit barcodes ===
    possible_barcodes = [''.join(map(str, combo)) for combo in product([1, 2], repeat=4)]
    # → ['1111', '1112', ..., '2222']

    # === NEW: Count barcodes and compute max proportion per gland ===
    barcode_counts = {}           # gland → {barcode: count}
    max_proportions = {}          # gland → float

    for gland in all_glands:
        cells = last_run_lineage_arrays.get(gland, [])
        count_dict = {bc: 0 for bc in possible_barcodes}
        total_cells = len(cells)

        for cell in cells:
            lineage = cell[0]
            if len(lineage) >= 4:
                barcode = ''.join(map(str, lineage[:4]))
                if barcode in count_dict:
                    count_dict[barcode] += 1

        barcode_counts[gland] = count_dict
        max_proportions[gland] = max(count_dict.values()) / total_cells if total_cells > 0 else 0.0

    # === Build header with new columns ===
    for gland in all_glands:
        for bc in sorted(possible_barcodes):
            header.append(f'gld{gland}_{bc}')

    for gland in all_glands:
        header.append(f'gld{gland}_maxprop_4div')

    # === Build row ===
    row = [params[k] for k in meta_keys] + \
          [gland_variances.get(g, 0) for g in all_glands] + [overall_var] + \
          p_a + p_b + p_cross + \
          [np.mean([gland_variances.get(g, 0) for g in a_glands]),
           np.mean([gland_variances.get(g, 0) for g in b_glands]),
           np.mean(list(gland_variances.values())),
           np.mean(p_a), np.mean(p_b), np.mean(p_cross)]

    # === Add barcode counts ===
    for gland in all_glands:
        counts = barcode_counts[gland]
        for bc in sorted(possible_barcodes):
            row.append(counts[bc])

    # === Add max proportions ===
    for gland in all_glands:
        row.append(max_proportions[gland])

    # === Write ===
    with open('summary_stats.csv', 'a', newline='') as csvfile:
        writer = csv.writer(csvfile)
        if run == 1 and first_file:
            writer.writerow(header)
        writer.writerow(row)

def main():
    num_runs, input_dir = 1, 'cellnumber_tables_master_phase5'
    all_files = sorted([f for f in os.listdir(input_dir) if f.endswith('.csv')])
    total_files = len(all_files)
    
    print(f"--- Simulation Started ---")
    print(f"Tracking {total_files} files...")
    array1, array2, p0_1, p1_0 = read_lookup_table('lookup5_table.csv')
    
    for idx, csv_filename in enumerate(all_files, 1):
        params = parse_filename_detailed(csv_filename)
        print(f"[{idx}/{total_files}] Processing: {csv_filename}...")
        
        cellnumber_table = read_cellnumber_table(os.path.join(input_dir, csv_filename))
        max_div = max(entry[0] for entry in cellnumber_table)
        init_mother = [(np.zeros(max_div + 1, dtype=int), np.array(array1), np.array(array2), np.array(p0_1), np.array(p1_0))]
        split_map = {entry[5]: i for i, entry in enumerate(cellnumber_table) if entry[5] != '0'}

        for run in range(1, num_runs + 1):
            g_avgs, g_vars, last_lineage = {}, {}, {}
            a_saved = simulate_from(0, init_mother, cellnumber_table, run, 'A1', last_lineage, (run==num_runs), g_avgs, g_vars, csv_filename)
            for br in ['12', '13', '14']:
                if br in a_saved:
                    simulate_from(split_map[br]+1, a_saved[br], cellnumber_table, run, f'A{br[1]}', last_lineage, (run==num_runs), g_avgs, g_vars, csv_filename)
            if 'AB' in a_saved:
                b_saved = simulate_from(split_map['AB']+1, a_saved['AB'], cellnumber_table, run, 'B1', last_lineage, (run==num_runs), g_avgs, g_vars, csv_filename)
                for br in ['12', '13', '14']:
                    if br in b_saved:
                        simulate_from(split_map[br]+1, b_saved[br], cellnumber_table, run, f'B{br[1]}', last_lineage, (run==num_runs), g_avgs, g_vars, csv_filename)

            write_summary_stats(run, g_avgs, g_vars, csv_filename, (idx == 1), params, last_lineage)
        
        print(f"    Done.")
            
    print(f"--- All {total_files} files processed successfully. ---")

if __name__ == "__main__":
    main()
