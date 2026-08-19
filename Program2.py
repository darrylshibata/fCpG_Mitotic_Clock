import csv
import os
import random


def batch_generate_cellnumber_tables(
    input_dirs=['variants_glandvar_master'],
    split_threshold=100,
    num_splits=18,
    output_dir='cellnumber_tables_master',
    phase5_options=None
):
    """
    Batch generate cellnumber_table_split csv files from qrs_flip variants
    produced by the gland-variance-accepting generator.

    Applies Phase 5 maintenance after the 18 splits, starting from
    the actual post-division population (produced = q + r*2) of the last split.

    For Phase 5A: forces r == s exactly after rounding the 0.05 fraction.
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    if phase5_options is None:
        phase5_options = [
            ('none', 0),
            ('A', 50),
            ('A', 100),
            ('A', 200),
            ('B', 100),
            ('B', 200),
        ]

    positions = {
        'E': (3, 4, 5),
        'M': (8, 9, 10),
        'L': (12, 13, 14)
    }

    for input_dir in input_dirs:
        print(f"Processing directory: {input_dir}")
        for file in os.listdir(input_dir):
            if file.endswith('.csv') and file.startswith('var'):
                qrs_flip_file = os.path.join(input_dir, file)

                # Extract trailing name including varXXXX_
                base_name = file.replace('.csv', '')
                trailing = base_name

                # Generate base table up to 18 splits
                temp_table = 'temp_cellnumber_table.csv'
                generate_cellnumber_table(split_threshold, num_splits, qrs_flip_file, temp_table)

                # Read base table
                with open(temp_table, 'r') as f:
                    reader = csv.DictReader(f)
                    base_rows = list(reader)

                # Find positions of split rows
                split_positions = {}
                for i, row in enumerate(base_rows):
                    if row['split'] != '0':
                        try:
                            split_num = int(row['split'].replace('split', ''))
                            split_positions[split_num] = i
                        except ValueError:
                            pass

                # For each Phase 5 choice
                for p5_mode, p5_divs in phase5_options:
                    # Extend the base table with Phase 5 if needed
                    extended_rows = [row.copy() for row in base_rows]

                    # Get the ACTUAL final population after last split
                    last_row = extended_rows[-1]
                    try:
                        last_q = int(last_row['q'])
                        last_r = int(last_row['r'])
                        current_mothers = last_q + last_r * 2
                    except (KeyError, ValueError):
                        # Fallback (should not happen in normal runs)
                        current_mothers = int(last_row['desired_population'])

                    if p5_divs > 0:
                        current_div = int(extended_rows[-1]['division']) + 1
                        for maint_i in range(1, p5_divs + 1):
                            desired = current_mothers
                            split_label = f'maint{maint_i}'

                            if p5_mode == 'A':
                                # Compute ONE rounded value for the 0.05 fraction
                                turnover_count = round(0.05 * desired)
                                # Force r == s
                                r_count = turnover_count
                                s_count = turnover_count
                                # q gets the remainder → exact total
                                q_count = desired - r_count - s_count
                            else:  # 'B'
                                q_count = desired
                                r_count = 0
                                s_count = 0

                            produced = q_count * 1 + r_count * 2
                            extended_rows.append({
                                'division': str(current_div),
                                'desired_population': str(desired),
                                'q': str(q_count),
                                'r': str(r_count),
                                's': str(s_count),
                                'split': split_label
                            })
                            current_mothers = produced
                            current_div += 1

                    # Now apply branch labeling to the (possibly extended) table
                    for pos, splits in positions.items():
                        new_rows = [row.copy() for row in extended_rows]
                        for row in new_rows:
                            row['split'] = '0'

                        # AB at split 1
                        if 1 in split_positions:
                            new_rows[split_positions[1]]['split'] = 'AB'

                        # Label 12,13,14
                        labels = ['12', '13', '14']
                        for j, s in enumerate(splits):
                            if s in split_positions:
                                new_rows[split_positions[s]]['split'] = labels[j]

                        # Build output filename with Phase 5 annotation
                        p5_tag = f"_p5{p5_mode}{p5_divs}" if p5_mode != 'none' else '_p5none'
                        output_file = os.path.join(
                            output_dir,
                            f'cellnumber_split_{trailing}{p5_tag}_{pos}.csv'
                        )

                        with open(output_file, 'w', newline='') as csvfile:
                            writer = csv.DictWriter(csvfile, fieldnames=new_rows[0].keys())
                            writer.writeheader()
                            writer.writerows(new_rows)

                        print(f"  Generated: {os.path.basename(output_file)}")

                # Clean up temp file
                if os.path.exists(temp_table):
                    os.remove(temp_table)


def generate_cellnumber_table(split_threshold, num_splits, qrs_flip_file='qrs_flip_random.csv', output_file=None):
    """
    Generate cellnumber table from qrs_flip file (unchanged core logic).
    """
    qrs_flip = {}
    with open(qrs_flip_file, 'r') as file:
        reader = csv.DictReader(file)
        reader.fieldnames = [name.strip() for name in reader.fieldnames if name is not None]
        for row in reader:
            split_num = int(row['split'].strip())
            q = float(row['q'])
            r = float(row['r'])
            s = float(row['s'])
            grow = int(row['growN'])
            qrs_flip[split_num] = (q, r, s, grow)

    table = []
    current_div = 0
    current_mothers = 1
    splits_count = 0
    current_phase = 0

    # Div 0
    desired = current_mothers
    q_p, r_p, s_p, grow_p = qrs_flip[0]
    split_label = '0'

    if abs(q_p + r_p + s_p - 1.0) < 1e-6:
        q_count = round(q_p * desired)
        r_count = round(r_p * desired)
        s_count = desired - q_count - r_count
        q_count = max(0, q_count)
        r_count = max(0, r_count)
        s_count = max(0, s_count)
        if q_count + r_count + s_count != desired:
            s_count = desired - q_count - r_count
    else:
        effective_grow = max(min(grow_p, desired), -desired)
        min_s = max(0, -effective_grow)
        max_s = (desired - effective_grow) // 2
        if min_s > max_s:
            raise ValueError(f"Cannot achieve growth {grow_p} for desired {desired}")
        s_count = random.randint(min_s, max_s)
        r_count = s_count + effective_grow
        q_count = desired - r_count - s_count

    produced = q_count * 1 + r_count * 2
    table.append((current_div, desired, q_count, r_count, s_count, split_label))
    current_mothers = produced
    current_div += 1

    while splits_count < num_splits:
        desired = current_mothers
        split_label = '0'

        if desired > split_threshold:
            q_count = desired // 2
            r_count = 0
            s_count = desired - q_count
            split_label = f'split{splits_count + 1}'
            splits_count += 1
            current_phase = splits_count
        else:
            if current_phase not in qrs_flip:
                raise ValueError(f"No qrs for phase {current_phase}")
            q_p, r_p, s_p, grow_p = qrs_flip[current_phase]

            if abs(q_p + r_p + s_p - 1.0) < 1e-6:
                q_count = round(q_p * desired)
                r_count = round(r_p * desired)
                s_count = desired - q_count - r_count
                q_count = max(0, q_count)
                r_count = max(0, r_count)
                s_count = max(0, s_count)
                if q_count + r_count + s_count != desired:
                    s_count = desired - q_count - r_count
                if s_count < 0:
                    adjust = -s_count
                    if q_count >= adjust:
                        q_count -= adjust
                    else:
                        r_count -= (adjust - q_count)
                        q_count = 0
                    s_count = 0
                if q_count < 0 or r_count < 0:
                    s_count = desired - max(0, q_count) - max(0, r_count)
            else:
                effective_grow = max(min(grow_p, desired), -desired)
                min_s = max(0, -effective_grow)
                max_s = (desired - effective_grow) // 2
                if min_s > max_s:
                    raise ValueError(f"Cannot achieve growth {grow_p} for desired {desired}")
                s_count = random.randint(min_s, max_s)
                r_count = s_count + effective_grow
                q_count = desired - r_count - s_count

        produced = q_count * 1 + r_count * 2
        table.append((current_div, desired, q_count, r_count, s_count, split_label))
        current_mothers = produced
        current_div += 1

    if output_file is None:
        output_file = f'cellnumber_table_split{split_threshold}_splits{num_splits}.csv'

    with open(output_file, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['division', 'desired_population', 'q', 'r', 's', 'split'])
        for row in table:
            writer.writerow(row)

    print(f"Generated {output_file} with {len(table)} rows ({len(table)-1} divisions) and {num_splits} splits.")


if __name__ == "__main__":
    # ── Screening mode: pick which Phase 5 options to generate ──
    use_none   = False    # no extra divisions
    use_A10    = True
    use_A50    = True
    use_A100   = True
    use_A200   = True
    use_B10    = False   # ← example: turn off B10 to reduce variants
    use_B50    = False
    use_B100   = False
    use_B200   = False

    # Build the list dynamically from your choices
    phase5_options = []
    if use_none:
        phase5_options.append(('none', 0))
    if use_A10:
        phase5_options.append(('A', 10))
    if use_A50:
        phase5_options.append(('A', 50))
    if use_A100:
        phase5_options.append(('A', 100))
    if use_A200:
        phase5_options.append(('A', 200))
    if use_B10:
        phase5_options.append(('B', 10))
    if use_B50:
        phase5_options.append(('B', 50))
    if use_B100:
        phase5_options.append(('B', 100))
    if use_B200:
        phase5_options.append(('B', 200))

    # If nothing selected, at least include 'none' as fallback
    if not phase5_options:
        phase5_options = [('none', 0)]

    batch_generate_cellnumber_tables(
        input_dirs=['variants_glandvar_master'],
        output_dir='cellnumber_tables_master_phase5',
        phase5_options=phase5_options
    )
