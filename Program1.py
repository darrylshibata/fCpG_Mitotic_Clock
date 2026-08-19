import csv
import random
import os

def compute_gland_variance_from_table(table_csv_path):
    """
    Placeholder: Compute the variance of average methylation across 8 glands.
    
    YOU NEED TO IMPLEMENT THIS FUNCTION based on your methylation model.
    Currently returns dummy value in range that roughly matches target.
    """
    # ------------------- REPLACE THIS WITH YOUR ACTUAL IMPLEMENTATION -------------------
    import numpy as np
    dummy_variance = round(random.uniform(0.045, 0.095), 6)  # biased toward acceptance range for testing
    print(f"  [placeholder] gland variance = {dummy_variance:.6f}")
    return dummy_variance
    # -------------------------------------------------------------------------------------


def generate_variant_qrs_flip(
    min_gland_var=0.048,
    max_gland_var=0.080,
    num_variants=100,
    split_threshold=100,
    num_splits=18,
    output_dir='variants_glandvar_master',
    max_attempts=300000
):
    """
    Generate qrs_flip variants using flexible phase lengths (0–18) and A/B subtypes.
    Acceptance based only on gland variance [min_gland_var, max_gland_var].
    
    Phases:
    - Phase 1: always fixed (split 0: q=0, r=1, s=0, growN=0)
    - Phase 2,3,4: length 0–18 each, sum lengths = 18
    - If length > 0: subtype A (fixed growN + possible loss) or B (low r, proportional)
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    variant_count = 0
    attempt = 0

    GROWN_OPTIONS = [1,2,3,4,5,6,7,8,9,10,12,15,20,30,50,70,100]
    R_OPTIONS = [0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.12,0.15,0.2,0.3,0.5]

    while variant_count < num_variants:
        if attempt > max_attempts:
            print(f"\nMax attempts ({max_attempts}) exceeded. Stopping with {variant_count} valid variants generated.")
            break

        attempt += 1

        # Sample lengths: 0 to 18 inclusive
        phase2_rows = random.randint(0, 18)
        phase3_rows = random.randint(0, 18)
        phase4_rows = random.randint(0, 18)

        if phase2_rows + phase3_rows + phase4_rows != num_splits:
            continue

        # Decide subtypes only for phases that exist
        phase2_type = random.choice(['A', 'B']) if phase2_rows > 0 else '0'
        phase3_type = random.choice(['A', 'B']) if phase3_rows > 0 else '0'
        phase4_type = random.choice(['A', 'B']) if phase4_rows > 0 else '0'

        # Sample parameters depending on type
        phase2_growN = random.choice(GROWN_OPTIONS) if phase2_type == 'A' and phase2_rows > 0 else 0
        phase2_r     = random.choice(R_OPTIONS) if phase2_type == 'B' and phase2_rows > 0 else 0.0

        phase3_growN = random.choice(GROWN_OPTIONS) if phase3_type == 'A' and phase3_rows > 0 else 0
        phase3_r     = random.choice(R_OPTIONS) if phase3_type == 'B' and phase3_rows > 0 else 0.0

        phase4_growN = random.choice(GROWN_OPTIONS) if phase4_type == 'A' and phase4_rows > 0 else 0
        phase4_r     = random.choice(R_OPTIONS) if phase4_type == 'B' and phase4_rows > 0 else 0.0

        # Write temporary qrs_flip.csv
        temp_file = 'temp_qrs_flip.csv'
        with open(temp_file, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(['split ', 'q', 'r', 's', 'growN'])

            split_num = 0
            writer.writerow([split_num, 0, 1, 0, 0])  # Phase 1 fixed
            split_num += 1

            # Phase 2
            for _ in range(phase2_rows):
                if phase2_type == 'A':
                    writer.writerow([split_num, 0, 0, 0, phase2_growN])
                else:  # B
                    writer.writerow([split_num, 1.0 - phase2_r, phase2_r, 0, 0])
                split_num += 1

            # Phase 3
            for _ in range(phase3_rows):
                if phase3_type == 'A':
                    writer.writerow([split_num, 0, 0, 0, phase3_growN])
                else:  # B
                    writer.writerow([split_num, 1.0 - phase3_r, phase3_r, 0, 0])
                split_num += 1

            # Phase 4
            for _ in range(phase4_rows):
                if phase4_type == 'A':
                    writer.writerow([split_num, 0, 0, 0, phase4_growN])
                else:  # B
                    writer.writerow([split_num, 1.0 - phase4_r, phase4_r, 0, 0])
                split_num += 1

        # Simulate cell numbers
        temp_table = 'temp_cellnumber_table.csv'
        generate_cellnumber_table(split_threshold, num_splits, temp_file, temp_table)

        # Compute gland variance
        gland_var = compute_gland_variance_from_table(temp_table)

        # Get division count (for logging / filename only)
        with open(temp_table, 'r') as f:
            num_rows = sum(1 for _ in f) - 1
        total_div = num_rows - 1

        accepted = gland_var is not None and min_gland_var <= gland_var <= max_gland_var

        if accepted:
            variant_count += 1

            # Build phase strings with parameter values (no extra leading dot)
            p2_part = f"p2{phase2_type}{phase2_rows}"
            if phase2_rows > 0:
                if phase2_type == 'A':
                    p2_part += f"_{phase2_growN}"
                else:
                    p2_part += f"_{phase2_r:.2f}"   # ← changed: no leading '_.' → just '_0.02'

            p3_part = f"p3{phase3_type}{phase3_rows}"
            if phase3_rows > 0:
                if phase3_type == 'A':
                    p3_part += f"_{phase3_growN}"
                else:
                    p3_part += f"_{phase3_r:.2f}"

            p4_part = f"p4{phase4_type}{phase4_rows}"
            if phase4_rows > 0:
                if phase4_type == 'A':
                    p4_part += f"_{phase4_growN}"
                else:
                    p4_part += f"_{phase4_r:.2f}"

            # absent phases
            if phase2_rows == 0: p2_part = "p2_00"
            if phase3_rows == 0: p3_part = "p3_00"
            if phase4_rows == 0: p4_part = "p4_00"

            # gland variance string (dot + 5 decimals)
            gvar_str = f"{gland_var:.5f}"

            output_file = os.path.join(
                output_dir,
                f"var{variant_count:04d}_{p2_part}_{p3_part}_{p4_part}_"
                f"div{total_div}_gvar{gvar_str}.csv"
            )

            if os.path.exists(output_file):
                os.remove(output_file)
            os.rename(temp_file, output_file)

            print(f"ACCEPTED  {variant_count:4d}   div={total_div:4d}   "
                  f"gvar={gvar_str}   {p2_part} {p3_part} {p4_part}")

        else:
            gvar_display = f"{gland_var:.5f}" if gland_var is not None else "None"
            print(f"rejected  attempt {attempt:6d}   div={total_div:4d}   "
                  f"gvar={gvar_display:<8}")

        # Cleanup temporary files
        for path in [temp_table, temp_file]:
            if os.path.exists(path):
                try:
                    os.remove(path)
                except OSError:
                    pass


def generate_cellnumber_table(split_threshold, num_splits, qrs_flip_file='qrs_flip_random.csv', output_file=None):
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
    generate_variant_qrs_flip(
        min_gland_var=0.048,
        max_gland_var=0.080,
        num_variants=100,
        output_dir='variants_glandvar_master',
        max_attempts=5000000
    )
