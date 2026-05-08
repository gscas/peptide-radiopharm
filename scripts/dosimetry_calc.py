#!/usr/bin/env python3
"""
Dosimetry Calculator for Peptide Radiopharmaceuticals.

Converts biodistribution (%ID/g) data to organ absorbed doses using MIRD formalism.

Usage:
  python3 dosimetry_calc.py --input biodist.csv --nuclide Lu-177 --injected_activity 7.4e9
  python3 dosimetry_calc.py --input biodist.csv --nuclide Ac-225 --organ_masses organ_masses.json

Input CSV format:
  Time_h, Organ, PercentID_per_g
  0.5, Kidney, 5.2
  1.0, Kidney, 4.1
  ...

S-values from ICRP 133 / OLINDA for reference adult male.
"""

import argparse
import csv
import json
import math
import sys
from typing import Dict, List, Tuple, Optional

# S-values (Gy/Bq·s) for reference adult male (ICRP 133 / OLINDA)
# Simplified subset — self-dose S(target ← target) for most relevant organs
# For full S-value matrices, see references/s_values.md

# S-values in Gy/Bq·s — converted from OLINDA/EXM mGy/MBq·h
# Conversion: 1 mGy/(MBq·h) = 2.7778e-13 Gy/(Bq·s)
# These are source→target self-dose values for reference adult male
# Key reference: Stabin & Siegel, RADAR; ICRP Publication 133

_S_SCALE = 2.7778e-13  # mGy/(MBq·h) → Gy/(Bq·s)

S_VALUES_SELF_MGy_MBq_h = {
    # Lu-177 self-dose S-values (mGy/MBq·h)
    "Lu-177": {
        "Kidneys": 2.88e-03,
        "Liver": 1.67e-03,
        "Spleen": 3.97e-03,
        "Red_marrow": 1.82e-03,
        "Bone_surfaces": 3.13e-03,
        "Thyroid": 2.83e-03,
        "Adrenals": 4.82e-03,
        "Salivary_glands": 2.92e-03,
        "Whole_body": 4.40e-05,
    },
    # Y-90 self-dose S-values (mGy/MBq·h)
    "Y-90": {
        "Kidneys": 1.18e-02,
        "Liver": 6.52e-03,
        "Spleen": 1.58e-02,
        "Red_marrow": 7.51e-03,
        "Bone_surfaces": 1.27e-02,
        "Thyroid": 1.11e-02,
        "Adrenals": 1.93e-02,
        "Whole_body": 1.75e-04,
    },
    # Ga-68 self-dose S-values (mGy/MBq·h)
    "Ga-68": {
        "Kidneys": 7.68e-03,
        "Liver": 3.99e-03,
        "Spleen": 1.05e-02,
        "Red_marrow": 3.58e-03,
        "Whole_body": 1.08e-04,
    },
    # Ac-225 self-dose S-values (mGy/MBq·h) — includes daughters
    "Ac-225": {
        "Kidneys": 9.35e-02,
        "Liver": 5.14e-02,
        "Spleen": 1.27e-01,
        "Red_marrow": 5.68e-02,
        "Bone_surfaces": 1.01e-01,
        "Whole_body": 1.43e-03,
    },
    # Tc-99m self-dose S-values (mGy/MBq·h)
    "Tc-99m": {
        "Kidneys": 2.55e-03,
        "Liver": 1.42e-03,
        "Spleen": 3.54e-03,
        "Red_marrow": 1.11e-03,
        "Thyroid": 2.42e-03,
        "Whole_body": 3.07e-05,
    },
    # In-111 self-dose S-values (mGy/MBq·h)
    "In-111": {
        "Kidneys": 6.02e-03,
        "Liver": 3.44e-03,
        "Spleen": 8.47e-03,
        "Red_marrow": 3.24e-03,
        "Whole_body": 7.97e-05,
    },
}

# Convert to Gy/(Bq·s)
S_VALUES_SELF = {}
for nuclide, organs in S_VALUES_SELF_MGy_MBq_h.items():
    S_VALUES_SELF[nuclide] = {organ: s * _S_SCALE for organ, s in organs.items()}

# Default organ masses (g) — ICRP 89 reference adult male
DEFAULT_ORGAN_MASSES = {
    "Kidneys": 310,
    "Liver": 1800,
    "Spleen": 150,
    "Red_marrow": 1170,
    "Bone_surfaces": 1000,
    "Thyroid": 20,
    "Adrenals": 14,
    "Salivary_glands": 82,
    "Whole_body": 73000,
}

# Physical half-lives (seconds)
HALF_LIVES = {
    "Lu-177": 6.65 * 86400,
    "Y-90": 2.67 * 86400,
    "Ga-68": 67.7 * 60,
    "Ac-225": 9.92 * 86400,
    "Tc-99m": 6.0 * 3600,
    "In-111": 2.80 * 86400,
    "Cu-67": 2.58 * 86400,
    "Cu-64": 12.7 * 3600,
    "F-18": 109.8 * 60,
    "Sc-47": 3.35 * 86400,
    "Tb-161": 6.90 * 86400,
}


def read_biodistribution(filepath: str) -> Dict[str, List[Tuple[float, float]]]:
    """Read biodistribution CSV into {organ: [(time_h, percentID_per_g), ...]}"""
    data = {}
    with open(filepath, "r", encoding="utf-8-sig") as f:
        reader = csv.reader(f)
        for row in reader:
            if len(row) < 3:
                continue
            try:
                t = float(row[0].strip())
                organ = row[1].strip()
                pidg = float(row[2].strip())
                if organ not in data:
                    data[organ] = []
                data[organ].append((t, pidg))
            except ValueError:
                continue
    # Sort by time
    for organ in data:
        data[organ].sort(key=lambda x: x[0])
    return data


def compute_cumulated_activity_per_organ(
    time_data: List[Tuple[float, float]],
    half_life_phys: float,
    organ_mass_g: float,
    injected_activity_Bq: float,
) -> float:
    """Calculate cumulated activity Ã (Bq·s) in an organ from %ID/g time course.
    
    Steps:
    1. Convert %ID/g(t) → activity in organ A(t) = (%ID/g / 100) × mass × A0
    2. Integrate A(t) using trapezoidal rule
    3. Add tail (physical decay beyond last time point)
    
    Returns cumulated activity in Bq·s.
    """
    if not time_data:
        return 0.0
    
    # Convert %ID/g to organ activity (Bq)
    activity_data = []
    for t_h, pidg in time_data:
        # A(t) = (pidg / 100) * organ_mass_g * injected_activity_Bq
        activity = (pidg / 100.0) * organ_mass_g * injected_activity_Bq
        activity_data.append((t_h * 3600, activity))  # convert to seconds
    
    if len(activity_data) == 1:
        # Single time point: assume only physical decay
        decay_const = math.log(2) / half_life_phys
        return activity_data[0][1] / decay_const
    
    # Trapezoidal integration
    cumulated = 0.0
    for i in range(1, len(activity_data)):
        dt = activity_data[i][0] - activity_data[i - 1][0]
        avg_a = 0.5 * (activity_data[i][1] + activity_data[i - 1][1])
        cumulated += avg_a * dt
    
    # Tail: physical decay beyond last measured point
    last_a = activity_data[-1][1]
    decay_const = math.log(2) / half_life_phys
    tail = last_a / decay_const
    
    return cumulated + tail


def calculate_doses(
    biodist: Dict[str, List[Tuple[float, float]]],
    nuclide: str,
    injected_activity_Bq: float,
    organ_masses: Dict[str, float] = None,
) -> Dict[str, dict]:
    """Calculate organ absorbed doses from biodistribution data."""
    if organ_masses is None:
        organ_masses = DEFAULT_ORGAN_MASSES
    
    s_values = S_VALUES_SELF.get(nuclide)
    if not s_values:
        return {"error": f"No S-values for {nuclide}. Available: {list(S_VALUES_SELF.keys())}"}
    
    half_life = HALF_LIVES.get(nuclide)
    if not half_life:
        return {"error": f"No half-life for {nuclide}. Available: {list(HALF_LIVES.keys())}"}
    
    results = {}
    
    for organ, time_course in biodist.items():
        organ_key = organ.replace(" ", "_")
        mass_g = organ_masses.get(organ_key, organ_masses.get(organ, None))
        
        cumulated_activity = None
        dose_Gy = None
        s_self = s_values.get(organ_key, None)
        
        if mass_g and s_self:
            cumulated_activity = compute_cumulated_activity_per_organ(
                time_course, half_life, mass_g, injected_activity_Bq
            )
            dose_Gy = cumulated_activity * s_self
        
        results[organ] = {
            "cumulated_activity_Bq_s": round(cumulated_activity, 1) if cumulated_activity else None,
            "S_value_self_Gy_per_Bq_s": s_self,
            "absorbed_dose_Gy": round(dose_Gy, 4) if dose_Gy else None,
            "organ_mass_g": mass_g,
            "organ_key": organ_key,
        }
    
    return results


def main():
    parser = argparse.ArgumentParser(description="Dosimetry calculator for peptide radiopharmaceuticals")
    parser.add_argument("--input", required=True, help="Input biodistribution CSV (Time_h, Organ, PercentID_per_g)")
    parser.add_argument("--nuclide", required=True, help="Radionuclide (e.g., Lu-177, Y-90, Ac-225, Ga-68)")
    parser.add_argument("--injected_activity", type=float, required=True, help="Injected activity in Bq")
    parser.add_argument("--organ-masses", default=None, help="JSON file with organ masses (g)")
    parser.add_argument("--output", default=None, help="Output JSON file (default: stdout)")
    args = parser.parse_args()
    
    biodist = read_biodistribution(args.input)
    if not biodist:
        print("Error: No valid data in input file", file=sys.stderr)
        sys.exit(1)
    
    organ_masses = DEFAULT_ORGAN_MASSES
    if args.organ_masses:
        with open(args.organ_masses) as f:
            organ_masses = json.load(f)
    
    results = calculate_doses(biodist, args.nuclide, args.injected_activity, organ_masses)
    
    output = json.dumps(results, indent=2)
    if args.output:
        with open(args.output, "w") as f:
            f.write(output)
        print(f"Dosimetry results written to {args.output}")
    else:
        print(output)


if __name__ == "__main__":
    main()
