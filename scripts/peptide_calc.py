#!/usr/bin/env python3
"""
Peptide Property Calculator for Radiopharmaceutical Design.

Usage:
  python3 peptide_calc.py --sequence "CFYWKLLSNC"
  python3 peptide_calc.py --sequence "D-Phe-Cys-Tyr-D-Trp-Lys-Thr-Cys-Thr-OH" --modifications "N-term:DOTA,C-term:amide"

Computes: MW, pI, net charge at pH 7.4, GRAVY score, instability index, cyclization sites.
"""

import argparse
import math
import sys
from typing import Dict, List, Tuple

# Amino acid data: (3-letter, 1-letter, monoisotopic MW, pKa_sidechain_or_N/C_term, hydrophobicity_KD)
# MW includes H2O loss (peptide bond). For terminal: add H (N) + OH (C) = 18.01056
AMINO_ACIDS = {
    "A": ("Ala", 71.03711, None, 1.8),
    "R": ("Arg", 156.10111, 12.48, -4.5),
    "N": ("Asn", 114.04293, None, -3.5),
    "D": ("Asp", 115.02694, 3.65, -3.5),
    "C": ("Cys", 103.00919, 8.18, 2.5),
    "E": ("Glu", 129.04259, 4.25, -3.5),
    "Q": ("Gln", 128.05858, None, -3.5),
    "G": ("Gly", 57.02146, None, -0.4),
    "H": ("His", 137.05891, 6.00, -3.2),
    "I": ("Ile", 113.08406, None, 4.5),
    "L": ("Leu", 113.08406, None, 3.8),
    "K": ("Lys", 128.09496, 10.53, -3.9),
    "M": ("Met", 131.04049, None, 1.9),
    "F": ("Phe", 147.06841, None, 2.8),
    "P": ("Pro", 97.05276, None, -1.6),
    "S": ("Ser", 87.03203, None, -0.8),
    "T": ("Thr", 101.04768, None, -0.7),
    "W": ("Trp", 186.07931, None, -0.9),
    "Y": ("Tyr", 163.06333, 10.07, -1.3),
    "V": ("Val", 99.06841, None, 4.2),
}

# D-amino acids have same MW and properties
D_PREFIX = "d-"

# Common modifications MW
MODIFICATIONS = {
    "DOTA": 404.17,
    "NOTA": 304.15,
    "NODAGA": 347.17,
    "DFO": 601.34,
    "HYNIC": 122.08,
    "PEG2k": 2000.0,
    "PEG4k": 4000.0,
    "amide": -0.98402,  # C-terminal amide: replace OH with NH2 (loss of O, gain of NH)
    "acetyl": 42.01056,  # N-terminal acetylation
    "fosfestyl": 79.96633,  # Phosphorylation
    "biotin": 226.07760,
    "fluorescein": 332.07,
    "Evans blue": 463.19,
    "Amino linker": 57.02146,
}

# pKa values for termini
PKA_N_TERM = 9.69
PKA_C_TERM = 2.34

# Instability index weights (dipeptide)
# Simplified: using a subset of the Guruprasad instability weight table
INSTABILITY_WEIGHTS = {
    ("W", "Y"): 1.0, ("Y", "W"): 1.0, ("C", "W"): 1.0, ("W", "C"): 1.0,
    ("M", "P"): 1.0, ("P", "M"): 1.0, ("G", "P"): 1.0, ("P", "G"): 1.0,
    # Default for unlisted dipeptides
    "_default": 1.0,
}


def parse_sequence(seq: str) -> List[str]:
    """Parse sequence string into list of single-letter amino acid codes.
    Handles: 'CFYWKLLSNC', 'Cys-Tyr-D-Trp...', 'c[CFYWKLLSNC]'
    """
    seq = seq.strip()
    
    # Bracket notation: c[CFYWKLLSNC] → cyclic
    if seq.startswith("c[") and seq.endswith("]"):
        inner = seq[2:-1]
        return list(inner)
    
    # Dash notation with 3-letter codes
    if "-" in seq and len(seq.split("-")[0]) > 1:
        residues = []
        for part in seq.split("-"):
            part = part.strip()
            is_d = False
            if part.upper().startswith("D-") or part.startswith("d-"):
                is_d = True
                part = part[2:]
            # 3-letter to 1-letter
            for code, (three, *_) in AMINO_ACIDS.items():
                if part.upper() == three.upper():
                    residues.append(code)
                    break
            else:
                if len(part) == 1 and part.upper() in AMINO_ACIDS:
                    residues.append(part.upper())
        return residues
    
    # Simple 1-letter code string
    residues = []
    i = 0
    while i < len(seq):
        if seq[i:i+2].lower() == "d-" and i+2 < len(seq):
            residues.append(seq[i+2].upper())
            i += 3
        elif seq[i].upper() in AMINO_ACIDS:
            residues.append(seq[i].upper())
            i += 1
        else:
            i += 1  # skip unknown chars
    
    return residues


def calculate_mw(residues: List[str], modifications: Dict[str, bool] = None) -> float:
    """Calculate molecular weight of peptide."""
    mw = 18.01056  # H2O for peptide backbone (H + OH termini)
    for res in residues:
        if res in AMINO_ACIDS:
            mw += AMINO_ACIDS[res][1]
    if modifications:
        for mod, active in modifications.items():
            if active and mod in MODIFICATIONS:
                mw += MODIFICATIONS[mod]
    return mw


def calculate_pi(residues: List[str]) -> float:
    """Estimate isoelectric point."""
    # Collect all ionizable groups
    pka_list = []
    # N-terminus (positive)
    pka_list.append(("pos", PKA_N_TERM))
    # C-terminus (negative)
    pka_list.append(("neg", PKA_C_TERM))
    
    for res in residues:
        if res == "K":
            pka_list.append(("pos", 10.53))
        elif res == "R":
            pka_list.append(("pos", 12.48))
        elif res == "H":
            pka_list.append(("pos", 6.00))
        elif res == "D":
            pka_list.append(("neg", 3.65))
        elif res == "E":
            pka_list.append(("neg", 4.25))
        elif res == "C":
            pka_list.append(("neg", 8.18))
        elif res == "Y":
            pka_list.append(("neg", 10.07))
    
    # Bisection method to find pH where net charge = 0
    ph_low, ph_high = 0.0, 14.0
    for _ in range(100):
        ph_mid = (ph_low + ph_high) / 2
        charge = 0.0
        for ion_type, pka in pka_list:
            if ion_type == "pos":
                charge += 1.0 / (1.0 + 10 ** (ph_mid - pka))
            else:
                charge -= 1.0 / (1.0 + 10 ** (pka - ph_mid))
        if charge > 0:
            ph_low = ph_mid
        else:
            ph_high = ph_mid
    return round((ph_low + ph_high) / 2, 2)


def net_charge_at_ph(residues: List[str], ph: float = 7.4) -> float:
    """Calculate net charge at given pH."""
    charge = 0.0
    # N-terminus
    charge += 1.0 / (1.0 + 10 ** (ph - PKA_N_TERM))
    # C-terminus
    charge -= 1.0 / (1.0 + 10 ** (PKA_C_TERM - ph))
    
    for res in residues:
        if res == "K":
            charge += 1.0 / (1.0 + 10 ** (ph - 10.53))
        elif res == "R":
            charge += 1.0 / (1.0 + 10 ** (ph - 12.48))
        elif res == "H":
            charge += 1.0 / (1.0 + 10 ** (ph - 6.00))
        elif res == "D":
            charge -= 1.0 / (1.0 + 10 ** (3.65 - ph))
        elif res == "E":
            charge -= 1.0 / (1.0 + 10 ** (4.25 - ph))
        elif res == "C":
            charge -= 1.0 / (1.0 + 10 ** (8.18 - ph))
        elif res == "Y":
            charge -= 1.0 / (1.0 + 10 ** (10.07 - ph))
    
    return round(charge, 3)


def gravy_score(residues: List[str]) -> float:
    """Calculate Grand Average of Hydropathy (GRAVY)."""
    if not residues:
        return 0.0
    total = sum(AMINO_ACIDS.get(r, ("", 0, None, 0.0))[3] for r in residues)
    return round(total / len(residues), 3)


def instability_index(residues: List[str]) -> float:
    """Estimate instability index (simplified).
    < 40 = stable, > 40 = unstable.
    """
    if len(residues) < 2:
        return 0.0
    total = 0.0
    for i in range(len(residues) - 1):
        dipep = (residues[i], residues[i + 1])
        total += INSTABILITY_WEIGHTS.get(dipep, INSTABILITY_WEIGHTS["_default"])
    ii = (10.0 / len(residues)) * total
    return round(ii, 2)


def cyclization_sites(residues: List[str]) -> Dict[str, List[int]]:
    """Find potential cyclization sites."""
    sites = {}
    # Cysteine pairs (disulfide)
    cys_positions = [i for i, r in enumerate(residues) if r == "C"]
    if len(cys_positions) >= 2:
        sites["disulfide"] = cys_positions
    
    # N-terminal amine (head-to-tail)
    if residues:
        sites["head_to_tail"] = [0, len(residues) - 1]
    
    # Lys + Asp/Glu (lactam)
    lys_pos = [i for i, r in enumerate(residues) if r == "K"]
    asp_glu_pos = [i for i, r in enumerate(residues) if r in ("D", "E")]
    if lys_pos and asp_glu_pos:
        sites["lactam"] = {"Lys": lys_pos, "Asp_Glu": asp_glu_pos}
    
    return sites


def main():
    parser = argparse.ArgumentParser(description="Peptide property calculator for radiopharmaceutical design")
    parser.add_argument("--sequence", required=True, help="Peptide sequence (1-letter or 3-letter codes)")
    parser.add_argument("--modifications", default=None, help="Comma-separated modifications (e.g., N-term:DOTA,C-term:amide,PEG4k)")
    parser.add_argument("--ph", type=float, default=7.4, help="pH for net charge calculation (default: 7.4)")
    args = parser.parse_args()
    
    residues = parse_sequence(args.sequence)
    if not residues:
        print("Error: Could not parse sequence", file=sys.stderr)
        sys.exit(1)
    
    mods = {}
    if args.modifications:
        for mod_str in args.modifications.split(","):
            mod_str = mod_str.strip()
            # Handle "N-term:DOTA" or just "DOTA"
            mod_name = mod_str.split(":")[-1] if ":" in mod_str else mod_str
            if mod_name in MODIFICATIONS:
                mods[mod_name] = True
    
    mw = calculate_mw(residues, mods)
    pi = calculate_pi(residues)
    charge = net_charge_at_ph(residues, args.ph)
    gravy = gravy_score(residues)
    instab = instability_index(residues)
    cyc_sites = cyclization_sites(residues)
    
    result = {
        "sequence_1letter": "".join(residues),
        "residue_count": len(residues),
        "molecular_weight_Da": round(mw, 2),
        "isoelectric_point_pI": pi,
        f"net_charge_at_pH_{args.ph}": charge,
        "GRAVY_score": gravy,
        "instability_index": instab,
        "stability_prediction": "stable" if instab < 40 else "unstable",
        "cyclization_sites": {k: str(v) for k, v in cyc_sites.items()},
        "cysteine_count": residues.count("C"),
        "charged_residues": sum(1 for r in residues if r in "KRHDE"),
        "aromatic_residues": sum(1 for r in residues if r in "FWY"),
    }
    
    if mods:
        result["modifications"] = list(mods.keys())
    
    import json
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
