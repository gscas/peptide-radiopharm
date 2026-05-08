---
name: peptide-radiopharm
description: |
  Peptide radiopharmaceutical drug design and development skill.
  Use when: (1) Designing peptide-based radiotracers or radiotherapeutics,
  (2) Selecting radionuclides (diagnostic vs therapeutic, half-life matching),
  (3) Choosing chelators for radiometal labeling (DOTA, NOTA, DFO, HYNIC, etc.),
  (4) Planning radiolabeling protocols and QC methods,
  (5) Evaluating in vitro properties (binding affinity, cell uptake, serum stability),
  (6) Analyzing biodistribution and dosimetry,
  (7) Designing theranostic pairs (diagnostic + therapeutic),
  (8) Optimizing peptide sequences for in vivo performance,
  (9) Calculating absorbed dose from biodistribution data,
  (10) Planning clinical translation of peptide radiopharmaceuticals.
  Covers peptide design, radiochemistry, preclinical evaluation, dosimetry, and clinical development.
---

# Peptide Radiopharmaceutical Drug Design & Development

This skill covers the full pipeline of peptide-based radiopharmaceutical development:
from peptide design and radionuclide selection through radiolabeling, preclinical
evaluation, dosimetry, and clinical translation.

## Prerequisites

- Python 3.8+ (for scripts)
- Basic understanding of nuclear medicine and peptide chemistry
- Access to OLINDA/EXM or IDAC-Dose for clinical dosimetry (optional)

## Workflow Overview

1. Target Selection & Peptide Design
2. Radionuclide Selection
3. Chelator / Prosthetic Group Selection
4. Radiolabeling & QC
5. In Vitro Evaluation
6. In Vivo Preclinical Evaluation
7. Dosimetry
8. Theranostic Pair Design
9. Clinical Translation
10. GMP Production & Regulatory

---

## 1. Target Selection & Peptide Design

### Target Selection Criteria

| Criterion | Consideration |
|---|---|
| Expression level | >10,000 receptors/cell for imaging; >50,000 for therapy |
| Tumor-to-normal ratio | >5:1 for imaging; >10:1 preferred for therapy |
| Internalization | Rapid internalization desired for radiotherapy |
| Shedding | Soluble receptor may cause background |
| Homology | Off-target binding to related receptors |

### Peptide Design Strategies

**Linear peptides:**
- Usually 6-20 amino acids
- High flexibility → often need stabilization
- Examples: RGD, bombesin(7-14), exendin-4 fragments

**Cyclic peptides:**
- Disulfide bridge (Cys-Cys): classical, reducible in vivo
- Head-to-tail cyclization: improved stability
- Side-chain-to-side-chain (lactam, triazole click): stable, versatile
- Hydrocarbon stapling: α-helix stabilization, cell penetration

**Stabilization modifications:**

| Modification | Effect | Example |
|---|---|---|
| D-amino acid substitution | Protease resistance | D-Phe, D-Tyr |
| N-methylation | Protease resistance + membrane permeability | N-Me-Ala |
| PEGylation (PEG2-4k) | Prolonged circulation, reduced kidney uptake | PEG2k at N-terminus |
| Glycosylation | Altered pharmacokinetics | sialic acid capping |
| Albumin binding moiety | Extended half-life | Evans blue derivative, fatty acid |
| C-terminal amide | Prevent carboxypeptidase degradation | -CONH2 |
| N-terminal acetylation | Prevent aminopeptidase degradation | Ac- |
| Truncation | Remove non-essential residues | Octreotide → [Tyr3]-octreotide |

**Key peptide scaffolds in radiopharmacy:**

| Scaffold | Target | Clinical examples |
|---|---|---|
| Somatostatin analog | SSTR2 | ⁶⁸Ga/¹⁷⁷Lu-DOTATATE, ⁶⁸Ga-DOTATOC |
| RGD | αvβ3 integrin | ⁶⁸Ga/¹⁸F-RGD, ⁹⁹ᵐTc-NC100692 |
| PSMA inhibitor | PSMA | ⁶⁸Ga/¹⁷⁷Lu-PSMA-617, ¹⁸F-DCFPyL |
| Bombesin analog | GRPR | ⁶⁸Ga-RM2, ¹⁷⁷Lu-AMBA |
| Exendin analog | GLP-1R | ⁶⁸Ga/¹⁷⁷Lu-Exendin-4 |
| CCK/gastrin analog | CCK2R | ⁶⁸Ga/¹⁷⁷Lu-PP-F11N |
| CXCR4 antagonist | CXCR4 | ⁶⁸Ga/¹⁷⁷Lu-Pentixafor |
| Melanocortin analog | MC1R | ⁶⁸Ga/¹⁷⁷Lu-DOTA-MSH |
| NT analog | NTR1 | ⁶⁸Ga/¹⁷⁷Lu-DOTA-NT-20.3 |
| Fibroblast activation | FAP | ⁶⁸Ga/¹⁷⁷Lu-FAPI-04, FAPI-46 |

### Peptide Property Calculation

Use `scripts/peptide_calc.py` to compute:
- Molecular weight (from sequence)
- Isoelectric point (pI)
- Net charge at pH 7.4
- Hydrophobicity index (GRAVY score)
- Instability index

```
python3 peptide_calc.py --sequence "CFYWKLLSNC" --modifications "N-term:DOTA,C-term:amide"
```

---

## 2. Radionuclide Selection

### Diagnostic Radionuclides (Imaging)

| Nuclide | Emission | Half-life | Energy (keV) | Imaging | Production | Key Use |
|---|---|---|---|---|---|---|
| ⁶⁸Ga | β⁺ | 67.7 min | 511 (β⁺), 1077 (γ) | PET | ⁶⁸Ge/⁶⁸Ga generator | SSTR, PSMA, FAP imaging |
| ¹⁸F | β⁺ | 109.8 min | 511 (β⁺) | PET | Cyclotron | FDG, FAPI, PSMA |
| ⁶⁴Cu | β⁺ | 12.7 h | 511 (β⁺) | PET | Cyclotron | Antibody/peptide, longer PK |
| ⁴⁴Sc | β⁺ | 3.97 h | 511 (β⁺) | PET | Cyclotron/⁴⁴Ti gen | Theranostic match with ⁴⁷Sc |
| ⁹⁹ᵐTc | γ | 6.0 h | 140 | SPECT | ⁹⁹Mo/⁹⁹ᵐTc gen | SPECT imaging, widely available |
| ¹¹¹In | γ | 2.80 d | 171, 245 | SPECT | Cyclotron | Longer PK tracers, dosimetry |
| ¹²³I | γ | 13.2 h | 159 | SPECT | Cyclotron | Iodinated peptides |
| ⁸⁹Zr | β⁺ | 78.4 h | 511 (β⁺) | PET | Cyclotron | Long-circulating (antibodies) |

### Therapeutic Radionuclides

| Nuclide | Emission | Half-life | Energy (MeV) | Range | Production | Key Use |
|---|---|---|---|---|---|---|
| ¹⁷⁷Lu | β⁻ | 6.65 d | 0.133 (max 0.497) | ~2 mm | Reactor | PRRT (SSTR, PSMA, FAP) |
| ⁹⁰Y | β⁻ | 2.67 d | 0.935 (max 2.28) | ~11 mm | Generator/Reactor | Large tumors, PRRT |
| ²²⁵Ac | α | 9.92 d | 5.8 (net 4 α) | 50-80 μm | Generator (²²⁹Th) | Micro-metastases, α-PRT |
| ²¹²Pb | α (via ²¹²Bi) | 10.6 h | 6.1 (via ²¹²Bi) | 50-80 μm | Generator (²²⁸Th) | α-PRT with in vivo generator |
| ²²³Ra | α | 11.4 d | 5.0-7.5 | 50-80 μm | Generator | Bone mets (Xofigo) |
| ¹³¹I | β⁻ | 8.02 d | 0.182 (max 0.606) | ~2 mm | Reactor | Iodinated peptides, thyroid |
| ⁴⁷Sc | β⁻ | 3.35 d | 0.143 (max 0.600) | ~0.5 mm | Cyclotron/Reactor | Theranostic with ⁴⁴Sc |
| ⁶⁷Cu | β⁻ | 2.58 d | 0.141 (max 0.476) | ~0.6 mm | Cyclotron | Theranostic with ⁶⁴Cu |
| ¹⁶¹Tb | β⁻ | 6.90 d | 0.154 (max 0.593) | ~0.7 mm | Reactor | Theranostic with ¹⁵²Tb |
| ²²⁷Th | α | 18.7 d | 5.9 | 50-80 μm | Generator | Antibody/peptide α-PRT |

### Theranostic Pair Matching

| Diagnostic | Therapeutic | Rationale |
|---|---|---|
| ⁶⁸Ga | ¹⁷⁷Lu | Same chelator (DOTA), widely used |
| ⁶⁸Ga | ²²⁵Ac | Same chelator (DOTA), emerging |
| ⁶⁸Ga | ⁹⁰Y | Same chelator, cross-fire effect |
| ⁶⁴Cu | ⁶⁷Cu | Same element, identical chemistry |
| ⁴⁴Sc | ⁴⁷Sc | Same element, ideal theranostic |
| ¹⁵²Tb | ¹⁶¹Tb | Same element, matched half-lives |
| ¹¹¹In | ¹⁷⁷Lu | Similar chemistry, dosimetry prediction |
| ⁹⁹ᵐTc | ¹⁸⁶Re | Same group, similar chemistry |
| ¹²⁵I | ¹³¹I | Same element, identical labeling |
| ⁸⁹Zr | ²²⁷Th | Long half-life pair for antibodies |

### Half-life Matching Principle

The radionuclide half-life should match the biological half-life of the peptide:

- **Fast-clearing peptides** (t½ < 2h): use ⁶⁸Ga, ¹⁸F, ⁹⁹ᵐTc
- **Medium-clearing peptides** (t½ 2-24h): use ⁶⁴Cu, ⁶⁷Cu, ⁴⁴Sc, ⁴⁷Sc
- **Slow-clearing peptides** (t½ > 24h): use ⁸⁹Zr, ¹¹¹In, ²²⁷Th

---

## 3. Chelator / Prosthetic Group Selection

### Chelators for Radiometals

| Chelator | Metal preference | Coordination | Labeling conditions | Stability |
|---|---|---|---|---|
| DOTA | ⁶⁸Ga, ¹⁷⁷Lu, ⁹⁰Y, ²²⁵Ac, ⁶⁴Cu, ⁴⁷Sc | Macrocyclic, 8-dentate | 80-100°C, pH 4-5, 10-15 min | High (macrocyclic effect) |
| NOTA | ⁶⁸Ga, ⁶⁴Cu, ⁴⁴Sc | Macrocyclic, 6-dentate | RT-40°C, pH 4-5, 5-10 min | Very high for Ga |
| NODAGA | ⁶⁸Ga, ⁶⁴Cu | Macrocyclic, 6+1 dentate | RT, pH 4-5, 5 min | Very high for Ga |
| DFO (desferrioxamine) | ⁸⁹Zr, ²²⁷Th | Linear, hexadentate | RT, pH 6-7, 30-60 min | Moderate (Zr8-coord preferred) |
| DFO* (DFO-macrocyclic) | ⁸⁹Zr | Macrocyclic, 8-dentate | RT, pH 6-7, 30 min | Superior to DFO |
| HYNIC | ⁹⁹ᵐTc | Bifunctional, needs coligand | RT, pH 5-7, 15 min | Good |
| MAG3 | ⁹⁹ᵐTc | Linear, N2S2 | RT-boiling, pH 7-9 | Moderate |
| Pycup | ⁶⁴Cu, ⁶⁷Cu | Cross-bridged macrocycle | RT, pH 5-6, 10 min | Excellent for Cu |
| CHX-A''-DTPA | ¹⁷⁷Lu, ⁹⁰Y, ²²⁵Ac | Linear, 8-dentate | RT, pH 5-6, 15 min | Good |
| Macropa | ²²⁵Ac, ²²⁷Th | Macrocyclic, 10-dentate | RT, pH 5-6, 5 min | Superior for Ac |
| 3p-C-NETA | ¹⁷⁷Lu, ⁹⁰Y, ²²⁵Ac | Hybrid macrocyclic | RT, pH 5.5, 5 min | Excellent |

### Chelator Selection Decision Tree

```
Which radiometal?
├── ⁶⁸Ga
│   ├── Fast labeling needed (<10 min)? → NOTA or NODAGA
│   └── Need ¹⁷⁷Lu theranostic? → DOTA (shared chelator)
├── ¹⁷⁷Lu
│   ├── Standard → DOTA
│   └── Need room-temp labeling? → CHX-A''-DTPA
├── ⁹⁰Y
│   └── DOTA or CHX-A''-DTPA
├── ²²⁵Ac
│   ├── Macropa (best stability)
│   └── DOTA (proven, but less stable for Ac)
├── ⁸⁹Zr
│   ├── DFO* (best stability)
│   └── DFO (standard, acceptable for short studies)
├── ⁶⁴Cu/⁶⁷Cu
│   ├── Pycup (best stability)
│   ├── NOTA/NODAGA (good)
│   └── DOTA (suboptimal for Cu)
├── ⁹⁹ᵐTc
│   ├── HYNIC (most common for peptides)
│   └── MAG3
└── ⁴⁷Sc
    └── DOTA or NOTA
```

### Prosthetic Groups for ¹⁸F and Radioiodine

**¹⁸F labeling (for peptides):**
| Method | Temperature | Time | Notes |
|---|---|---|---|
| Al¹⁸F-NOTA | RT-100°C | 15 min | Direct, one-step with NOTA chelate |
| ¹⁸F-SiFA | RT | 5 min | Isotopic exchange, very fast |
| ¹⁸F-SFB (prosthetic) | 2-step | 60 min | Activated ester, conjugates to Lys |
| ¹⁸F-Click (azide-alkyne) | RT | 30 min | Bioorthogonal, modular |
| ¹⁸F-TFP ester | RT-40°C | 20 min | Improved SFB variant |

**Radioiodine labeling:**
| Method | Nuclides | Oxidant | Conditions |
|---|---|---|---|
| Direct iodination (Iodogen) | ¹²³I, ¹²⁵I, ¹³¹I | Iodogen | RT, 10 min, pH 7-8 |
| Direct iodination (Chloramine-T) | same | Chloramine-T | RT, 1-2 min |
| Iodination via prosthetic (SIB) | same | Pre-labeled | RT, conjugation 30 min |
| Cu-mediated radioiodination | same | Cu²⁺ | 45°C, 30 min |

---

## 4. Radiolabeling & QC

### Typical Radiolabeling Protocols

**⁶⁸Ga-DOTA-peptide:**
```
1. Elute ⁶⁸Ga from generator (0.1M HCl, 5-7 mL)
2. Pre-purify (optional): fractionation or anion exchange
3. Adjust pH to 3.5-4.0 with sodium acetate buffer
4. Add peptide precursor (5-20 nmol)
5. Heat at 95°C for 10-15 min
6. Cool, adjust pH to 7
7. Purify: C18 Sep-Pak or HPLC
8. Formulate in saline with ascorbic acid
```

**¹⁷⁷Lu-DOTA-peptide:**
```
1. Add ¹⁷⁷LuCl₃ (no-carrier-added) to reaction vial
2. Add sodium acetate buffer (pH 5.0)
3. Add peptide precursor (10-50 nmol)
4. Heat at 95°C for 15-30 min (or RT for 30-60 min for some chelators)
5. Cool, quench with DTPA (scavenge free Lu)
6. Purify if needed (C18 Sep-Pak)
7. Formulate in saline + ascorbic acid (prevent radiolysis)
```

**²²⁵Ac-DOTA-peptide (or Macropa):**
```
1. Add ²²⁵Ac in HCl to reaction vial
2. Add ammonium acetate buffer (pH 5.5-6.0)
3. Add peptide precursor (20-100 nmol, higher due to low Ac specific activity)
4. Heat at 70-80°C for 15-30 min (DOTA) or RT for 5-10 min (Macropa)
5. Cool, quench with DTPA
6. Purify: C18 Sep-Pak or HPLC
7. Formulate in saline + ascorbic acid (critical for Ac due to α-radiolysis)
8. QC: check ²²¹Fr/²¹⁷At daughters (α spectrum)
```

### QC Requirements

| Test | Method | Specification |
|---|---|---|
| Radiochemical purity (RCP) | radio-HPLC or radio-TLC | >95% (clinical), >90% (preclinical) |
| Radiochemical identity | HPLC co-injection with reference | Retention time match |
| Radionuclidic purity | γ-spectroscopy | As per Ph. Eur./USP |
| Chemical purity | HPLC-UV | < specified limits |
| Specific activity | Calculation or measurement | MBq/nmol |
| pH | pH strip/meter | 4-8 (injectable) |
| Sterility | Filter (0.22 μm) / sterility test | Sterile |
| Endotoxin | LAL test | <175 EU/V (clinical) |
| Free radiometal | ITLC (with DTPA) | <2% |
| Colloids | ITLC or filter trap | <5% |

---

## 5. In Vitro Evaluation

### Binding Affinity

| Assay | Method | Readout | Typical Target |
|---|---|---|---|
| Saturation binding | Increasing conc. on cells/tissue | Kd, Bmax | Kd < 10 nM (imaging), < 1 nM (therapy) |
| Competition binding | Compete with labeled standard | IC50, Ki | Ki < 100 nM |
| Kinetics (on/off rate) | Real-time (Biacore/SPR) or cell-based | kon, koff | Slow off-rate preferred |

**Key formulas:**
```
Ki = IC50 / (1 + [L]/Kd_label)    (Cheng-Prusoff equation)
Bmax: pmol/mg protein or sites/cell
```

### Cell-Based Assays

| Assay | Purpose | Method |
|---|---|---|
| Cell uptake | Receptor-mediated internalization | Incubate, wash, lyse, count |
| Internalization | Surface vs internalized fraction | Acid wash (strip surface) |
| Efflux/blocked washout | Retention of radioactivity | Time-course after washout |
| Serum stability | In vitro metabolic stability | Incubate in serum, HPLC at time points |
| Serum protein binding | Free fraction | Size exclusion or ultrafiltration |
| Blocking studies | Receptor specificity | Co-incubate with excess cold compound |

### Serum Stability Targets

| Nuclide | Target stability | Time point |
|---|---|---|
| ⁶⁸Ga | >90% intact | 2 h in human serum, 37°C |
| ¹⁷⁷Lu | >85% intact | 24 h |
| ²²⁵Ac | >80% intact | 7 d |
| ⁹⁹ᵐTc | >90% intact | 4 h |
| ¹⁸F | >95% intact | 2 h |

---

## 6. In Vivo Preclinical Evaluation

### Animal Models

| Model | Use | Notes |
|---|---|---|
| Xenograft (SCID/nu) | Tumor targeting | Subcutaneous or orthotopic |
| PDX | More clinically relevant | Patient-derived |
| Syngeneic | Immunocompetent | Immune tumor microenvironment |
| Transgenic | Spontaneous tumors | GEMM |

### Biodistribution Study Design

```
Groups: n=3-5 per time point
Time points: 0.5, 1, 2, 4, 24 h (adjust per nuclide half-life)
Dose: 0.1-1 nmol, 1-10 MBq per animal (for SPECT/PET also imaging)
Blocking: +50-fold excess cold compound (confirm receptor specificity)
Organs: blood, tumor, heart, lung, liver, spleen, kidney, intestine,
        muscle, bone, brain, stomach, pancreas, adrenals, thyroid
Readout: %ID/g (percent injected dose per gram)
```

### Key Metrics

| Metric | Formula | Target |
|---|---|---|
| Tumor uptake | %ID/g tumor | >5 %ID/g at 1h |
| Tumor-to-blood | Tumor %ID/g / Blood %ID/g | >5 |
| Tumor-to-muscle | Tumor %ID/g / Muscle %ID/g | >10 |
| Tumor-to-kidney | Tumor %ID/g / Kidney %ID/g | >0.5 (therapy), >1 (imaging) |
| Tumor-to-liver | Tumor %ID/g / Liver %ID/g | >1 |

### Kidney Uptake Reduction Strategies

Renal uptake is the dose-limiting factor for peptide radiotherapy.

| Strategy | Mechanism | Typical Reduction |
|---|---|---|
| Lysine/arginine infusion | Block proximal tubule reabsorption | 30-50% |
| Gelofusine (Gelatin) | Competitive reabsorption inhibition | 20-40% |
| Succinylated gelatin | Same as above | 30-50% |
| Albumin fragment | Block megalin-mediated uptake | 40-60% |
| PEGylation | Increase size, reduce filtration | 30-70% (depends on PEG size) |
| Cleavable linker (GMG/AG) | Enzymatic cleavage in kidney | 50-80% |
| Negatively charged linker | Reduce tubular reabsorption | 20-40% |
| Reducible disulfide (Cys-Cys) | Intracellular reduction + efflux | Variable |

---

## 7. Dosimetry

### Preclinical Dosimetry (from Biodistribution)

Use `scripts/dosimetry_calc.py` to convert biodistribution data to absorbed doses.

**MIRD formalism:**
```
D_target = Σ_s Ã_s × S(target ← source)

Where:
  D_target = absorbed dose to target organ [Gy]
  Ã_s = cumulated activity in source organ [Bq·s]
  S(target ← source) = S-value [Gy/Bq·s]
```

**Cumulated activity calculation:**
```
Ã = ∫ A(t) dt

From biodistribution %ID/g data:
A(t) = %ID/g(t) × organ_mass × injected_activity

Fitting: mono- or bi-exponential, or trapezoidal integration
```

### Clinical Dosimetry Methods

| Method | Data Required | Accuracy | Complexity |
|---|---|---|---|
| Single-time-point (STP) | 1 SPECT/PET + population kinetics | Moderate | Low |
| Multi-time-point | 3-5 SPECT/PET scans | High | High |
| Whole-body counting | Gamma probe measurements | Moderate | Low |
| Blood-based | Serial blood samples | Moderate (red marrow) | Medium |

### S-Value Tables

See `references/s_values.md` for organ S-values for common nuclides (⁶⁸Ga, ¹⁷⁷Lu, ⁹⁰Y, ²²⁵Ac, ⁹⁹ᵐTc).

### Dose Limits for Therapy

| Organ | Dose limit (Gy) | Basis |
|---|---|---|
| Kidneys | 23 (⁹⁰Y), 28 (¹⁷⁷Lu) | BED-adjusted |
| Red marrow | 2 | Hematologic toxicity |
| Liver | 30 | Radiation hepatitis |
| Salivary glands | 30 (¹⁷⁷Lu) | Xerostomia |
| Bowel | 35 | Ulceration |
| Bone surfaces | 50 | Fracture risk |
| Whole body | 2 (non-stochastic) | ICRP |

### BED (Biologically Effective Dose)

```
BED = D × (1 + D/(α/β × d))

For kidney (α/β = 2-3 Gy):
  BED = D × (1 + D/(2.5 × d))

Where D = total dose, d = dose per cycle
BED limit for kidney: ~40-45 Gy (corresponds to 28 Gy physical for ¹⁷⁷Lu 4-cycle)
```

---

## 8. Theranostic Pair Design

### Design Principles

1. **Same targeting vector** — identical peptide-chelator conjugate
2. **Matched pharmacokinetics** — diagnostic predicts therapeutic biodistribution
3. **Chelator compatibility** — same chelator for both nuclides (ideal)
4. **Half-life matching** — nuclide t½ ≈ biological t½ of peptide

### Workflow

```
Step 1: Patient imaging with diagnostic (e.g., ⁶⁸Ga-DOTATATE PET/CT)
Step 2: Quantify tumor and organ uptake (SUV → %ID/g)
Step 3: Dosimetry calculation (predict ¹⁷⁷Lu organ doses)
Step 4: Verify tumor uptake sufficient (>SUVmax 20 for SSTR)
Step 5: Calculate maximum tolerable activity
Step 6: Treat with ¹⁷⁷Lu-DOTATATE (typically 7.4 GBq × 4 cycles)
Step 7: Post-therapy SPECT for verification dosimetry
Step 8: Restage with ⁶⁸Ga PET after 2-4 cycles
```

### Theranostic Concierge Criteria

Before therapy, diagnostic scan must confirm:
- [ ] Sufficient tumor uptake (Krenning score ≥3 for SSTR)
- [ ] Acceptable kidney uptake (can be managed)
- [ ] Bone marrow reserve adequate
- [ ] No contraindications (renal failure, severe myelosuppression)
- [ ] Hepatic metastases not dominating clearance

---

## 9. Clinical Translation

### Phase I Design for Peptide Radiopharmaceuticals

```
Design: 3+3 dose escalation
Starting dose: Based on preclinical NOAEL with 100× safety margin
Dose levels: e.g., 1.85, 3.7, 5.55, 7.4 GBq (for ¹⁷⁷Lu)
Cohort size: 3-6 patients per level
DLT window: 4-8 weeks post-injection
DLT criteria: Grade 4 hematologic (4 wk), Grade 3+ non-hematologic
Expansion cohort: At RP2D, 10-15 patients
Primary endpoint: Safety (CTCAE), MTD/RP2D
Secondary endpoints: PK, dosimetry, tumor response (RECIST 1.1)
```

### Key Regulatory References

- FDA Guidance: "Clinical Non-Oncology Indications for Radiopharmaceuticals"
- EMA Guideline: "Radiopharmaceuticals — non-clinical studies"
- ICRP Publication 128: Radiation Dose to Patients from Radiopharmaceuticals
- Ph. Eur. Monographs for ⁶⁸Ga, ¹⁷⁷Lu, ⁹⁹ᵐTc radiopharmaceuticals

---

## 10. GMP Production & Regulatory

### GMP Radiopharmacy Requirements

| Aspect | Requirement |
|---|---|
| Facility | ISO Class 5 hot cell, ISO Class 7 background |
| Synthesis module | GMP-qualified (e.g., Gallia, Modular Lab) |
| Generator | GMP-certified (⁶⁸Ge/⁶⁸Ga) |
| Starting materials | GMP-grade peptide precursor |
| Aseptic processing | Terminal sterilization by 0.22 μm filter |
| QC release | RCP, sterility, endotoxin, pH, radionuclidic purity |
| Documentation | Batch record, SOP, deviation handling |
| Shelf life | Based on RCP stability (typically 2-4 h for ⁶⁸Ga) |

### Regulatory Pathways

| Region | Pathway | Typical Timeline |
|---|---|---|
| US | IND → NDA/BLA | 6-12 mo IND, 8-12 mo review |
| EU | Clinical Trial Application → MAA | 6-12 mo CTA, 12-18 mo review |
| China | IND → NDA (NMPA) | Similar timeline |
| Japan | Clinical Trial Notification → J-NDA | Similar timeline |

---

## Quick Reference: Decision Checklist

### New Peptide Radiopharmaceutical Project

- [ ] Target validated? (expression, tumor-to-normal ratio)
- [ ] Lead peptide selected? (affinity, stability, internalization)
- [ ] Radionuclide chosen? (diagnostic first → theranostic)
- [ ] Chelator matched? (metal-chelator compatibility)
- [ ] Labeling optimized? (RCP >95%, SA sufficient)
- [ ] In vitro profiled? (Kd, IC50, serum stability, cell uptake)
- [ ] Biodistribution acceptable? (tumor >5%ID/g, kidney manageable)
- [ ] Dosimetry calculated? (organs within limits)
- [ ] Toxicology done? (single/multi-dose, no-carrier-added)
- [ ] GMP production ready? (module, QC methods, batch record)
- [ ] Clinical protocol approved? (Phase I design, ethics)
