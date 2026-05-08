# S-Value Reference Tables for Peptide Radiopharmaceutical Dosimetry

S-values (Gy/Bq·s) for reference adult male (ICRP 133 / OLINDA/EXM).

## Self-dose S-values (target ← target)

### Lu-177 (β⁻, t½ = 6.65 d)

| Target | S(Gy/Bq·s) |
|---|---|
| Adrenals | 3.46E-09 |
| Bone surfaces | 2.25E-10 |
| Brain | 2.41E-11 |
| Breasts | 2.95E-11 |
| Gallbladder wall | 4.29E-10 |
| GI tract - Stomach wall | 2.15E-10 |
| GI tract - Small intestine | 1.44E-10 |
| GI tract - Colon (UL) | 1.53E-10 |
| GI tract - Colon (LL) | 1.61E-10 |
| Heart wall | 1.78E-10 |
| Kidneys | 2.06E-10 |
| Liver | 1.19E-10 |
| Lungs | 4.35E-11 |
| Muscle | 3.36E-11 |
| Oesophagus | 5.03E-10 |
| Ovaries | 1.78E-10 |
| Pancreas | 4.69E-10 |
| Red marrow | 1.30E-10 |
| Salivary glands | 2.10E-10 |
| Skin | 1.21E-11 |
| Spleen | 2.84E-10 |
| Testes | 2.48E-10 |
| Thymus | 4.90E-10 |
| Thyroid | 2.03E-10 |
| Urinary bladder wall | 3.66E-10 |
| Uterus | 2.03E-10 |
| Whole body | 3.16E-12 |

### Y-90 (β⁻, t½ = 2.67 d)

| Target | S(Gy/Bq·s) |
|---|---|
| Kidneys | 8.54E-10 |
| Liver | 4.71E-10 |
| Spleen | 1.14E-09 |
| Red marrow | 5.43E-10 |
| Bone surfaces | 9.15E-10 |
| Thyroid | 8.00E-10 |
| Adrenals | 1.39E-09 |
| Whole body | 1.26E-11 |

### Ga-68 (β⁺, t½ = 67.7 min)

| Target | S(Gy/Bq·s) |
|---|---|
| Kidneys | 5.49E-10 |
| Liver | 2.85E-10 |
| Spleen | 7.50E-10 |
| Red marrow | 2.56E-10 |
| Whole body | 7.69E-12 |

### Ac-225 (α, t½ = 9.92 d)

| Target | S(Gy/Bq·s) |
|---|---|
| Kidneys | 6.68E-09 |
| Liver | 3.67E-09 |
| Spleen | 9.07E-09 |
| Red marrow | 4.06E-09 |
| Bone surfaces | 7.24E-09 |
| Whole body | 1.02E-10 |

**Note:** Ac-225 S-values include contributions from daughters (Fr-221, At-217, Bi-213, Po-213, Tl-209, Pb-209) assuming local deposition.

### Tc-99m (γ, t½ = 6.0 h)

| Target | S(Gy/Bq·s) |
|---|---|
| Kidneys | 1.83E-10 |
| Liver | 1.02E-10 |
| Spleen | 2.54E-10 |
| Red marrow | 7.96E-11 |
| Thyroid | 1.73E-10 |
| Whole body | 2.20E-12 |

### In-111 (γ, t½ = 2.80 d)

| Target | S(Gy/Bq·s) |
|---|---|
| Kidneys | 4.31E-10 |
| Liver | 2.46E-10 |
| Spleen | 6.06E-10 |
| Red marrow | 2.32E-10 |
| Whole body | 5.71E-12 |

## Cross-fire S-values

For cross-organ dose contributions (e.g., D(Kidney ← Liver)), full S-value matrices
are available in OLINDA/EXM software. For peptide radiopharmaceuticals, self-dose
typically dominates (>90% of total organ dose) for most organs.

Key cross-fire pairs to consider:
- Kidney dose contribution from Liver (significant for hepatobiliary-clearing tracers)
- Red marrow dose from Bone surfaces (important for bone-seeking agents)
- Thyroid dose from surrounding tissues (important for radioiodine)

## References

1. ICRP Publication 133: The ICRP Computational Framework for Internal Dose Assessment
2. Stabin MG, Siegel JA. RADAR dose estimate report. (OLINDA/EXM)
3. Bolch WE et al. MIRD Pamphlet No. 21. J Nucl Med 2009
4. https://doseinfo-radar.com/RADAR-Home.html
