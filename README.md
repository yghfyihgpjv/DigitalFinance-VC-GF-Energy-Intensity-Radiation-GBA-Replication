
# Replication Package README

This repository contains all replication materials for the paper  
*Asymmetric Digital Finance Exposure and Forecasting Energy Transitions in Polycentric Urban Systems: Beyond a Tale of Two Cores*.

Running the main script reproduces all empirical results reported in the manuscript and is intended to support the journal’s peer review and reproducibility assessment.

---

## Data Sources and References

For replication purposes, the main data sources follow established and peer-reviewed methodologies:

- **Digital Finance (DF)**:  
  Guo et al. (2020), *China Economic Quarterly*  
  DOI: 10.13821/j.cnki.ceq.2020.03.12

- **Energy Intensity (EI)**:  
  Estimated using nighttime light intensity, cross-calibrated following  
  Wu et al. (2014), *Geographical Research* (DOI: 10.11821/dlyj201404006) and  
  Yang et al. (2023), *China Industrial Economy* (DOI: 10.19581/j.cnki.ciejournal.2023.05.004)

- **Green Total Factor Productivity (GTFP)**:  
  Measured via the SBM-ML method following  
  Tone (2001), *European Journal of Operational Research* (DOI: 10.1016/S0377-2217(99)00407-5) and  
  Oh (2010), *Journal of Productivity Analysis* (DOI: 10.1007/s11123-010-0178-y)

- **Environmental Regulation Intensity (EnvReg)**:  
  Word-frequency intensity of environmental keywords in government work reports, following  
  Shao et al. (2024), *Management World*  
  DOI: 10.19744/j.cnki.11-1235/f.2024.0093
  
**Historical Information Infrastructure index（PCA）**
   Chao, X. J., Lian, Y. M., & Luo, L. K. (2021). Impact of new digital infrastructure on the high-quality development of manufacturing industry. Finance and Trade Research, 32(10), 1-13. https://doi.org/10.19337/j.cnki.34-1093/f.2021.10.001
---

## 1. File Structure

```text
/data
  /rawdata
    df_index.xls
    gf_index.xls
    vc_index.xls
    vc_micro_firm.xls
   Historical Information Infrastructure index（PCA）.xls
    测度中国数字普惠金融发展_数据编制与空间特....pdf
    gf_measurement_method.docx
    绿色金融促进城市经济高质量...272个地级市....pdf
    绿色金融对农业高质量发展的...槛效应和中介效....pdf
    绿色金融对数字经济绿色发展影响效应研究_谢非.pdf
    绿色金融与新质生产力：促进...技术创新与环境....pdf

  /processeddata
    29.dta
    processeddata.xls

/analysis
  code.do
  wconstructingcode.do

/results
  a.txt

README.md
````

---

## 2. Data Description

### (1) Raw Data (`/data/rawdata/`)

All raw datasets (except literature and methodology documents) are provided in **Excel (.xls)** format for transparency and traceability.

* **DF data** (`df_index.xls`):
  Total digital finance index and sub-indicators (coverage breadth, usage depth, digitization level).
  Measurement logic and original data sources are documented in the accompanying DF methodology PDF.

* **GF data** (`gf_index.xls`):
  Core green finance indicators. Variable definitions and construction rules are documented in `gf_measurement_method.docx`.

* **VC data**:

  * `vc_micro_firm.xls`: Micro firm-level VC records (basic aggregation unit).
  * `vc_index.xls`: Prefecture-level VC panel aggregated from micro data, following the rules described in the manuscript.

### (2) Processed Data (`/data/processeddata/`)

* **29.dta** (Stata format, core operational file):
  The only dataset directly used in regression analysis. It integrates outcome variables, digital finance indicators, institutional variables (VC, GF), control variables, fixed-effect identifiers, and pre-constructed radiation weights.

* **processeddata.xls** (Excel format, auxiliary file):
  Provided for traceability only and not used in estimation.

---

## 3. Code Description

### (1) Main Script (`/analysis/code.do`)

This is the only script required to replicate all empirical results. It:

1. Loads `29.dta` using relative paths.
2. Generates descriptive statistics.
3. Estimates baseline fixed-effects regressions with clustered standard errors.
4. Conducts spatial heterogeneity and institutional mechanism analyses.
5. Performs robustness checks and counterfactual forecasting simulations.
6. Writes all outputs to `/results/a.txt`.

### (2) Weight Construction Script (`/analysis/wconstructingcode.do`)

This script documents the construction of directed radiation weights.
All weights are already embedded in `29.dta`.
The script is provided for transparency only and does not need to be executed.

---

## 4. How to Replicate

### Data Placement

Ensure that `29.dta` remains in `data/processeddata/`.
All scripts use relative paths; no files need to be moved.

### Execution

From the repository root directory, run in Stata:

```stata
do analysis/code.do
```

### Output

All results will be saved to:

```text
/results/a.txt
```

---

## 5. Software Requirements

* **Stata 15 or later** (tested on Stata 17)
* **Microsoft Excel / WPS** (optional, for viewing raw data)
* No additional user-written Stata packages are required.

---

## 6. Additional Notes

This replication package is provided for the journal’s peer review process.

* Publication metadata will be updated upon acceptance.
* A permanent access link will be provided for long-term reproducibility.
* Only `29.dta` is required to run the analysis; other files are for documentation.

**All scripts use relative paths. Please run the code from the repository root directory.**

```


