# Replication Package README

This repository contains all replication materials for the paper *Asymmetric Digital Finance Exposure and Forecasting Energy Transitions in Polycentric Urban Systems: Beyond a Tale of Two Cores*. Running the main script reproduces all empirical results to support the journal's peer review process.
## Data Sources and References
For replication purposes, please refer to the following authoritative sources:
- **Digital Finance(DF)**: Guo et al. (2020) [DOI: 10.13821/j.cnki.ceq.2020.03.12]
-**Energy Intensity (EI)**: Estimated using nighttime light intensity, cross-calibrated following Wu et al. (2014) [DOI: 10.11821/dlyj201404006] and Yang et al. (2023) [DOI: 10.19581/j.cnki.ciejournal.2023.05.004].
- **Green Total Factor Productivity (GTFP)**: Measured via the SBM-ML method, following Tone (2001) [DOI: 10.1016/S0377-2217(99)00407-5] and Oh (2010) [DOI: 10.1007/s11123-010-0178-y].
- **Environmental Regulation Intensity (EnvReg)**: Word-frequency intensity of environmental keywords in government work reports, following Shao et al. (2024) [DOI: 10.19744/j.cnki.11-1235/f.2024.0093].


## 1. File Structure
```
/data
  /rawdata
    df_index.xls                    # Digital Finance (DF) index & sub-indicators (Excel)
    gf_index.xls                    # Green Finance (GF) index (Excel)
    vc_index.xls                    # VC prefecture-level panel (Excel)
    vc_micro_firm.xls               # VC micro firm data (basic unit of VC panel, Excel)
    测度中国数字普惠金融发展_数据编制与空间特....pdf  # DF source literature (includes DF measurement)
    gf_measurement_method.docx      # GF measurement methodology
    绿色金融促进城市经济高质量...272个地级市....pdf  # GF source literature
    绿色金融对农业高质量发展的...槛效应和中介效....pdf  # Supplementary GF literature
    绿色金融对数字经济绿色发展影响效应研究_谢非.pdf  # Supplementary GF literature
    绿色金融与新质生产力：促进...技术创新与环境....pdf  # Supplementary GF literature
  /processeddata
    29.dta                          # Core processed data for regression (actual operational file, Stata)
    processeddata.xls               # Auxiliary processed panel data (Excel, for data traceability)
/analysis
  code.do                           # Main script (run this to replicate results)
  wconstructingcode.do              # Radiation weight construction (for inspection only)
/results
  a.log                            # Full regression result log
README.md                           # This document
```


## 2. Data Description
### (1) Raw Data (/data/rawdata/)
All raw data files (except literature and methodology documents) are in **Excel (.xls)** format, ensuring easy access to original data sources:
- **DF data**: `df_index.xls` (DF total index & sub-indicators: Coverage Breadth, Usage Depth, Digitization Level). Its measurement logic and data source are fully documented in `测度中国数字普惠金融发展_数据编制与空间特....pdf`.
- **GF data**: `gf_index.xls` (GF core indicators). Its measurement methodology (variable definitions, calculation formulas) is in `gf_measurement_method.docx`; raw data sources are referenced in the GF-related literature in this folder.
- **VC data**:
  - `vc_micro_firm.xls`: Micro firm-level VC data (the "fundamental unit" for aggregating the prefecture-level panel, Excel).
  - `vc_index.xls`: Prefecture-level VC panel data (aggregated from `vc_micro_firm.xls`; aggregation rules are detailed in the paper’s Methodology section, Excel).


### (2) Processed Data (/data/processeddata/)
This folder contains two key files, with distinct formats and roles:
1. **29.dta** (Stata format, core operational file):
   - The only Stata dataset used for regression, integrating core variables (DF/GF/VC indicators, Carbon Intensity (CI) — the paper’s key outcome variable)、control variables (economic scale: lnGDP, industrial structure: IndUp/SecShare, etc.)、fixed effect identifiers (city ID, year) and pre-built radiation weights (output of `wconstructingcode.do`).
   - **Critical for replication**: This is the actual file called by the main script; see Section 4 for placement requirements.
2. **processeddata.xls** (Excel format, auxiliary file):
   - Records intermediate processing results (e.g., raw data cleaning logs, weight calculation process) for traceability. It is not directly used in regression and is provided only for transparency.


## 3. Code Description
### (1) Main Script (/analysis/code.do)
This is the only script required to replicate all empirical results, with logic aligned to the core file `29.dta`:
1. Import the core processed data (`29.dta`; note: follow Section 4 to place the file correctly).
2. Generate descriptive statistics for core variables (CI, DF, VC, GF, etc.).
3. Estimate baseline regressions (xtreg fixed effects with cluster-robust standard errors, following the FDDC framework).
4. Conduct heterogeneity tests (Market Core/G1, Policy Core/G2, Periphery/G3) and mechanism validation (DF×VC, DF×GF interactions).
5. Export all results (coefficients, SE, p-values, model fit) to `/result/a.log`.


### (2) Weight Construction Script (/analysis/wconstructingcode.do)
This script documents the full process of constructing **directed radiation weights** (nested weights, geographic/economic functional proximity weights) for the FDDC framework.  
**Key Note**: Radiation weight results are already embedded in `29.dta`. This script is provided exclusively for transparency (e.g., peer review inspection) and does not need to be run to replicate the paper’s results.


## 4. How to Replicate
### Critical Preparatory Step
Before running the script, ensure the processed data file 29.dta remains in the repository at data/processeddata/. The Stata script uses relative paths so you do not need to move any files outside the repository.
### Execution Command
In Stata (Windows/macOS compatible), execute the following command:
```stata
do analysis/code.do
```

### Result Output
After running, all empirical results (baseline regressions, heterogeneity/mechanism tests, robustness checks) will be automatically saved to:
```
/result/a.log
```


## 5. Software Requirements
- Stata 15/16/17/18 (no version compatibility issues; required to run `29.dta` and `code.do`).
- Microsoft Excel (or WPS): Required to view raw data files (.xls) if needed for data verification.
- No additional Stata modules required: All used commands are built-in, or will be automatically installed via `ssc install` in `code.do` (e.g., for regression diagnostics).


## 6. Additional Note
This replication package is provided exclusively for the journal’s peer review process.  
- Upon the paper’s acceptance, the citation will be updated with formal publication details (volume, issue, pages, DOI).  
- A permanent access link for the replication materials (e.g., archived GitHub repository) will be added to ensure long-term reproducibility.  
- For data format clarification: Only `29.dta` is in Stata format (core operational file); all other data files are in Excel (.xls) format for easy access to raw information.
