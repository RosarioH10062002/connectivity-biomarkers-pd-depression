# Functional Connectivity Biomarkers of Depression in Parkinson’s Disease (Impact Scholars Program – Neuromatch Academy)

---

## Project Overview

Depression is one of the most prevalent and disabling non-motor symptoms of Parkinson’s disease (PD), affecting up to 50% of patients and significantly worsening both motor and cognitive outcomes. Despite its clinical relevance, the neural mechanisms underlying depression in PD remain incompletely understood.

This project investigates whether **resting-state functional connectivity (rs-fMRI)** can be used to derive **candidate neuroimaging biomarkers of depression in Parkinson’s disease**, leveraging large-scale neuroimaging data and machine learning techniques.

---

## Research Question

**Can resting-state functional connectivity be used to identify reliable biomarkers of depression in Parkinson’s disease?**

---

## Objectives

1. Extract functional connectivity features from rs-fMRI data in PD and control participants  
2. Identify network-level connectivity patterns associated with depression in PD  
3. Develop machine learning models to classify:
   - PD with depression  
   - PD without depression  
4. Evaluate model performance and interpretability  

---

## Dataset

- **Source:** Parkinson’s Progression Markers Initiative (PPMI)
- **Modality:** Resting-state fMRI
- **Participants:** Parkinson’s disease patients and healthy controls
- **Depression labeling:**
  - Depressed: Geriatric Depression Scale (GDS) ≥ 5  
  - Non-depressed: GDS < 5
---

## Methodology

### Preprocessing
- Standard fMRI pipelines:
  - CONN toolbox  
- Motion correction, nuisance regression, band-pass filtering, normalization  
- Brain parcellation and ROI time-series extraction  

## Tools & Libraries

- Python  
- Nilearn  
- Scikit-learn  
- PyTorch  
- NetworkX  
- fMRIPrep  


