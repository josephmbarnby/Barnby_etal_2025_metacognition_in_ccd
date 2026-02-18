# Corpus Callosum Dysgenesis impairs metacognition: evidence from multi-modality and multi-cohort replications

**Authors:** Joseph M. Barnby*¹², Reinder Dean*³, H. Burgess³, Peter Dayan ⌇ ⁴⁵, Linda J. Richards ⌇ ³

\* Joint first author  
⌇ Joint senior author

¹ Institute of Psychiatry, Psychology and Neuroscience, King’s College London, UK  
² Centre for AI and Machine Learning, Edith Cowan University, WA, Australia  
³ Department of Neuroscience, Washington University in St. Louis Medical School, St Louis, MO, USA  
⁴ Max Planck Institute for Biological Cybernetics, Tübingen, Germany  
⁵ University of Tübingen, Germany

---

## 📌 Overview

This repository contains the data and analysis code for the paper **"Corpus Callosum Dysgenesis impairs metacognition: evidence from multi-modality and multi-cohort replications"**.

**Abstract:**
The corpus callosum is the largest commissure in the mammalian brain, supporting cognitive processes required for adapting to complex environments. Individuals with Corpus Callosum Dysgenesis (CCD) often exhibit deficits in social navigation and decision-making, yet metacognition—a key process supporting these functions—has not been comprehensively tested. 

Across three experiments (Online, Lab, and VR), we tested perceptual accuracy, confidence judgements, and metacognitive efficiency in individuals with CCD using variants of a Random Dot Kinematogram (RDK) task. We found that while individuals with CCD typically displayed normal perceptual accuracy, they failed to adjust confidence judgements in line with task difficulty. Computational modelling revealed this was driven by consistently lower metacognitive sensitivity (meta-d'). These results suggest the corpus callosum plays a crucial role in supporting metacognition.

## 🧪 Experiments

This study utilized three distinct experimental setups to test metacognition across different cohorts and modalities:

| Experiment | Modality | Participants | Task Variant | Key Features |
| :--- | :--- | :--- | :--- | :--- |
| **Exp 1** | **Online** | 10 CCD, 85 NT | Classic RDK | Hosted on Gorilla.sc; focused on baseline visual perception and confidence. |
| **Exp 2** | **Lab** | 11 CCD, 7 NT | Classic RDK | Conducted inside an MRI scanner (fMRI data reported separately). |
| **Exp 3** | **VR** | 10 CCD, 13 NT | Modified RDK | Meta Quest Pro; tested Binocular, Monocular, and Lateralized (single hemifield) presentations. |

## 📂 Repository Structure

The repository is organized to facilitate the reproduction of the analysis and figures presented in the manuscript.

```text
.
├── data/                   # Raw and processed behavioral data
├── analysis/               # R scripts for behavioral and computational analysis
│   ├── 01_preprocessing.R  # Data cleaning and exclusion criteria application
│   ├── 02_behavioral.R     # Accuracy, RT, and confidence analysis
│   ├── 03_modelling.R      # HMeta-d computational modelling (Hierarchical Bayesian)
│   └── 04_figures.R        # Code to generate manuscript figures
├── models/                 # Stan models for HMeta-d
└── figures/                # Output figures
