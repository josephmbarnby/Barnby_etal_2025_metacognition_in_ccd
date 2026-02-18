# Corpus Callosum Dysgenesis impairs metacognition: evidence from multi-modality and multi-cohort replications

**Authors:** Joseph M. Barnby*¹², Ryan Dean*³, H. Burgess³, Peter Dayan ⌇ ⁴⁵, Linda J. Richards ⌇ ³

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

## 🛠️ Usage & Reproducibility

### Prerequisites
[cite_start]The analysis was conducted using **R (v4.3.3)**[cite: 101]. Key dependencies include:
* `tidyverse`
* `rstan` (for HMeta-d modelling)
* `brms` / `lme4` (for mixed-effects models)
* `ggplot2` (for visualization)

### Running the Analysis

Clone this repository:
   ```bash
   git clone [https://github.com/josephmbarnby/Barnby_etal_2025_metacognition_in_ccd.git](https://github.com/josephmbarnby/Barnby_etal_2025_metacognition_in_ccd.git)
  ```

### Citation
If you use this data or code, please cite our paper:

Barnby, J.M., Dean, R., Burgess, H., Dayan, P., & Richards, L.J. (2026). Corpus Callosum Dysgenesis impairs metacognition: evidence from multi-modality and multi-cohort replications. _Nature Commmunications_

Barnby, J.M., Dean, R., Burgess, H., Dayan, P., & Richards, L.J. (2025). Corpus Callosum Dysgenesis impairs metacognition: evidence from multi-modality and multi-cohort replications. [Preprint/Journal Reference].
