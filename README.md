# Statistical Analyses – Repeated Stressor Exposure in Wild Crickets

This repository contains the R script used to perform the statistical analyses for the manuscript:
**Gilford, E.R., Li, R., Rodríguez-Muñoz, R., Tregenza, T., & Kuijper, B.**, *"Trait-specific behavioural plasticity in response to repeated stressor exposure in wild crickets"*

## 📄 About this script

The analysis explores trait-specific behavioural plasticity in response to repeated simulated predator exposure. Crickets (*Gryllus campestris*) were exposed to vibrational stimuli of differing intensities (weak vs strong), and their responses across trials were measured.

The script performs:

1. **Dependencies and setup**
   - Loads R packages and sets working directory

2. **Data handling**
   - Reads in data file
   - Reverses post-emergence distance (to align directionality across traits)
   - Centres behavioural responses by burrow and stimulus type

3. **Statistical analysis**
   - Fits a multivariate Bayesian mixed-effects model (`rsmodel`) using `brms`
   - Extracts posterior predictions and individual plasticity estimates

4. **Visualisation**
   - Produces:
     - Posterior prediction plots (Figure 1)
     - Forest plot of fixed effects (Figure 2)
     - Paired individual slope plots across stimuli (Figure 3)

## 🧪 Study System

- Species: *Gryllus campestris*
- Location: Northern Spain
- Study period: March–May 2024
- Behavioural responses measured:
  - Emergence time/latency
  - Escape speed
  - Post-emergence distance

## 📂 Data

The data file (`data.csv`) may be provided upon reasonable request.

## 📝 Attribution

If you use or adapt this code, please cite the associated manuscript (currently in preparation):

**Gilford, E.R., Li, R., Rodríguez-Muñoz, R., Tregenza, T., & Kuijper, B.**  
*"Trait-specific behavioural plasticity in response to repeated stressor exposure in wild crickets."* In preparation, 2025.
