# Post-activation-Performance-Enhancement-of-Ballistic-Military-Press-Simulation-Data-Analysis

# Power Analysis for Post-Activation Performance Enhancement in Military Press 

## 📌 Project Overview
This project uses simulation-based power analysis to estimate whether a typical sports science study design is sufficiently powered to detect Post-Activation Performance Enhancement (PAPE) effects in upper-body performance.

The analysis focuses on the seated military press (SMP) as a conditioning activity and its impact on ballistic performance (SBMPT).

---

## ❓ Research Question
Is a sample size of N = 30 sufficient to reliably detect a small-to-moderate PAPE effect (Cohen's d ≈ 0.3) in military press performance?

---

## 📚 Background
Post-Activation Performance Enhancement (PAPE) refers to acute improvements in muscular performance following high-intensity contractions.

Previous literature suggests:
- Typical effect sizes around **d ≈ 0.3**
- Most studies focus on bench press → bench press throw
- Limited evidence exists for vertical pushing movements (e.g., military press)

---

## 🧪 Study Design (Simulated)
- Planned sample size: **N = 30 participants**
- Within-subject design:
  - Baseline vs PAPE condition
  - Multiple time points: 0, 5, 7, 10 minutes

### Measured outcomes:
- Peak Power Output (PPO)
- Peak Velocity Output (PVO)

---

## 📊 Simulation Approach
- Monte Carlo simulation (10,000 iterations)
- Effect size assumption: **d ≈ 0.3**
- Data generated based on:
  - pilot measurements (single participant)
  - literature-based estimates

### Statistical model:
- Linear mixed-effects models (LMM)
- Fixed effects:
  - condition (baseline vs PAPE)
  - time point
  - interaction (condition × time)
- Random effects:
  - participant ID
  - repeated attempts

---

## 📈 Key Results

- Estimated statistical power for detecting PAPE effect (d ≈ 0.3):
  - **Power ≈ 0.93 (93%) with N = 30**

- Power estimation based on:
  - proportion of simulations where 95% CI excluded zero

---

## 🧠 Interpretation

The simulation suggests that:

- A sample size of **N = 30 is sufficient** to detect small-to-moderate PAPE effects under the assumed conditions
- Properly designed within-subject studies can achieve high statistical power even with moderate sample sizes

This has important implications for:
- experimental design in sports science
- improving reproducibility of PAPE research
- avoiding underpowered studies common in the literature

---

## 🧰 Tools
R, Monte Carlo simulation, linear mixed-effects models (lme4)


See the preregistration of the project: https://osf.io/cu27d/overview
