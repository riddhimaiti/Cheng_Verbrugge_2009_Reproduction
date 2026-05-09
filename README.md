# Diffusion-Induced Stresses in Spherical Electrode Particles

Numerical reproduction of **Cheng & Verbrugge (2009)** : *"Evolution of stress within a spherical insertion electrode particle under potentiostatic and galvanostatic operation"*, Journal of Power Sources **190**, 453–460.

---

## Overview

This project solves the solid-state diffusion equation inside a spherical insertion electrode particle using the **Finite Difference Method (FDM)** with a **Crank-Nicolson** scheme, then computes the resulting diffusion-induced **radial**, **tangential**, and **shear stresses**, as well as the **total elastic strain energy**.

Both **potentiostatic** (constant surface concentration) and **galvanostatic** (constant surface flux) boundary conditions are implemented, reproducing **Figures 1–4** from the original paper.

## Key Results

### Potentiostatic: Concentration & Stress Profiles
Concentration increases with time from the surface inward. Radial stress is tensile at the center and zero at the surface. Tangential stress is compressive at the surface (max at $t = 0$) and tensile at the center. Stresses peak transiently then decay.

<p align="center">
  <img src="fig1_potentiostatic.png" width="85%"/>
</p>

### Potentiostatic: Strain Energy vs Time
Total dimensionless strain energy peaks early then decays to zero as concentration equilibrates. Lower Poisson ratios produce higher strain energy.

<p align="center">
  <img src="fig2_strain_energy_potentiostatic.png" width="65%"/>
</p>

### Galvanostatic: Concentration & Stress Profiles
Under constant flux, stresses increase monotonically and approach a steady-state given by Eqs. 25–26 of the paper (shown as dashed lines). This is fundamentally different from the potentiostatic case.

<p align="center">
  <img src="fig3_galvanostatic.png" width="85%"/>
</p>

### Galvanostatic: Strain Energy vs Time
Unlike the potentiostatic case, strain energy increases monotonically to a finite steady-state value.

<p align="center">
  <img src="fig4_strain_energy_galvanostatic.png" width="65%"/>
</p>

## Verified Against Paper Predictions

| Paper Equation | Prediction | FDM Result |
|----------------|-----------|------------|
| Eq. 16: Max $\tilde{\sigma}_r$ at center | ≈ 0.4 at $\tau$ = 0.0574 | **0.3856 at $\tau$ = 0.0574** |
| Eq. 17: $\tilde{\sigma}_\theta$ at surface ($t \to 0$) | −1.0 | **−0.9894** |
| Eq. 25: Galvanostatic steady-state $\tilde{\sigma}_r(0)$ | 0.2 | **0.2001** |
| Eq. 26: Galvanostatic steady-state $\tilde{\sigma}_\theta(1)$ | −0.2 | **−0.1999** |

## Mathematical Formulation

### Governing Equation
1D spherical diffusion in dimensionless form ($x = r/R$, $\tau = Dt/R^2$):

$$\frac{\partial c}{\partial \tau} = \frac{1}{x^2}\frac{\partial}{\partial x}\left(x^2 \frac{\partial c}{\partial x}\right) = \frac{\partial^2 c}{\partial x^2} + \frac{2}{x}\frac{\partial c}{\partial x}$$

### Boundary Conditions

| Condition | Center ($x=0$) | Surface ($x=1$) |
|-----------|----------------|-----------------|
| **Potentiostatic** | $\partial c/\partial x = 0$ — **Neumann** (symmetry) | $c = 1$ — **Dirichlet** (fixed concentration) |
| **Galvanostatic** | $\partial c/\partial x = 0$ — **Neumann** (symmetry) | $\partial c/\partial x = 1$ — **Neumann** (fixed flux) |

> **Dirichlet BC** prescribes the *value* of the unknown ($c = \text{const}$).
> Used at the surface for potentiostatic control, where the electrode voltage fixes the surface concentration.
>
> **Neumann BC** prescribes the *derivative* (flux) of the unknown ($\partial c/\partial x = \text{const}$).
> Used at the center (zero flux by symmetry) and at the surface for galvanostatic control, where the applied current fixes the surface flux via $D\,\partial C/\partial r|_{r=R} = I/F$.

### Stress Formulas (Paper Eq. 3)
Normalized by $E\Omega\Delta c / [3(1-\nu)]$:

$$\tilde{\sigma}_r = \frac{2}{3}[\bar{c}(1) - \bar{c}(x)], \quad \tilde{\sigma}_\theta = \frac{1}{3}[2\bar{c}(1) + \bar{c}(x) - 3c(x)]$$

where $\bar{c}(x) = (3/x^3)\int_0^x x'^2 c(x') dx'$ is the volume-averaged concentration.

### Numerical Method
- **Crank-Nicolson** time discretization (unconditionally stable, 2nd-order)
- **Central finite differences** in space (2nd-order)
- **L'Hôpital's rule** at $x = 0$ to handle the coordinate singularity
- **Ghost-point method** for Neumann BC (galvanostatic surface flux)

## Requirements

- Python 3.8+
- NumPy
- Matplotlib
- Jupyter (to run the notebook interactively)

```bash
pip install numpy matplotlib jupyter
```

## References

1. **Y.-T. Cheng, M.W. Verbrugge**, "Evolution of stress within a spherical insertion electrode particle under potentiostatic and galvanostatic operation," *J. Power Sources* **190** (2009) 453–460: https://doi.org/10.1016/j.jpowsour.2009.01.021

2. Class Notes on **ECEG-6201 Analytical & Computational Methods by Murad Ridwan** from Addis Ababa University: http://ndl.ethernet.edu.et/bitstream/123456789/87841/6/chapter4.pdf
3. **Crank-Nikolson method (Wikipedia)**: https://en.wikipedia.org/wiki/Crank-Nicolson_method

## License

This project is for educational and research purposes. The original paper is by Y.-T. Cheng and M.W. Verbrugge (General Motors R&D Center).
