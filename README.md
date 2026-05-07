# Complete Mathematical Derivation: FDM for Spherical Diffusion-Induced Stresses

## Reproducing Cheng & Verbrugge (2009), *J. Power Sources* 190, 453–460

---

## 1. Problem Statement

We solve for the concentration profile $C(r,t)$ inside a sphere of radius $R$ undergoing solid-state diffusion, then compute the resulting mechanical stresses.

### 1.1 Governing PDE — Fick's Second Law in Spherical Coordinates

$$\frac{\partial C}{\partial t}=\frac{D}{r^2}\frac{\partial}{\partial r}\left(r^2 \frac{\partial C}{\partial r}\right) \tag{1}$$

where $D$ is the diffusion coefficient (constant), $r$ is the radial coordinate, and $C(r,t)$ is the solute concentration.

**Expanding the spatial operator:**

$$\frac{1}{r^2}\frac{\partial}{\partial r}\left(r^2 \frac{\partial C}{\partial r}\right) = \frac{\partial^2 C}{\partial r^2} + \frac{2}{r}\frac{\partial C}{\partial r} \tag{2}$$

*Proof:* Apply the product rule to $\frac{\partial}{\partial r}(r^2 \frac{\partial C}{\partial r}) = 2r\frac{\partial C}{\partial r} + r^2\frac{\partial^2 C}{\partial r^2}$, then divide by $r^2$.

### 1.2 Dimensionless Variables

Define:
- $x = r/R$ (normalized radius, $0 \le x \le 1$)
- $\tau = Dt/R^2$ (dimensionless time)

**Potentiostatic:** $c = (C - C_0)/(C_R - C_0)$

**Galvanostatic:** $c = (C - C_0)/(IR/FD)$

where $C_0$ is the initial concentration, $C_R$ is the surface concentration (potentiostatic), $I$ is the current density, and $F$ is Faraday's constant.

The PDE becomes:

$$\frac{\partial c}{\partial \tau} = \frac{\partial^2 c}{\partial x^2} + \frac{2}{x}\frac{\partial c}{\partial x} \tag{3}$$

### 1.3 Boundary Conditions

**At the center ($x = 0$):** Symmetry requires $\frac{\partial c}{\partial x}\big|_{x=0} = 0$

**At the surface ($x = 1$):**
- Potentiostatic: $c(1, \tau) = 1$
- Galvanostatic: $\frac{\partial c}{\partial x}\big|_{x=1} = 1$

**Initial condition:** $c(x, 0) = 0$

### 1.4 Singularity at $x = 0$

The term $\frac{2}{x}\frac{\partial c}{\partial x}$ is singular at $x = 0$. We resolve this using **L'Hôpital's rule**:

$$\lim_{x \to 0} \frac{2}{x}\frac{\partial c}{\partial x} = 2\lim_{x \to 0}\frac{\partial^2 c}{\partial x^2} = 2\frac{\partial^2 c}{\partial x^2}\bigg|_{x=0}$$

Therefore, at $x = 0$:

$$\frac{\partial c}{\partial \tau}\bigg|_{x=0} = 3\frac{\partial^2 c}{\partial x^2}\bigg|_{x=0} \tag{4}$$

---

## 2. Finite Difference Discretization

### 2.1 Grid Setup

Divide $[0, 1]$ into $N$ equal intervals:
- $x_i = i \cdot \Delta x$ for $i = 0, 1, \ldots, N$
- $\Delta x = 1/N$
- $c_i^n \approx c(x_i, n\cdot\Delta\tau)$

### 2.2 Standard Finite Difference Approximations

**Second derivative (central difference):**

$$\frac{\partial^2 c}{\partial x^2}\bigg|_{x_i} \approx \frac{c_{i+1} - 2c_i + c_{i-1}}{\Delta x^2} \tag{5}$$

**First derivative (central difference):**

$$\frac{\partial c}{\partial x}\bigg|_{x_i} \approx \frac{c_{i+1} - c_{i-1}}{2\Delta x} \tag{6}$$

Both are second-order accurate: $O(\Delta x^2)$.

### 2.3 Spatial Operator at Interior Points ($i = 1, \ldots, N-1$)

Substituting (5) and (6) into (3):

$$\mathcal{L}(c)_i = \frac{c_{i+1} - 2c_i + c_{i-1}}{\Delta x^2} + \frac{2}{x_i}\cdot\frac{c_{i+1} - c_{i-1}}{2\Delta x}$$

Since $x_i = i\cdot\Delta x$:

$$\frac{2}{x_i \cdot 2\Delta x} = \frac{1}{i\cdot\Delta x^2}$$

Therefore:

$$\mathcal{L}(c)_i = \frac{1}{\Delta x^2}\left[\frac{i-1}{i}\,c_{i-1} - 2c_i + \frac{i+1}{i}\,c_{i+1}\right] \tag{7}$$

### 2.4 Spatial Operator at the Center ($i = 0$)

From Eq. (4): $\frac{\partial c}{\partial \tau}\big|_{x=0} = 3\frac{\partial^2 c}{\partial x^2}\big|_{x=0}$

Using the symmetry condition $c_{-1} = c_1$ (ghost point):

$$\frac{\partial^2 c}{\partial x^2}\bigg|_{x=0} \approx \frac{c_1 - 2c_0 + c_{-1}}{\Delta x^2} = \frac{2(c_1 - c_0)}{\Delta x^2} \tag{8}$$

Therefore:

$$\mathcal{L}(c)_0 = \frac{6(c_1 - c_0)}{\Delta x^2} \tag{9}$$

### 2.5 The Crank-Nicolson Scheme

The Crank-Nicolson method averages the spatial operator at times $n$ and $n+1$:

$$\frac{c_i^{n+1} - c_i^n}{\Delta\tau} = \frac{1}{2}\left[\mathcal{L}(c^{n+1})_i + \mathcal{L}(c^n)_i\right] \tag{10}$$

**Properties:**
- **Unconditionally stable** (no restriction on $\Delta\tau/\Delta x^2$)
- **Second-order accurate** in both space and time: $O(\Delta x^2 + \Delta\tau^2)$
- **Implicit** — requires solving a linear system at each time step

Define $s = \Delta\tau / (2\Delta x^2)$.

#### 2.5.1 Interior Points ($i = 1, \ldots, N-1$)

Substituting (7) into (10) and rearranging:

**Left-hand side (unknowns at $n+1$):**

$$-s\frac{i-1}{i}\,c_{i-1}^{n+1} + (1 + 2s)\,c_i^{n+1} - s\frac{i+1}{i}\,c_{i+1}^{n+1}$$

**Right-hand side (knowns at $n$):**

$$s\frac{i-1}{i}\,c_{i-1}^n + (1 - 2s)\,c_i^n + s\frac{i+1}{i}\,c_{i+1}^n$$

**Derivation step-by-step:**

Starting from Eq. (10):

$$c_i^{n+1} - c_i^n = s\left[\frac{i-1}{i}(c_{i-1}^{n+1} + c_{i-1}^n) - 2(c_i^{n+1} + c_i^n) + \frac{i+1}{i}(c_{i+1}^{n+1} + c_{i+1}^n)\right]$$

Move all $n+1$ terms to the left:

$$c_i^{n+1} - s\frac{i-1}{i}c_{i-1}^{n+1} + 2s\,c_i^{n+1} - s\frac{i+1}{i}c_{i+1}^{n+1} = c_i^n + s\frac{i-1}{i}c_{i-1}^n - 2s\,c_i^n + s\frac{i+1}{i}c_{i+1}^n$$

This gives the tridiagonal system $\mathbf{A}\,\mathbf{c}^{n+1} = \mathbf{B}\,\mathbf{c}^n + \mathbf{b}$.

#### 2.5.2 Center Point ($i = 0$)

From Eq. (9):

$$(1 + 6s)\,c_0^{n+1} - 6s\,c_1^{n+1} = (1 - 6s)\,c_0^n + 6s\,c_1^n \tag{11}$$

#### 2.5.3 Summary of Matrix Coefficients

| Location | Sub-diagonal | Diagonal | Super-diagonal |
|----------|-------------|----------|---------------|
| $i = 0$ | — | $1 + 6s$ | $-6s$ |
| $i = 1\ldots N-1$ | $-s(i-1)/i$ | $1 + 2s$ | $-s(i+1)/i$ |

The **B matrix** (RHS) has the same structure with opposite signs on $s$:

| Location | Sub-diagonal | Diagonal | Super-diagonal |
|----------|-------------|----------|---------------|
| $i = 0$ | — | $1 - 6s$ | $6s$ |
| $i = 1\ldots N-1$ | $s(i-1)/i$ | $1 - 2s$ | $s(i+1)/i$ |

---

## 3. Boundary Condition Implementation

### 3.1 Potentiostatic ($c(1, \tau) = 1$)

Simply set:
- $A_{N,N} = 1$, all other entries in row $N$ are zero
- RHS: $b_N = 1$

### 3.2 Galvanostatic ($\partial c/\partial x|_{x=1} = 1$)

**Ghost-point approach:** Introduce a fictitious point $c_{N+1}$ beyond the boundary:

$$\frac{\partial c}{\partial x}\bigg|_{x=1} = \frac{c_{N+1} - c_{N-1}}{2\Delta x} = 1 \quad \Rightarrow \quad c_{N+1} = c_{N-1} + 2\Delta x \tag{12}$$

Apply the **interior stencil at $i = N$** and substitute $c_{N+1}$. Let $a_\text{sub} = -s(N-1)/N$ and $a_\text{sup} = -s(N+1)/N$.

After substitution:
- Coefficient of $c_{N-1}$: $a_\text{sub} + a_\text{sup} = -2s$
- Diagonal: $1 + 2s$ (unchanged)
- Constant: $-a_\text{sup} \cdot 2\Delta x$ moves to the RHS

For the B matrix similarly with $b_\text{sub} = s(N-1)/N$, $b_\text{sup} = s(N+1)/N$:

$$\text{RHS constant} = (-a_\text{sup} + b_\text{sup}) \cdot 2\Delta x = \frac{2s(N+1)}{N} \cdot 2\Delta x = \frac{2(N+1)\Delta\tau}{N\Delta x}$$

---

## 4. Stress Calculations

### 4.1 Volume-Averaged Concentration

$$\bar{c}(x) = \frac{3}{x^3}\int_0^x {x'}^2\,c(x')\,dx' \tag{14}$$

At $x = 0$: $\bar{c}(0) = c(0)$ (by L'Hôpital). Evaluated numerically via the trapezoidal rule.

### 4.2 Dimensionless Stresses (from Eq. 3 of paper)

Normalized by $E\Omega\Delta c / [3(1-\nu)]$:

$$\tilde{\sigma}_r = \frac{2}{3}[\bar{c}(1) - \bar{c}(x)], \quad \tilde{\sigma}_\theta = \frac{1}{3}[2\bar{c}(1) + \bar{c}(x) - 3c(x)] \tag{15}$$

$$\tilde{\sigma}_\text{shear} = \frac{\tilde{\sigma}_r - \tilde{\sigma}_\theta}{2} \tag{16}$$

*Derivation of (15):* From paper Eq. (3): $\sigma_r = \frac{2E\Omega}{9(1-\nu)}[C_\text{av}(R) - C_\text{av}(r)]$. Dividing by $\frac{E\Omega\Delta c}{3(1-\nu)}$ gives:

$$\tilde{\sigma}_r = \frac{2}{9} \cdot 3 \cdot \frac{C_\text{av}(R)-C_\text{av}(r)}{\Delta c} = \frac{2}{3}[\bar{c}(1)-\bar{c}(x)]$$

---

## 5. Strain Energy

### 5.1 Strain Energy Density (Eq. 6 of paper)

$$e(r) = \frac{\sigma_r^2 + 2\sigma_\theta^2 - 2\nu\sigma_\theta(2\sigma_r + \sigma_\theta)}{2E} \tag{17}$$

### 5.2 Dimensionless Total Strain Energy (Eqs. 19, 28)

$$\tilde{E}_T = \frac{E_T}{2\pi R^3 E[\Omega\Delta c/(3(1-\nu))]^2} = \int_0^1 [\tilde{\sigma}_r^2 + 2\tilde{\sigma}_\theta^2 - 2\nu\tilde{\sigma}_\theta(2\tilde{\sigma}_r + \tilde{\sigma}_\theta)]x^2\,dx \tag{18}$$

---

## 6. Analytical Solutions for Validation

### 6.1 Potentiostatic (Eq. 11 — Carslaw & Jaeger [32])

$$c(x, \tau) = 1 + 2\sum_{n=1}^{\infty}\frac{(-1)^n}{n\pi x}\sin(n\pi x)\,e^{-n^2\pi^2\tau} \tag{19}$$

### 6.2 Galvanostatic (Eq. 20)

$$c(x, \tau) = 3\tau + \frac{x^2}{2} - \frac{3}{10} - \frac{2}{x}\sum_{n=1}^{\infty}\frac{\sin(\alpha_n x)}{\alpha_n^2\sin(\alpha_n)}e^{-\alpha_n^2\tau} \tag{20}$$

where $\alpha_n$ are positive roots of $\tan(\alpha) = \alpha$.

---

## 7. Implementation Summary

| Aspect | Choice | Rationale |
|--------|--------|-----------|
| Time scheme | Crank-Nicolson | Unconditionally stable, 2nd-order |
| Spatial scheme | Central differences | 2nd-order, symmetric |
| Center singularity | L'Hôpital + ghost point | Avoids division by zero |
| Neumann BC | Ghost point into interior stencil | Maintains 2nd-order accuracy |
| Linear solve | Dense `numpy.linalg.solve` | Simple; tridiagonal solver possible |
| Integration | Trapezoidal rule | Sufficient for smooth profiles |
| Grid | $N = 200$ | Errors $\sim O(10^{-5})$ vs analytical |

## References

1. Y.-T. Cheng, M.W. Verbrugge, *J. Power Sources* 190 (2009) 453–460
2. H.S. Carslaw, J.C. Jaeger, *Conduction of Heat in Solids*, 2nd ed. (1959)
3. S.P. Timoshenko, J.N. Goodier, *Theory of Elasticity*, 3rd ed. (1970)
4. R.J. LeVeque, *Finite Difference Methods for ODEs and PDEs* (2007)
5. J. Crank, *The Mathematics of Diffusion*, 2nd ed. (1975)
