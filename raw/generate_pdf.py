"""Generate a PDF derivation document using matplotlib's LaTeX rendering."""
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import textwrap

plt.rcParams.update({
    'text.usetex': False,
    'mathtext.fontset': 'cm',
    'font.family': 'serif',
    'font.size': 11,
})

def add_page(pdf, texts, title=None):
    """Add a page with text blocks. Each text is (y_pos, content, fontsize, style)."""
    fig = plt.figure(figsize=(8.5, 11))
    ax = fig.add_axes([0.08, 0.05, 0.84, 0.90])
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    ax.axis('off')
    if title:
        ax.text(0.5, 0.97, title, fontsize=14, fontweight='bold',
                ha='center', va='top', transform=ax.transAxes)
    for y, txt, fs, style in texts:
        fw = 'bold' if style == 'bold' else 'normal'
        fst = 'italic' if style == 'italic' else 'normal'
        ax.text(0.02, y, txt, fontsize=fs, fontweight=fw, fontstyle=fst,
                va='top', transform=ax.transAxes, wrap=True,
                family='serif')
    pdf.savefig(fig)
    plt.close(fig)

def eq_page(pdf, texts):
    """Page with mixed text and equations."""
    fig = plt.figure(figsize=(8.5, 11))
    ax = fig.add_axes([0.06, 0.04, 0.88, 0.92])
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    ax.axis('off')
    for y, txt, fs, ha, style in texts:
        fw = 'bold' if 'bold' in style else 'normal'
        fst = 'italic' if 'italic' in style else 'normal'
        ax.text(0.5 if ha=='center' else 0.02, y, txt, fontsize=fs,
                fontweight=fw, fontstyle=fst, ha=ha, va='top',
                transform=ax.transAxes, family='serif')
    pdf.savefig(fig)
    plt.close(fig)

outpath = 'FDM_Derivation_Document.pdf'
with PdfPages(outpath) as pdf:

    # ---- PAGE 1: Title + Problem Statement ----
    eq_page(pdf, [
        (0.96, 'Complete Mathematical Derivation', 16, 'center', 'bold'),
        (0.92, 'FDM for Spherical Diffusion-Induced Stresses', 13, 'center', 'bold'),
        (0.88, 'Reproducing Cheng & Verbrugge (2009), J. Power Sources 190, 453–460', 10, 'center', 'italic'),
        (0.83, '1.  Problem Statement', 13, 'left', 'bold'),
        (0.80, '1.1  Governing PDE — Fick\'s Second Law in Spherical Coordinates', 11, 'left', 'bold'),
        (0.77, 'We solve for concentration C(r,t) inside a sphere of radius R with constant diffusivity D:', 10, 'left', ''),
        (0.73, r'$\frac{\partial C}{\partial t} = \frac{D}{r^2}\frac{\partial}{\partial r}\left(r^2\frac{\partial C}{\partial r}\right)$          ......(1)', 12, 'center', ''),
        (0.68, 'Expanding the spatial operator via the product rule:', 10, 'left', ''),
        (0.64, r'$\frac{1}{r^2}\frac{\partial}{\partial r}\left(r^2\frac{\partial C}{\partial r}\right) = \frac{\partial^2 C}{\partial r^2} + \frac{2}{r}\frac{\partial C}{\partial r}$          ......(2)', 12, 'center', ''),
        (0.59, '1.2  Dimensionless Variables', 11, 'left', 'bold'),
        (0.56, r'Define:  $x = r/R$  (normalized radius),   $\tau = Dt/R^2$  (dimensionless time)', 10, 'left', ''),
        (0.53, r'Potentiostatic:  $c = (C - C_0)/(C_R - C_0)$', 10, 'left', ''),
        (0.50, r'Galvanostatic:   $c = (C - C_0)/(IR/FD)$', 10, 'left', ''),
        (0.47, 'The dimensionless PDE becomes:', 10, 'left', ''),
        (0.43, r'$\frac{\partial c}{\partial \tau} = \frac{\partial^2 c}{\partial x^2} + \frac{2}{x}\frac{\partial c}{\partial x}$          ......(3)', 12, 'center', ''),
        (0.38, '1.3  Boundary Conditions', 11, 'left', 'bold'),
        (0.35, r'Center ($x = 0$):   $\partial c / \partial x = 0$   (symmetry)', 10, 'left', ''),
        (0.32, r'Surface ($x = 1$):  Potentiostatic: $c(1,\tau)=1$  |  Galvanostatic: $\partial c/\partial x|_{x=1}=1$', 10, 'left', ''),
        (0.29, r'Initial condition:  $c(x,0) = 0$', 10, 'left', ''),
        (0.24, '1.4  Singularity at x = 0 (L\'Hôpital\'s Rule)', 11, 'left', 'bold'),
        (0.21, r'The term $2x^{-1}\partial c/\partial x$ is singular at $x=0$. By L\'Hôpital:', 10, 'left', ''),
        (0.17, r'$\lim_{x \to 0}  \frac{2}{x}\frac{\partial c}{\partial x} = 2\frac{\partial^2 c}{\partial x^2}|_{x=0}$', 12, 'center', ''),
        (0.12, 'Therefore at x = 0 the PDE becomes:', 10, 'left', ''),
        (0.08, r'$\frac{\partial c}{\partial \tau}|_{x=0} = 3\,\frac{\partial^2 c}{\partial x^2}|_{x=0}$          ......(4)', 12, 'center', ''),
    ])

    # ---- PAGE 2: FDM Discretization ----
    eq_page(pdf, [
        (0.96, '2.  Finite Difference Discretization', 13, 'left', 'bold'),
        (0.93, '2.1  Grid Setup', 11, 'left', 'bold'),
        (0.90, r'Divide $[0,1]$ into $N$ equal intervals:  $x_i = i\cdot\Delta x$,  $\Delta x = 1/N$,  $i = 0,1,...,N$', 10, 'left', ''),
        (0.86, '2.2  Standard Finite Difference Approximations', 11, 'left', 'bold'),
        (0.83, 'Second derivative (central, 2nd-order accurate):', 10, 'left', ''),
        (0.79, r'$\frac{\partial^2 c}{\partial x^2}|_{x_i} \approx \frac{c_{i+1} - 2c_i + c_{i-1}}{\Delta x^2}$          ......(5)', 12, 'center', ''),
        (0.74, 'First derivative (central, 2nd-order accurate):', 10, 'left', ''),
        (0.70, r'$\frac{\partial c}{\partial x}|_{x_i} \approx \frac{c_{i+1} - c_{i-1}}{2\Delta x}$          ......(6)', 12, 'center', ''),
        (0.65, '2.3  Spatial Operator at Interior Points (i = 1, ..., N−1)', 11, 'left', 'bold'),
        (0.62, r'Substituting (5) and (6) into (3), and using $x_i = i\cdot\Delta x$:', 10, 'left', ''),
        (0.58, r'$\mathcal{L}(c)_i = \frac{1}{\Delta x^2}\left[\frac{i-1}{i}\,c_{i-1}\; -\; 2c_i\; +\; \frac{i+1}{i}\,c_{i+1}\right]$          ......(7)', 12, 'center', ''),
        (0.53, r'Note: the factors $(i\pm 1)/i$ arise from combining the $1/\Delta x^2$ and $1/(i\Delta x^2)$ terms.', 10, 'left', 'italic'),
        (0.48, '2.4  Spatial Operator at the Center (i = 0)', 11, 'left', 'bold'),
        (0.45, r'From Eq.(4), using symmetry ghost point $c_{-1} = c_1$:', 10, 'left', ''),
        (0.41, r'$\frac{\partial^2 c}{\partial x^2}|_{x=0} \approx \frac{2(c_1 - c_0)}{\Delta x^2}$   $\Rightarrow$   $\mathcal{L}(c)_0 = \frac{6(c_1 - c_0)}{\Delta x^2}$          ......(8)', 12, 'center', ''),
        (0.35, '2.5  Crank-Nicolson Time Discretization', 11, 'left', 'bold'),
        (0.32, 'Average the spatial operator at time levels n and n+1 (unconditionally stable, 2nd-order):', 10, 'left', ''),
        (0.28, r'$\frac{c_i^{n+1} - c_i^{n}}{\Delta\tau} = \frac{1}{2}\left[\mathcal{L}(c^{n+1})_i + \mathcal{L}(c^{n})_i\right]$          ......(9)', 12, 'center', ''),
        (0.23, r'Define $s = \Delta\tau\,/\,(2\Delta x^2)$.  Rearranging gives:  $\mathbf{A}\,\mathbf{c}^{n+1} = \mathbf{B}\,\mathbf{c}^{n} + \mathbf{b}$', 10, 'left', ''),
        (0.18, 'where A and B are tridiagonal matrices (see next page for coefficients).', 10, 'left', ''),
    ])

    # ---- PAGE 3: Matrix Coefficients ----
    eq_page(pdf, [
        (0.96, '2.6  Tridiagonal Matrix Coefficients', 13, 'left', 'bold'),
        (0.93, 'Interior points (i = 1, ..., N−1) — A matrix (LHS, implicit):', 10, 'left', 'bold'),
        (0.90, r'Sub-diagonal:  $-s(i-1)/i$', 10, 'left', ''),
        (0.87, r'Diagonal:          $1 + 2s$', 10, 'left', ''),
        (0.84, r'Super-diagonal: $-s(i+1)/i$', 10, 'left', ''),
        (0.80, 'B matrix (RHS, explicit) — same structure, opposite sign on s:', 10, 'left', 'bold'),
        (0.77, r'Sub-diagonal:  $+s(i-1)/i$', 10, 'left', ''),
        (0.74, r'Diagonal:          $1 - 2s$', 10, 'left', ''),
        (0.71, r'Super-diagonal: $+s(i+1)/i$', 10, 'left', ''),
        (0.66, 'Center point (i = 0):', 10, 'left', 'bold'),
        (0.62, r'$(1+6s)\,c_0^{n+1} - 6s\,c_1^{n+1} = (1-6s)\,c_0^{n} + 6s\,c_1^{n}$          ......(10)', 12, 'center', ''),
        (0.56, '3.  Boundary Condition Implementation', 13, 'left', 'bold'),
        (0.53, '3.1  Potentiostatic:  c(1, τ) = 1', 11, 'left', 'bold'),
        (0.50, r'Set $A_{N,N}=1$, all other entries in row N zero.  RHS$_N = 1$.', 10, 'left', ''),
        (0.46, r'3.2  Galvanostatic:  $\partial c/\partial x|_{x=1} = 1$  (Ghost-Point Method)', 11, 'left', 'bold'),
        (0.43, r'Introduce ghost point $c_{N+1}$ beyond the boundary:', 10, 'left', ''),
        (0.39, r'$\frac{c_{N+1} - c_{N-1}}{2\Delta x} = 1$   $\Rightarrow$   $c_{N+1} = c_{N-1} + 2\Delta x$          ......(11)', 12, 'center', ''),
        (0.34, 'Substitute (11) into the interior stencil (7) at i = N:', 10, 'left', ''),
        (0.31, r'Let $a_{sub} = -s(N\!-\!1)/N$,  $a_{sup} = -s(N\!+\!1)/N$.  After substitution:', 10, 'left', ''),
        (0.28, r'$A_{N,N-1} = a_{sub}+a_{sup} = -2s$,     $A_{N,N} = 1+2s$', 10, 'left', ''),
        (0.25, r'RHS constant $= (-a_{sup}+b_{sup})\cdot 2\Delta x = \frac{2(N+1)\Delta\tau}{N\,\Delta x}$', 10, 'left', ''),
        (0.20, 'This preserves 2nd-order accuracy at the boundary.', 10, 'left', 'italic'),
    ])

    # ---- PAGE 4: Stress Calculations ----
    eq_page(pdf, [
        (0.96, '4.  Stress Calculations', 13, 'left', 'bold'),
        (0.93, '4.1  Volume-Averaged Concentration', 11, 'left', 'bold'),
        (0.89, r'$\bar{c}(x) = \frac{3}{x^3}\int_0^x {x\,}^{\prime\,2}\,c(x\,^\prime)\,dx\,^\prime$          ......(12)', 12, 'center', ''),
        (0.84, r'At $x=0$:  $\bar{c}(0) = c(0)$  (by L\'Hôpital).  Evaluated via trapezoidal rule.', 10, 'left', ''),
        (0.80, '4.2  Dimensionless Stresses  (from Paper Eq. 3)', 11, 'left', 'bold'),
        (0.77, r'Normalized by  $E\Omega\Delta c\,/\,[3(1-\nu)]$:', 10, 'left', ''),
        (0.73, r'Radial:       $\tilde{\sigma}_r = \frac{2}{3}\left[\bar{c}(1) - \bar{c}(x)\right]$          ......(13)', 12, 'center', ''),
        (0.68, r'Tangential:  $\tilde{\sigma}_\theta = \frac{1}{3}\left[2\bar{c}(1) + \bar{c}(x) - 3c(x)\right]$          ......(14)', 12, 'center', ''),
        (0.63, r'Shear:         $\tilde{\sigma}_{shear} = (\tilde{\sigma}_r - \tilde{\sigma}_\theta)\,/\,2$          ......(15)', 12, 'center', ''),
        (0.57, 'Derivation of (13):', 10, 'left', 'bold'),
        (0.54, r'From paper Eq.(3):  $\sigma_r = \frac{2E\Omega}{9(1-\nu)}[C_{av}(R) - C_{av}(r)]$', 10, 'left', ''),
        (0.50, r'Dividing by $E\Omega\Delta c/[3(1-\nu)]$:', 10, 'left', ''),
        (0.46, r'$\tilde{\sigma}_r = \frac{2}{9}\cdot\frac{3}{\Delta c}[C_{av}(R)-C_{av}(r)] = \frac{2}{3}[\bar{c}(1)-\bar{c}(x)]$   $(verified)$', 11, 'center', ''),
        (0.40, '4.3  Key Stress Properties', 11, 'left', 'bold'),
        (0.37, r'• At $x=0$:  $\bar{c}(0)=c(0)$  $\Rightarrow$  $\tilde{\sigma}_r = \tilde{\sigma}_\theta$  (hydrostatic)', 10, 'left', ''),
        (0.34, r'• At $x=1$:  $\tilde{\sigma}_r(1)=0$  always (free surface)', 10, 'left', ''),
        (0.31, r'• Potentiostatic max radial: $\tilde{\sigma}_{r,max}\approx 0.4$ at $\tau=0.0574$  (Eq. 16)', 10, 'left', ''),
        (0.28, r'• Potentiostatic max tangential: $\tilde{\sigma}_\theta(1,0) = -1$  (Eq. 17)', 10, 'left', ''),
        (0.25, r'• Galvanostatic steady-state: $\tilde{\sigma}_r = \frac{1}{5}(1-x^2)$  (Eq. 25)', 10, 'left', ''),
        (0.22, r'• Galvanostatic steady-state: $\tilde{\sigma}_\theta = \frac{1}{5}(1-2x^2)$  (Eq. 26)', 10, 'left', ''),
    ])

    # ---- PAGE 5: Strain Energy + Analytical ----
    eq_page(pdf, [
        (0.96, '5.  Strain Energy', 13, 'left', 'bold'),
        (0.93, '5.1  Strain Energy Density  (Paper Eq. 6)', 11, 'left', 'bold'),
        (0.90, 'For an isotropic sphere with spherical symmetry:', 10, 'left', ''),
        (0.86, r'$e(r) = \frac{1}{2E}\left[\sigma_r^2 + 2\sigma_\theta^2 - 2\nu\sigma_\theta(2\sigma_r+\sigma_\theta)\right]$          ......(16)', 12, 'center', ''),
        (0.81, '5.2  Dimensionless Total Strain Energy  (Paper Eqs. 19, 28)', 11, 'left', 'bold'),
        (0.77, r'$\tilde{E}_T = \int_0^1\left[\tilde{\sigma}_r^2 + 2\tilde{\sigma}_\theta^2 - 2\nu\tilde{\sigma}_\theta(2\tilde{\sigma}_r+\tilde{\sigma}_\theta)\right]x^2\,dx$          ......(17)', 12, 'center', ''),
        (0.72, r'Normalized by  $2\pi R^3 E\,[\Omega\Delta c\,/\,(3(1-\nu))]^2$.  Evaluated via trapezoidal rule.', 10, 'left', ''),
        (0.66, '6.  Analytical Solutions for Validation', 13, 'left', 'bold'),
        (0.63, '6.1  Potentiostatic  (Paper Eq. 11, Carslaw & Jaeger)', 11, 'left', 'bold'),
        (0.59, r'$c(x,\tau) = 1 + 2\sum_{n=1}^{\infty}\frac{(-1)^n}{n\pi x}\sin(n\pi x)\,e^{-n^2\pi^2\tau}$          ......(18)', 12, 'center', ''),
        (0.54, r'At $x=0$: use $\lim_{x\to 0}\sin(n\pi x)/(n\pi x)=1$', 10, 'left', ''),
        (0.50, '6.2  Galvanostatic  (Paper Eq. 20)', 11, 'left', 'bold'),
        (0.46, r'$c(x,\tau) = 3\tau + \frac{x^2}{2} - \frac{3}{10} - \frac{2}{x}\sum_{n=1}^{\infty}\frac{\sin(\alpha_n x)}{\alpha_n^2\sin(\alpha_n)}e^{-\alpha_n^2\tau}$          ......(19)', 11, 'center', ''),
        (0.41, r'where $\alpha_n$ are positive roots of  $\tan(\alpha)=\alpha$.', 10, 'left', ''),
        (0.35, '7.  Implementation Summary', 13, 'left', 'bold'),
        (0.32, 'Time scheme:           Crank-Nicolson (unconditionally stable, 2nd-order)', 10, 'left', ''),
        (0.29, 'Spatial scheme:        Central differences (2nd-order, symmetric)', 10, 'left', ''),
        (0.26, "Center singularity:   L'Hôpital + symmetry ghost point", 10, 'left', ''),
        (0.23, 'Neumann BC:           Ghost point substituted into interior stencil', 10, 'left', ''),
        (0.20, 'Linear solve:           numpy.linalg.solve (dense)', 10, 'left', ''),
        (0.17, 'Integration:             Trapezoidal rule', 10, 'left', ''),
        (0.14, r'Grid:                       $N = 200$  (errors $\sim 10^{-5}$ vs analytical)', 10, 'left', ''),
        (0.09, 'References:', 10, 'left', 'bold'),
        (0.06, '1. Cheng & Verbrugge, JPS 190 (2009)  2. Carslaw & Jaeger (1959)  3. Timoshenko (1970)', 9, 'left', ''),
        (0.03, '4. LeVeque, FDM for ODEs & PDEs (2007)  5. Crank, Mathematics of Diffusion (1975)', 9, 'left', ''),
    ])

print(f"PDF created: {outpath}")
