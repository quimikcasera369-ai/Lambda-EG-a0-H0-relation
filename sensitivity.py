"""
Sensitivity analysis of the fixed point Phi = 0 with respect to
the Taylor coefficients of mu(x) = x + c2 x^2 + c3 x^3 + c4 x^4 + c5 x^5 + ...

Question: does perturbing each c_k change the topology of the attractor,
or only its rate?  We compute beta(Phi) symbolically to O(Phi^2) with all
coefficients kept generic and then read off which coefficients enter which
terms.
"""

import sympy as sp
from sympy import (symbols, sqrt, series, simplify, expand, diff,
                   Rational, oo, Poly, S)

u, Phi = symbols('u Phi', positive=True, real=True)
c2, c3, c4, c5 = symbols('c2 c3 c4 c5', real=True)

# --------------------------------------------------------------------------
# Generic mu to order 5
# --------------------------------------------------------------------------
mu  = u + c2*u**2 + c3*u**3 + c4*u**4 + c5*u**5
mup = diff(mu, u)

beta_u = -u*(mu - u*mup) / (mu*(mu + u*mup))

# Change variable u = sqrt(Phi(2+Phi))
beta_Phi = beta_u.subs(u, sqrt(Phi*(2 + Phi)))

print("=" * 74)
print("BETA(PHI) EXPANSION KEEPING c2, c3, c4, c5 GENERIC")
print("=" * 74)

ser = series(beta_Phi, Phi, 0, 3).removeO()
ser = expand(ser)
print(f"\nbeta(Phi) = {ser}\n")

# Extract terms
print("Coefficients of beta(Phi) at each power of Phi:")
for pw_label, pw in [("sqrt(Phi)", sp.Rational(1,2)),
                     ("Phi",       1),
                     ("Phi^(3/2)", sp.Rational(3,2)),
                     ("Phi^2",     2)]:
    if pw == sp.Rational(1,2):
        coef = ser.coeff(sqrt(Phi))
    elif pw == sp.Rational(3,2):
        coef = ser.coeff(Phi**sp.Rational(3,2))
    else:
        coef = ser.coeff(Phi, pw)
    coef = simplify(coef)
    print(f"  {pw_label:<12s} = {coef}")

# --------------------------------------------------------------------------
# Structural analysis
# --------------------------------------------------------------------------
print("\n" + "=" * 74)
print("STRUCTURAL SENSITIVITY OF THE FIXED POINT PHI = 0")
print("=" * 74)

print("""
When c2 = 0 (the case relevant to surviving interpolants), the flow near
Phi = 0 becomes

    beta(Phi) = 2 c3 * Phi + a2(c3, c5) * Phi^2 + O(Phi^3),

so the *topology* of the fixed point depends on the sign of c3:

    c3 < 0  : linear attractor, exponential relaxation Phi ~ Phi_0 e^{2 c3 lambda}
    c3 = 0  : degenerate, marginal fixed point (quadratic approach at best)
    c3 > 0  : linear repulsor, exponential departure

Perturbations c2 -> c2 + delta with delta != 0 (breaking c2 = 0):
    Introduce a sqrt(Phi) term that DOMINATES the linear term near Phi = 0.
    Qualitative change: finite-time collapse (if delta < 0) or unbounded
    growth (if delta > 0) replaces the exponential relaxation.
    ==> c2 = 0 is a KNIFE-EDGE condition; any nonzero c2, no matter how
    small, changes the topology.

Perturbations c3 -> c3 + eps (with c2 = 0 preserved):
    Smooth change in beta'(0) = 2c3.
    Preserves attractor topology iff sign(c3 + eps) = sign(c3).
    For the standard c3 = -1/2, perturbations |eps| < 1/2 preserve
    the attractor; larger perturbations can flip it.

Perturbations c4 -> c4 + eps (with c2 = 0):
    From the expansion, c4 enters ONLY the Phi^(3/2) coefficient,
    hence subdominant near Phi = 0. Does not change fixed-point topology
    to leading order.

Perturbations c5 -> c5 + eps (with c2 = 0):
    c5 enters the Phi^2 coefficient (hence kappa). Changes curvature
    but not linear behaviour or topology.
""")

# --------------------------------------------------------------------------
# Numerical illustration for the surviving interpolants
# --------------------------------------------------------------------------
print("=" * 74)
print("CHECK: MARGINAL CASES mu_alpha = x/(1+x^alpha)^(1/alpha) FOR alpha=3,4")
print("=" * 74)

for alpha_val in [2, 3, 4, 5]:
    mu_a = u / (1 + u**alpha_val)**sp.Rational(1, alpha_val)
    ser_a = series(mu_a, u, 0, 6).removeO()
    c2_a = ser_a.coeff(u, 2)
    c3_a = ser_a.coeff(u, 3)
    c4_a = ser_a.coeff(u, 4)
    c5_a = ser_a.coeff(u, 5)
    print(f"  alpha = {alpha_val}:  c2 = {c2_a}, c3 = {c3_a}, c4 = {c4_a}, c5 = {c5_a}")
    if c3_a == 0:
        # degenerate case: c3 = 0
        print(f"           => beta'(0) = 2*c3 = 0  (marginal/degenerate)")
    else:
        print(f"           => beta'(0) = 2*c3 = {2*c3_a}")

# --------------------------------------------------------------------------
# Formal derivatives / linear sensitivities of observable quantities
# --------------------------------------------------------------------------
print("\n" + "=" * 74)
print("LINEAR SENSITIVITIES: HOW OBSERVABLES RESPOND TO PARAMETER CHANGES")
print("=" * 74)

# For sigma_Phi(z) = sigma_Phi(0) * (H/H0)^{beta'(0)}, take derivatives
# with respect to c3 (which controls beta'(0) = 2*c3 when c2 = 0):
print("""
For an interpolant with c2 = 0, the BTFR scatter evolves as
    sigma_Phi(z) = sigma_Phi(0) * (H(z)/H0)^{beta'(0)} = sigma_0 * xi^{2 c3},
with xi = H(z)/H0. Linear sensitivity to c3:

    d ln sigma / d c3 = 2 * ln xi.

At z = 2 (xi = 2.966 for Omega_m = 0.3), d ln sigma / d c3 = 2 * 1.087 = 2.17.
A 10 percent shift in c3 (i.e., delta c3 = 0.05 from c3 = -1/2 -> -0.55)
produces a fractional shift in sigma_Phi(z=2) of
    delta ln sigma = 2.17 * 0.05 = 0.109  (about 11 %).

Sensitivity to c2 (breaking the c2 = 0 condition) is NOT linear; it is
DISCONTINUOUS: any c2 != 0 introduces the sqrt-branch and produces
finite-time collapse.  This is precisely why the criterion is sharp.

Sensitivity to c4 and c5: does not enter beta'(0), hence does not affect
sigma_Phi to leading order.  Enters only via higher-Phi corrections
(intermediate-scale phenomenology, not the fixed-point structure).
""")
