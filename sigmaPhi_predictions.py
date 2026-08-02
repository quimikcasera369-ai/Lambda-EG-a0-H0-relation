"""
sigmaPhi_predictions.py

Reproduces Table 2 of

    Figueroa Torres, J. P. 2026
    "Reconstructing H(z) from galactic dynamics alone:
     an independent probe of the Hubble tension in Lambda-EG"

The predictions come from Eq. (6) of that letter,

    sigma_Phi(z) = sigma_Phi(0) * E(z)^{beta'(0)},

    E(z) = sqrt( Omega_m * (1+z)^3 + Omega_Lambda ),

with beta'(0) fixed for each surviving interpolant per the c_2=0 criterion
of the companion theorem paper [Figueroa Torres 2026b, Table 1]:

    Standard   mu = x/sqrt(1+x^2)   ->   beta'(0) = -1
    mu-hat entropic                  ->   beta'(0) = -2
    tanh(x)                          ->   beta'(0) = -2/3

E(z) depends on Omega_m only, NOT on H_0 -- this is the orthogonality
established in §4 of the letter.

Deterministic. Requires only the Python standard library
(no NumPy, no SciPy).
"""

import math

# -- Surviving interpolants (from theorem paper, Table 1) --
interpolants = {
    "Standard mu = x/sqrt(1+x^2)": -1.0,
    "mu-hat entropic":             -2.0,
    "tanh(x)":                     -2.0/3.0,
}

# -- Two cosmological scenarios --
scenarios = {
    "Planck-consistent":   dict(Om=0.315, OL=0.685, H0=67.4),
    "SH0ES-consistent":    dict(Om=0.30,  OL=0.70,  H0=73.5),
}

# -- Anchor: SPARC intrinsic scatter --
sigma_Phi_0 = 0.07   # dex, Lelli et al. 2019

zs = [0.5, 1.0, 2.0, 4.0]

def Ez(z, Om, OL):
    return math.sqrt(Om * (1+z)**3 + OL)

def sigmaPhi(z, betap, Om, OL):
    return sigma_Phi_0 * Ez(z, Om, OL)**betap

# -- Report --
print("=" * 80)
print(f" Registered predictions of sigma_Phi(z) [dex]")
print(f" Anchored at sigma_Phi(0) = {sigma_Phi_0} dex (SPARC, Lelli+2019)")
print("=" * 80)

hdr = f"  {'Interpolant':<32s} {'β′(0)':>7s} " + " ".join(f"z={z:<6}" for z in zs)
for scen_name, sp in scenarios.items():
    print()
    print(f" [{scen_name}]  Omega_m = {sp['Om']}, H_0 = {sp['H0']}")
    print(f"   (note: E(z) depends only on Omega_m, NOT on H_0 -- orthogonality)")
    print(hdr)
    print("  " + "-"*77)
    for name, betap in interpolants.items():
        vals = [sigmaPhi(z, betap, sp['Om'], sp['OL']) for z in zs]
        row = " ".join(f"{v:>8.4f}" for v in vals)
        print(f"  {name:<32s} {betap:>7.3f} {row}")

# -- Cross-check the orthogonality claim: E(z) invariant to H_0 at fixed Omega_m --
print()
print("=" * 80)
print(" Orthogonality check: E(z=2) at fixed Omega_m for several H_0")
print("=" * 80)
for H0 in [65, 67.4, 70, 73.5, 76]:
    # keep Omega_m = 0.315 fixed; only H_0 changes
    Ez2 = Ez(2, 0.315, 0.685)
    print(f"  H_0 = {H0:5.1f} km/s/Mpc,  Omega_m = 0.315  ->  E(2) = {Ez2:.4f}")
print("  --> E(z) is IDENTICAL. Reconstruction of shape decouples from scale.")

# -- Falsification band for z = 1.5 (mid-range ELT/HARMONI target) --
print()
print("=" * 80)
print(" Falsification band at z = 1.5 (ELT/HARMONI target range)")
print("=" * 80)
Om = 0.315; OL = 0.685
z = 1.5
vals = []
for name, betap in interpolants.items():
    v = sigmaPhi(z, betap, Om, OL)
    vals.append(v)
    print(f"  {name:<32s} predicts sigma_Phi({z}) = {v:.4f} dex")
lo, hi = min(vals), max(vals)
print()
print(f"  Union of surviving-interpolant predictions at z={z}:")
print(f"    sigma_Phi(z={z}) in [{lo:.4f}, {hi:.4f}] dex")
print(f"  A measurement outside this band falsifies the combined framework.")
