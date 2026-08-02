"""
H0_from_a0.py

Reproduces Eq. (3) and Table 1 of

    Figueroa Torres, J. P. 2026
    "Reconstructing H(z) from galactic dynamics alone:
     an independent probe of the Hubble tension in Lambda-EG"

The evaluation is a direct substitution of the SPARC-canonical a_0 into

    H_0 = pi * sqrt(3) * a_0 / c

which is the immediate consequence of the two identities
  (i)   g_crit = c * H_0 / (4 * pi^2 * sqrt(3))       [framework paper]
  (ii)  a_0    = 4 * pi * g_crit                       [companion paper, §4]

Deterministic. Requires only the Python standard library
(no NumPy, no SciPy).

NOTE: the value of H_0 returned here is NOT an H_0-independent measurement.
See a0_H0_selfconsistency.py and Section 4 of the letter: 55% of the SPARC
sample carries distances computed under an assumed H_0 = 73 km/s/Mpc.
"""

import math

# -- physical constants (SI) --
c = 2.99792458e8        # m/s
Mpc_m = 3.0856775814913673e22
kmMpc_to_s = 1.0e3 / Mpc_m  # (km/s/Mpc) -> 1/s

# -- SPARC-canonical a_0 measurement --
# Lelli, McGaugh, Schombert 2016; McGaugh, Lelli, Schombert 2016.
# Stat uncertainty ~2%, systematic (M/L, inclination) ~20%.
a0_central = 1.20e-10          # m/s^2
a0_stat    = 0.02e-10          # m/s^2, ~2% statistical
a0_sys     = 0.24e-10          # m/s^2, ~20% systematic (dominant)

# -- Eq. (3): H_0 = pi * sqrt(3) * a_0 / c --
def H0_from_a0(a0):
    """Return H_0 in km/s/Mpc given a_0 in m/s^2."""
    H0_si = math.pi * math.sqrt(3) * a0 / c     # 1/s
    return H0_si / kmMpc_to_s                    # km/s/Mpc

H0_central = H0_from_a0(a0_central)
H0_stat    = H0_from_a0(a0_stat)       # linear propagation
H0_sys     = H0_from_a0(a0_sys)

# -- Reference determinations --
H0_Planck = 67.4;  H0_Planck_err = 0.5   # Planck Collaboration 2020
H0_SH0ES  = 73.5;  H0_SH0ES_err  = 0.7   # Riess et al. 2024

# -- Report --
print("=" * 62)
print(" Eq. (3):  H_0  =  pi * sqrt(3) * a_0 / c")
print("=" * 62)
print(f"  pi * sqrt(3)              = {math.pi*math.sqrt(3):.6f}")
print(f"  a_0 (SPARC canonical)     = {a0_central:.3e} m/s^2")
print(f"      +/- {a0_stat:.2e} (stat)  +/- {a0_sys:.2e} (sys)")
print()
print(f"  H_0 (central)             = {H0_central:.2f} km/s/Mpc")
print(f"      +/- {H0_stat:.2f} (stat)  +/- {H0_sys:.2f} (sys)")
print()
print("=" * 62)
print(" Table 1 comparison")
print("=" * 62)
print(f"  {'Method':<26s} {'H_0 [km/s/Mpc]':<20s} Anchor")
print(f"  {'-'*26} {'-'*20} {'-'*22}")
print(f"  {'Planck 2018 CMB':<26s} {H0_Planck:5.1f} +/- {H0_Planck_err:.1f}   "
      "recombination")
print(f"  {'SH0ES 2024':<26s} {H0_SH0ES:5.1f} +/- {H0_SH0ES_err:.1f}   "
      "Cepheid ladder")
print(f"  {'This letter (Eq. 3)':<26s} {H0_central:5.1f} (central)  "
      "SPARC a_0")
print()
print(f"  Delta vs. Planck   = {(H0_central - H0_Planck)/H0_Planck*100:+.2f}%")
print(f"  Delta vs. SH0ES    = {(H0_central - H0_SH0ES )/H0_SH0ES *100:+.2f}%")
print()
print("=" * 62)
print(" Sensitivity: what a_0 would each reference H_0 imply?")
print("=" * 62)
# invert Eq. (3):  a_0 = c * H_0 / (pi*sqrt(3))
def a0_from_H0(H0_kmMpc):
    return c * H0_kmMpc * kmMpc_to_s / (math.pi * math.sqrt(3))
for label, H0 in [("Planck", H0_Planck), ("SH0ES", H0_SH0ES)]:
    a0_req = a0_from_H0(H0)
    delta = (a0_req - a0_central) / a0_central * 100
    print(f"  {label} H_0 = {H0}    -->  a_0 = {a0_req:.3e} m/s^2  "
          f"(Delta = {delta:+.1f}% vs SPARC)")
