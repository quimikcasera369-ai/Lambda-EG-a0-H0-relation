"""
a0_H0_selfconsistency.py

Quantifies the H_0 content of the SPARC a_0 measurement, for Section 4 and
Table 2 of

    Figueroa Torres, J. P. 2026
    "Reconstructing H(z) from galactic dynamics alone:
     an independent probe of the Hubble tension in Lambda-EG"

The letter uses Eq. (3), H_0 = pi*sqrt(3)*a_0/c, to infer H_0 from a galactic
measurement of a_0.  That inference is only independent of the distance ladder
to the extent that the *distances* entering a_0 are themselves independent of
H_0.  They are not, for a majority of the sample.

SPARC Table 1 [Lelli, McGaugh & Schombert 2016] carries a distance-method
flag f_D whose value 1 is defined in the header as

    "1 = Hubble-Flow assuming H0=73 km/s/Mpc and correcting for
         Virgo-centric infall"

This script counts how much of the standard BTFR sample carries that flag and
propagates the resulting sensitivity of a_0 to the assumed H_0.

Scaling argument.  For a galaxy whose distance comes from the Hubble flow,
D ∝ 1/H_0^assumed.  Baryonic mass is inferred photometrically, so M ∝ D^2,
while V_flat is a velocity and carries no distance dependence.  Since
a_0 = V_flat^4/(G M) in the deep-MOND/BTFR normalisation,

    a_0 ∝ 1/M ∝ 1/D^2 ∝ (H_0^assumed)^2      [f_D = 1 galaxies only]

so d ln a_0 / d ln H_0 = 2 on that sub-sample, and equals the sample fraction
times 2 overall.

Requires only the Python standard library and the file SPARC_Lelli2016c.mrt
(machine-readable Table 1, distributed with this repository).
"""

import math
import statistics as st
from collections import Counter

MRT = "SPARC_Lelli2016c.mrt"

# -- physical constants (SI) --
c = 2.99792458e8
G = 6.674e-11
Msun = 1.98892e30
Mpc_m = 3.0856775814913673e22
kmMpc_to_s = 1.0e3 / Mpc_m

# -- SPARC conventions --
UPSILON_STAR = 0.5      # M/L at [3.6] um, canonical SPARC disk value
H0_ASSUMED = 73.0       # km/s/Mpc, hard-wired into the f_D = 1 distances

DISTANCE_METHOD = {
    1: "Hubble flow (assumes H0=73)",
    2: "TRGB",
    3: "Cepheids",
    4: "Ursa Major cluster",
    5: "Supernovae Ia",
}

# Three tiers of independence from H_0, in decreasing strictness.
#
#   PRIMARY  (f_D in {2,3,5}) : the distance is measured to the individual
#            galaxy by a stellar/standard-candle indicator. These are
#            "redshift-independent" in the NED-D sense. They are not free of
#            the distance-ladder zero point (TRGB and Cepheids share anchors
#            with SH0ES: NGC 4258, LMC detached eclipsing binaries, Gaia
#            parallaxes), but no value of H_0 enters the inference.
#
#   CLUSTER  (f_D = 4) : Ursa Major members. SPARC assigns ALL of them the
#            single cluster distance 18.0 Mpc; it is not an individual
#            measurement, and the cluster distance itself rests on a
#            Tully-Fisher calibration, which carries a distance-scale
#            dependence. We therefore treat this tier as ambiguous and
#            report results both ways.
#
#   FLOW     (f_D = 1) : D = v/H_0 with H_0 = 73 hard-wired. Unambiguously
#            H_0-dependent.
PRIMARY = {2, 3, 5}
CLUSTER = {4}
FLOW = {1}


def load_table(path=MRT):
    """Parse the machine-readable SPARC Table 1 into dicts."""
    lines = [l.rstrip("\n") for l in open(path)]
    # data begins after the last rule of the byte-by-byte header
    start = max(i for i, l in enumerate(lines) if l.startswith("-----"))
    out = []
    for line in lines[start + 1:]:
        tok = line.split()
        if len(tok) < 18:
            continue
        try:
            out.append(dict(
                name=tok[0], T=int(tok[1]), D=float(tok[2]), e_D=float(tok[3]),
                f_D=int(tok[4]), L36=float(tok[7]), MHI=float(tok[13]),
                Vflat=float(tok[15]), e_Vflat=float(tok[16]), Q=int(tok[17]),
            ))
        except ValueError:
            continue
    return out


def a0_galaxy(g):
    """a_0 = V_flat^4 / (G M_bar), with M_bar = Upsilon_* L_[3.6] + 1.33 M_HI."""
    M_bar = (UPSILON_STAR * g["L36"] + 1.33 * g["MHI"]) * 1e9 * Msun
    return (g["Vflat"] * 1e3) ** 4 / (G * M_bar)


def H0_from_a0(a0):
    """Eq. (3) of the letter, in km/s/Mpc."""
    return math.pi * math.sqrt(3) * a0 / c / kmMpc_to_s


gals = load_table()
# standard BTFR sample: quality flag 1-2, measured flat velocity
sample = [g for g in gals if g["Q"] <= 2 and g["Vflat"] > 0]
dep = [g for g in sample if g["f_D"] in FLOW]
ind = [g for g in sample if g["f_D"] not in FLOW]
primary = [g for g in sample if g["f_D"] in PRIMARY]
cluster = [g for g in sample if g["f_D"] in CLUSTER]

# f = fraction of the sample whose distance carries an assumed H_0.
# LOOSE  counts only the Hubble-flow galaxies.
# STRICT also counts the Ursa Major cluster (see the tier comment above).
f_loose = len(dep) / len(sample)
f_strict = (len(dep) + len(cluster)) / len(sample)
f_dep = f_loose      # the value quoted in the letter's headline number

print("=" * 74)
print(" 1. How much of the SPARC BTFR sample assumes a value of H_0?")
print("=" * 74)
print(f"  full SPARC Table 1                        : {len(gals)} galaxies")
print(f"  standard BTFR sample (Q <= 2, V_flat > 0) : {len(sample)} galaxies")
print()
counts = Counter(g["f_D"] for g in sample)
print(f"  {'distance method':<34s} {'N':>4s}  {'fraction':>9s}")
print(f"  {'-'*34} {'-'*4}  {'-'*9}")
for k in sorted(counts):
    print(f"  {DISTANCE_METHOD[k]:<34s} {counts[k]:4d}  {counts[k]/len(sample):8.1%}")
print()
print()
print("  Tiers of independence (see header comment):")
print(f"    PRIMARY  individual standard-candle distance : {len(primary):3d}")
print(f"    CLUSTER  Ursa Major, one shared distance     : {len(cluster):3d}")
print(f"    FLOW     D = v/H_0 at H_0 = 73               : {len(dep):3d}")
print()
print(f"  f (loose:  FLOW only)          = {len(dep)}/{len(sample)} = {f_loose:.1%}")
print(f"  f (strict: FLOW + CLUSTER)     = {len(dep)+len(cluster)}/{len(sample)} = {f_strict:.1%}")
print(f"  genuinely redshift-independent = {len(primary)} galaxies")
print()
print("  CHECK on the Ursa Major tier: distinct distance values assigned")
for k in sorted(set(g["f_D"] for g in sample)):
    D = sorted(set(g["D"] for g in sample if g["f_D"] == k))
    n = len([g for g in sample if g["f_D"] == k])
    tag = f"{len(D)} distinct value(s)" + (f"  = {D[0]} Mpc for all" if len(D) == 1 else "")
    print(f"    {DISTANCE_METHOD[k]:<30s} N={n:3d}  {tag}")
print("  -> the 25 Ursa Major galaxies share a single assigned distance; it is")
print("     a cluster membership assignment, not an individual measurement.")

print()
print("=" * 74)
print(" 2. Propagated sensitivity of a_0 to the assumed H_0")
print("=" * 74)
print("  On f_D = 1 galaxies:  D ∝ 1/H_0,  M ∝ D^2,  a_0 = V^4/(GM) ∝ H_0^2")
print()
print("  Why the exponent is exactly 2 and not fitted: in the deep-MOND regime")
print("  a^2 = a_0 a_N, so V^4/r^2 = a_0 G M/r^2 and the RADIUS CANCELS. The")
print("  BTFR has no length scale, so distance enters through M alone (V_flat")
print("  is a Doppler velocity, distance-free). Rescaling every distance by a")
print("  common lambda must therefore give a_0 ∝ lambda^-2 identically:")
print()
print(f"  {'lambda':>8s} {'median a_0':>13s} {'ratio':>10s} {'lambda^-2':>10s} {'exponent':>10s}")
base = st.median([a0_galaxy(g) for g in sample])
for lam in (0.90, 0.95, 1.05, 1.10, 1.20):
    def a0_scaled(g, L=lam):
        M = (UPSILON_STAR * g["L36"] + 1.33 * g["MHI"]) * L**2 * 1e9 * Msun
        return (g["Vflat"] * 1e3) ** 4 / (G * M)
    v = st.median([a0_scaled(g) for g in sample])
    print(f"  {lam:8.2f} {v:13.4e} {v/base:10.5f} {lam**-2:10.5f} "
          f"{math.log(v/base)/math.log(lam):10.5f}")
print()
print("  Dimensional note: in M = V^4/(G a_0) the powers of time cancel exactly")
print("  ([V^4] = m^4 s^-4, [G a_0] = m^4 kg^-1 s^-4), so the BTFR carries no")
print("  net time dimension. H_0 is an inverse time and therefore CANNOT enter")
print("  the relation structurally -- only through the flux-to-mass calibration,")
print("  which needs a distance. The circularity is a property of the data, not")
print("  of the relation.")
print()
print(f"  => d ln a_0 / d ln H_0 = 2f (sample-weighted):")
print(f"       loose  (f = {f_loose:.3f}) : {2*f_loose:.2f}")
print(f"       strict (f = {f_strict:.3f}) : {2*f_strict:.2f}")
print("     Either way the coefficient exceeds unity: a 1% error in the")
print("     assumed H_0 moves a_0 by more than 1%.")
print()
for H0_true in (67.4, 70.0, 73.0, 76.0):
    # a_0 rescales only on the f_D = 1 sub-sample
    factor = f_dep * (H0_true / H0_ASSUMED) ** 2 + (1 - f_dep)
    a0_corr = 1.20e-10 * factor
    print(f"  if the true H_0 were {H0_true:5.1f}:  a_0 -> {a0_corr:.3e} m/s^2 "
          f"({factor-1:+6.1%}),  Eq.(3) returns H_0 = {H0_from_a0(a0_corr):5.1f}")
print()
print("  The shift induced across the Planck-SH0ES range is comparable in size")
print("  to the tension itself. a_0 is therefore NOT an H_0-free observable at")
print("  present data quality, and Eq. (3) must be read as a consistency")
print("  relation rather than as an independent measurement.")

print()
print("-" * 74)
print("  Self-consistent solution")
print("-" * 74)
print("  Eq. (3) and the data-calibration coupling must hold simultaneously:")
print()
print("     H_0 = K * [ f*(H_0/73)^2 + (1-f) ],   K = pi*sqrt(3)*a_0^SPARC/c")
print()
K = H0_from_a0(1.20e-10)


def roots(f):
    """Roots of H = K[f(H/73)^2 + (1-f)]."""
    qa = K * f / H0_ASSUMED ** 2
    disc = 1.0 - 4 * qa * K * (1 - f)
    if disc < 0:
        return None
    return ((1 - math.sqrt(disc)) / (2 * qa), (1 + math.sqrt(disc)) / (2 * qa))


print(f"  K = {K:.3f} km/s/Mpc")
for lab, f in (("loose ", f_loose), ("strict", f_strict)):
    r = roots(f)
    print(f"  f = {f:.3f} ({lab}): roots H_0 = {r[0]:.1f} and {r[1]:.1f} km/s/Mpc")
print()
print("  Neither root is 67.2. The value quoted in the letter is the solution")
print("  of Eq. (3) at FIXED a_0, i.e. it holds a_0 constant while varying the")
print("  very quantity that a_0 was calibrated against. Reading 67.2 as an")
print("  H_0-independent determination is therefore not supported by the data")
print("  as currently distributed.")
print()
print()
print("-" * 74)
print("  Conditioning of the loop: how hard does H_0 respond to a_0?")
print("-" * 74)
print("  Writing H = K g(H) with g(H) = f(H/73)^2 + (1-f), the response of the")
print("  inferred H_0 to the published a_0 (which sets K) is")
print()
print("     dlnH_0/dln a_0 = (K/H) * g(H) / (1 - K g'(H)),   g' = 2fH/73^2")
print()
print(f"  {'f':>7s} {'H_0':>7s} {'1 - K g_prime':>15s} {'dlnH0/dln a0':>14s}")
print(f"  {'-'*7} {'-'*7} {'-'*15} {'-'*14}")
for f in (f_loose, f_strict):
    for H in (67.2, 70.0, 72.0, 73.0):
        g = f * (H / H0_ASSUMED) ** 2 + (1 - f)
        gp = 2 * f * H / H0_ASSUMED ** 2
        den = 1 - K * gp
        resp = (K / H) * g / den
        print(f"  {f:7.3f} {H:7.1f} {den:+15.4f} {resp:+14.2f}")
    H_sing = H0_ASSUMED ** 2 / (2 * f * K)
    print(f"  {'':7s} singular (denominator = 0) at H_0 = {H_sing:.2f}")
    print()
print("  With the loose f the loop is SINGULAR at H_0 = 72.1, i.e. inside the")
print("  range under dispute: the response changes sign and diverges there, so")
print("  no stable value of H_0 is determined at all. With the strict f the")
print("  singularity moves out to 53.3 and the response is a finite -3 to -2,")
print("  still large and NEGATIVE (a higher published a_0 implies a LOWER")
print("  self-consistent H_0, the opposite of the naive reading of Eq. 3).")
print("  Under either accounting the loop fails to determine H_0. That failure,")
print("  not any particular root, is the result of this section.")
print()
print("  CAVEAT: this fixed point uses a linear sample-weighted propagation.")
print("  The published a_0 comes from a global fit to the radial acceleration")
print("  relation, not from a median over Table 1, and the rescaling of the")
print("  f_D = 1 distances is not exactly a rigid multiplication (the")
print("  Virgo-centric infall correction is nonlinear). The roots above are")
print("  an order-of-magnitude statement about the size of the coupling,")
print("  NOT a competing determination of H_0.")

print()
print("=" * 74)
print(" 3. Direct test: a_0 on the two sub-samples")
print("=" * 74)
print(f"  estimator: a_0 = V_flat^4/(G M_bar), Upsilon_* = {UPSILON_STAR}")
print(f"  (this is the BTFR normalisation, numerically distinct from the RAR")
print(f"   scale g_dagger = 1.20e-10 quoted in the letter; only the RATIO")
print(f"   between sub-samples is used here)")
print()
print(f"  {'sub-sample':<30s} {'N':>4s} {'median a_0':>12s} {'scatter':>9s}")
print(f"  {'-'*30} {'-'*4} {'-'*12} {'-'*9}")
for label, S in (("f_D = 1 (assumes H0=73)", dep),
                 ("H_0-independent distance", ind),
                 ("all", sample)):
    v = [a0_galaxy(g) for g in S]
    lg = [math.log10(x) for x in v]
    print(f"  {label:<30s} {len(S):4d} {st.median(v):12.3e} "
          f"{st.pstdev(lg):8.3f} dex")
md = st.median([a0_galaxy(g) for g in dep])
mi = st.median([a0_galaxy(g) for g in ind])
print()
print(f"  offset (independent vs Hubble-flow) = {(mi-md)/md:+.1%}")

print()
print("  CONTROL: the same offset, in bins of V_flat (removes mass selection)")
print(f"  {'V_flat [km/s]':<16s} {'N dep':>6s} {'N ind':>6s} {'offset':>9s}")
print(f"  {'-'*16} {'-'*6} {'-'*6} {'-'*9}")
for lo, hi in ((0, 80), (80, 130), (130, 200), (200, 400)):
    d = [a0_galaxy(g) for g in dep if lo <= g["Vflat"] < hi]
    i = [a0_galaxy(g) for g in ind if lo <= g["Vflat"] < hi]
    if len(d) >= 3 and len(i) >= 3:
        off = f"{(st.median(i)-st.median(d))/st.median(d):+.1%}"
    else:
        off = "N insuf."
    print(f"  {f'{lo}-{hi}':<16s} {len(d):6d} {len(i):6d} {off:>9s}")
print()
print("  The offset does not hold a consistent sign across mass bins, so the")
print("  global +9.5% is NOT a clean measurement of the H_0 bias: it is")
print("  degenerate with mass-dependent selection between the sub-samples.")
print("  We report it as an upper bound on what this test can currently")
print("  resolve, not as a detection. The clean route is to re-derive a_0 on")
print(f"  the {len(ind)} H_0-independent galaxies with full mass modelling,")
print("  which requires the rotation curves rather than Table 1 alone.")
