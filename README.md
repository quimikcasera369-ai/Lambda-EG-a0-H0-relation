# The a₀–H₀ relation in Λ-EG

**A consistency test and its distance-ladder circularity**

J. P. Figueroa Torres — Independent Researcher
ORCID: [0009-0005-5297-8777](https://orcid.org/0009-0005-5297-8777)

---

## What this is

The Λ-EG framework fixes two acceleration scales in a rigid geometric relation:
the screening scale `g_crit = cH₀/(4π²√3)` and the MOND acceleration
`a₀ = 4π·g_crit`. Eliminating `g_crit` between them gives, in one line of algebra,

```
H₀ = π√3 · a₀ / c
```

Substituting the SPARC-canonical `a₀ = 1.20 × 10⁻¹⁰ m/s²` returns
**H₀ = 67.21 km/s/Mpc**, within 0.3 % of the Planck 2018 central value.

**That agreement cannot currently be read as an independent determination of H₀,
and the main result of this letter is the quantification of why.**

## The cautionary result

SPARC Table 1 carries a distance-method flag whose first value is defined as
*"Hubble-Flow assuming H0 = 73 km/s/Mpc"*. Counting the standard BTFR sample
(Q ≤ 2 with measured V_flat, 129 galaxies):

| Distance method | N | fraction |
|---|---|---|
| Hubble flow (assumes H₀ = 73) | 71 | 55.0 % |
| Ursa Major cluster (all assigned 18.0 Mpc) | 25 | 19.4 % |
| TRGB | 28 | 21.7 % |
| Cepheids | 3 | 2.3 % |
| Supernovae Ia | 2 | 1.6 % |

Between **55 % and 74 %** of the sample (depending on whether the cluster tier is
counted) carries distances computed under an assumed H₀ — and that assumed value,
73, sits at the SH0ES end of the tension. Only **33 galaxies** have individually
measured redshift-independent distances.

Because the BTFR contains no length scale (the radius cancels identically in
`V⁴/r² = a₀GM/r²`), the distance enters through the baryonic mass alone, with the
single exact power `M ∝ D²`. Therefore

```
a₀ ∝ (H₀ assumed)²        d ln a₀ / d ln H₀ = 2f = 1.10 (loose) to 1.49 (strict)
```

Imposing the relation and the calibration simultaneously gives a loop that is
**singular at H₀ ≈ 72 km/s/Mpc** — inside the range under dispute. It determines
no value of H₀ at all.

The companion observable `σ_Φ(z)` is immune: being a ratio, H₀ cancels identically.

## Reproducing the results

All scripts are deterministic and require **only the Python standard library**.

```bash
python3 H0_from_a0.py              # Eq. (3) and Table 1
python3 sigmaPhi_predictions.py    # Table 3, registered σ_Φ(z) predictions
python3 a0_H0_selfconsistency.py   # Section 4: Table 2, the coupling, the loop
```

`sensitivity.py` (Appendix A.11, symbolic fixed-point analysis) additionally
requires SymPy.

`a0_H0_selfconsistency.py` reads `SPARC_Lelli2016c.mrt`, included here for
reproducibility. That table is **not ours**: it is machine-readable Table 1 of
Lelli, McGaugh & Schombert (2016), AJ 152, 157, distributed by the SPARC
project. Please cite them for any use of the data.

## Contents

| File | What it is |
|---|---|
| `Lambda-EG_H0_reconstruction_letter.tex` / `.pdf` | The letter |
| `H0_from_a0.py` | Evaluation of `H₀ = π√3 a₀/c` and comparison table |
| `sigmaPhi_predictions.py` | Registered `σ_Φ(z)` predictions for the surviving interpolants |
| `a0_H0_selfconsistency.py` | The distance-provenance audit and the self-consistency loop |
| `AppendixA11_sensitivity.tex`, `sensitivity.py` | Appendix A.11 — fixed-point sensitivity to the Taylor coefficients |
| `SPARC_Lelli2016c.mrt` | SPARC Table 1 (Lelli et al. 2016), for reproducibility |

## What is claimed, and what is not

**Claimed.** That Λ-EG fixes a rigid algebraic relation between a₀ and H₀ with no
free parameters; that evaluating it on current data returns a value 0.3 % from
Planck; that this is suggestive but cannot count as evidence until the
circularity is removed; that the route to removing it is concrete (the 33
redshift-independent galaxies); and that the `σ_Φ(z)` channel is genuinely
H₀-independent.

**Not claimed.** A resolution of the Hubble tension. A third independent
determination of H₀. An earlier version of this letter asserted that the SPARC a₀
is independent of any assumption on H₀; that assertion was wrong and is retracted
in §6.

Eq. (3) also presupposes the branch `a₀ ∝ H(z)`. The branch `a₀ ∝ √Λ` returns
H₀ = 81.3 for the same a₀ and is excluded only at 0.85σ of the systematic. The
value 67.2 is conditional on the branch.

## Related work

Milgrom (2009, ApJ 698, 1630) derives the deep-MOND limit from space-time scale
invariance and conjectures that local gravitational physics would take exactly
the deep-MOND form in an exact de Sitter universe. This letter is agnostic on
that mechanism: Eq. (3) is used only as an algebraic relation between two
measured scales.

## Companion papers

- Figueroa Torres, J. P. 2026, *Lambda EG: Emergent Gravitational Coherence and
  Galactic Dynamics Without Particle Dark Matter*, Zenodo preprint,
  [doi:10.5281/zenodo.19784004](https://doi.org/10.5281/zenodo.19784004)
- Figueroa Torres, J. P. 2026b, *A geometric characterization of the standard
  MOND interpolating function*

## License

Text and figures: CC BY 4.0. Code: MIT. See `LICENSE`.

The included SPARC table is the property of its authors and is redistributed here
under the terms of its original publication; it is not covered by the above.
