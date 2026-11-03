# How the finite checks enter the integral proof

The reference is the supplied *Improved Algorithms in the Random-Order Model:
Knapsack, GAP, and Beyond*. Appendix B constructs the one-resource rules and
proves why finitely many checks suffice. Appendix C transfers the resulting
coefficient to integral GPAP. The repository reproduces the integral numerical
work described in Appendix E.1 and the relevant scope statement in E.3.

## Five proof obligations

| Obligation | Requirement | Implementation | Relevant manuscript material |
| --- | --- | --- | --- |
| V1: distribution | Target tails define probability distributions and are monotone in load and time. | Tail range and monotonicity checks in `push027_cross_strict_replay.py`; block endpoints share the same stored controls. | B.5.2 |
| V2: availability | The accepted measure never exceeds the available target measure and is supported where the item fits. | Full replay and `exact_geometry.availability_controls`; truncation in the rule definition. | B.5.2-B.5.3 |
| V3: single-rule crossing | Every stored envelope bounds the true crossing probability for all real sizes and thresholds. | `push027_cross_geometry.spatial_checks` plus `eplus1_audit_true_spatial.py`. | B.6.1, Lemma B.4 |
| V4: packing | Every observed part of a feasible configuration has an eligible rule assignment within the threshold budget. | `push027_strict_geometry.packing_values`, compared with tail-derivative controls plus reserves in the replay. | B.6.2-B.6.3, Lemmas B.5-B.6 |
| V5: value | Every distinguished item has the claimed comparison coefficient. | Lower-function tests and signed integral sums in the replay, checked separately by `verify_service_independent.py`. | B.7-B.8, especially equations B.38-B.48 |

The scalar constants for passing from the limiting calculation to finite inputs
are also checked. The inequalities relating them to online performance are
proved in B.9-B.11.

## Availability over continuous time

Target masses are cubic Bernstein polynomials and accepted masses are
quadratic. The quadratic is elevated to degree three. For cubic available
controls `m0,...,m3` and quadratic accepted controls `a0,a1,a2`, the difference
controls are:

```text
d0 = m0-a0
d1 = m1-(a0+2*a1)/3
d2 = m2-(2*a1+a2)/3
d3 = m3-a2
```

The verifier tests seven weighted sums corresponding to subdivision of the
cubic onto its two half-blocks. Their nonnegativity proves availability at
every time, because the Bernstein basis is nonnegative. Testing the two
size-endpoint atom masses covers intermediate sizes by interpolation.
Truncation at `1-s` only removes accepted density.

## Crossings over continuous size and load

In scaled coordinates, put `sigma=K*s`, `z=K*physical_threshold`,
`x=sigma-(j-1)`, and `y=z-l` in a size/load cell. The conservative crossing
expression is:

```text
B(z,sigma) = 1_{sigma>0 and z<=sigma}
             * ((j-sigma)*AL + (sigma-j+1)*A)
             + sum_r G_r * max(0,
                 min(r+1,z,K-sigma) - max(r,z-sigma,0))
```

At an atom boundary `z=sigma`, including the atom gives a conservative bound
for the actual half-open crossing event. The density expression changes across
the diagonal `x=y` and, in the final fitting cell, `x+y=1`.
Inside the resulting triangles, the crossing expression is affine; the stored
envelope is bilinear. The manuscript shows that their difference is controlled
by boundary pieces. Axis-parallel edges are affine; diagonal pieces are
quadratic. The main replay uses exact Bernstein inequalities on subdivided
diagonals, applied coefficientwise in time.

The separate spatial audit uses a different calculation on those same proven
polynomial pieces. It reconstructs each open quadratic from the interior
parameters `1/4`, `1/2`, `3/4`, checks its identity again at `1/3`, and evaluates
its exact minimum, including an interior stationary point when relevant. It
checks actual closed-boundary values separately from one-sided polynomial
limits so an atom jump is not interpolated through an open segment.

This is not an argument from a few sampled points. The proof of the piecewise
polynomial structure is what makes finite polynomial checks sufficient.

## Packing every configuration

Large items have size greater than `1/3`; there are at most two. The seven full
large-class patterns and all compatible observed patterns can therefore be
enumerated. A label records the observed pattern, a lower-endpoint reserve for
unseen large items, and a cap for observed small mass.

For a fixed label, threshold, and time control, crossing cost is affine in item
sizes. The endpoint argument in B.6 reduces its maximum to a grid of spacing
`1/(K*factor)`. A dynamic program computes the maximum cost of any multiset of
positive-weight small types under the capacity cap. Large-item offsets are
enumerated explicitly. Zero-weight class-1 endpoint types have zero cost and
can be omitted, so no bound on the number of small items is needed.

The code handles strict lower endpoints explicitly. For example, two actual
class-2 items both exceed `1/6`, so they cannot fill a small-mass cap of exactly
`1/3`. The sharp boundary cases exclude only endpoint patterns that the
manuscript proves unrealizable. Full-cap labels use the coarser endpoint grid
and its corresponding strict-boundary cases.

Every resulting maximum is compared with the appropriate quadratic coefficient
of the cubic target-tail derivative plus the nonnegative reserve slope.
Affine interpolation in load and the nonnegative Bernstein basis extend these
inequalities over all thresholds and times. Thus the finite check covers all
feasible configurations, rather than selected test packings.

## Value coefficients and signed integrals

The acceptance probability is lower-bounded by
`beta(s)+v*S+xi*S^2`, with `v<=0` and `xi<=0`. The full replay checks this
inequality for every compatible small-mass label. It tests endpoint and middle
Bernstein controls after splitting at the relevant size kinks.

For small companion sizes bounded by `delta=1/3` and total `C`, the largest
possible sum of squares is:

```text
J(C) = floor(C/delta)*delta^2
       + (C-floor(C/delta)*delta)^2
```

This bounds the first two moments of observed small mass. Because `v` and `xi`
are nonpositive, upper moment bounds give a valid lower value bound. Large
companions are enumerated, and the resulting expression is integrated over
the prefix fraction `p`, with logarithmic construction time `t=log(p)+T`.
Concavity in true size reduces each eligible pattern/class pair to its two
endpoints: 19 pairs, hence 38 coefficients. Actual zero-size items use a
separate free rule; the additional check is `1-exp(-T) >= gamma_star`.

The first integration method expands the powers of `1-p` and uses exact
endpoint antiderivatives, with rational enclosures of exponentials. The
independent method enumerates large companions by individual identity and
uses the positive Taylor/Beta expansion from B.8.2. It truncates the integral
series after degree 80, encloses the exponential prefactor through degree 220,
and rounds rational intervals outward at scale `10^70`.

Both methods preserve coefficient signs in interval sums. For a negative
coefficient, a lower product bound uses the **upper** integral endpoint.
Using all lower endpoints would be unsound. Decimal values in the output are
display aids; mathematical acceptance uses integers and `Fraction` objects.

The independent program reuses the selected lower-function data, but imports
none of the replay's code. It provides a second calculation of their integrals
and finite constants, not a substitute for V1-V4 or the lower-function checks.

## Reserves and finite-input bounds

The first reserve moves extra mass from zero to one throughout the construction.
The exact normalization factor is `alpha=1/(1+H_star)`, where
`H_star=(T/M)*sum(R_b)`. Both programs recompute it from the rational reserve
slopes. This normalization is already included in `gamma_star`.

For finite inputs, B.9 additionally accounts for a rational harmonic clock,
sampling without replacement, rounded block boundaries, discarded boundary
rounds, and a finite-drift reserve. The verifier checks bounds on:

```text
d_max  = 2
H_coef = 1200000000001/100000000000
C_tail = 8549714592/14111328125
E      = 71 > exp(T)
```

Here `H_coef` bounds the sum of absolute lower-function controls, and `C_tail`
bounds second time derivatives of the original target tails. The finite-loss
expression is:

```text
A_required = H_coef*(4*(d_max+2)^2 + (20*E+20)*M + 20)
             + (E/2)*T*C_tail
```

For this table the reported constant
`A_fin=4718849996128173/17000000000` equals that expression. The replay and
independent audit check that it is sufficient. The new root runner substitutes
both horizon/failure pairs and verifies strict rational margins above the
displayed targets.

The code checks these constants and substitutions. The manuscript proves the
uniform loss `A_fin/N`, the dummy-insertion law, the implementation of padding
and capped draws, and the global sampling-failure loss. None of those arguments
is established merely by finding a `PASS` field in a stored JSON file.

## Transfer to GPAP

Appendix C proposes one action per item and applies one-resource rules to its
components. The resource-deletion property turns the successful component tests
into an available action; linearity of expectation combines their value
contributions. No independence of resource outcomes is required for that value
sum. This analytical step extends the checked one-resource coefficient to
integral GPAP without an extra factor depending on the number of resources.

The benchmark remains the offline **integral** optimum. No fractional-GPAP or
hardness claim is verified by this repository.
