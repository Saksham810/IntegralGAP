# Rational table format

The single input is
`verification/integral/research/runs/eplus1_search_expanded5_cert.json.gz`.
It is a gzip-compressed UTF-8 JSON object. Python's `gzip` and `json` modules
read it directly. Its compressed-byte SHA-256 is:

```text
ed6fe9bcdcb3f24cacbc6b5a02cb823fb6f57bb124066c1e7d7e5f7d6ab4d264
```

The compressed file is 1,982,562 bytes; its uncompressed JSON is 10,823,421 bytes.
Preserve the compressed bytes: recompression can change the hash even if the
JSON has the same meaning.

## Parameters and indexing

The selected `summary` specifies:

```text
schema         = finecap-moment-cross2-strict-v1
K              = 6
h              = 2
factor         = 4
M              = 16
T              = 17/4
denominator Q  = 100000000000
target_degree  = 3
availability   = subdiv
variance_bound = exact
tight_B        = 1
```

All arrays listed below are lists of `[key, integer_numerator]` rows. Nested
list keys are converted to tuples by the programs. The value of a row is its
integer numerator divided by `Q`. Positive-mass and crossing arrays must be
nonnegative; the lower-function arrays permit signed integers.

All indices are zero-based except size class `j`, which runs from `1` to `6`.

| Symbol | Meaning |
| --- | --- |
| `b` | Time block, `0 <= b < M`; block length is `Delta = T/M = 17/64`. |
| `e` | Quadratic Bernstein control index, `0`, `1`, or `2`. |
| `i` | Global cubic control index for target tails, `0 <= i <= 3*M`. |
| `r` | Positive-load cell index, physical interval `(r/K, (r+1)/K)`. |
| `l` | Load-grid index `0 <= l <= K`. |
| `j` | Size class, physical sizes `((j-1)/K, j/K]`. |
| `ep` | Distinguished size endpoint: `0` for `(j-1)/K`, `1` for `j/K`. |
| `f` | Full large-item pattern: four counts for classes `3,4,5,6`. |
| `c` | Observed large-item pattern, in the same four-count order. |
| `label` | Nested key `[c, u, a]`, described next. |

Define `W(c) = sum((j-1)*c_j for j in 3,4,5,6)`. A full pattern satisfies
`W(f) <= K-1`. The seven possible patterns are the empty pattern, one item from
any large class, two class-3 items, and a class-3/class-4 pair.

For a full pattern `f` and observed pattern `c <= f`, the unseen-large reserve is
`u = W(f)-W(c)`. Write `R0 = K-W(c)-u`. The final label coordinate obeys
`0 <= a <= factor*R0`. For observed small total `S`, the comparison chooses
`a = ceil(K*factor*S)`, using `a=0` for `S=0`. A class-`j` large rule is eligible
when `c_j > 0`; a small rule is eligible when `factor*(j-1) < a`.

These labels are comparison indices in the analysis. They do not give the online
algorithm information about unseen items.

## Arrays

| Field | Key | Rows | Meaning after division by `Q` |
| --- | --- | ---: | --- |
| `P` | `[i, l]` | 343 | Cubic-in-time control values for target tails at load-grid points. Block `b` uses `P[3*b+q,l]`, `q=0,1,2,3`, so adjacent blocks share their endpoint controls. |
| `AL` | `[b, e, label, j]` | 13,920 | Accepted atom-at-zero mass at the lower endpoint of size class `j`. |
| `A` | `[b, e, label, j]` | 13,920 | Accepted atom-at-zero mass at the upper size endpoint. |
| `G` | `[b, e, label, j, r]` | 69,312 | Accepted mass of the full load cell `r`, before truncation at the true fitting boundary; `0 <= r <= K-j`. |
| `L` | `[b, e, label, j, l]` | 97,440 | Crossing-envelope value at the lower size endpoint and load-grid point `l`. |
| `U` | `[b, e, label, j, l]` | 97,440 | Crossing-envelope value at the upper size endpoint and load-grid point `l`. |
| `minor` | `[f, c, j, b, e, ep]` | 3,456 | Two size-endpoint coefficients of the affine term `beta(s)` in a lower function. |
| `slope` | `[f, c, j, b, e]` | 1,728 | Nonpositive linear small-mass coefficient `v`. |
| `curve` | `[f, c, j, b, e]` | 1,728 | Nonpositive quadratic small-mass coefficient `xi`. |

The full replay reconstructs each required key set from definitions and checks
for missing, extra, or repeated keys. It does not merely accept these row counts.

### Tail probabilities

At a given time, the target distribution is determined by:

```text
mass at zero             = 1 - P_0
mass in load cell r      = P_r - P_(r+1)
mass at one              = P_K
```

The density is constant within each open cell. The tail is affine between load
grid points. `P_K` represents the mass at one, equivalently the left limit of the
tail there; the strict tail beyond load one is zero.

### Accepted masses and crossing bounds

For true size `s` in class `j`, let `x = K*s-(j-1)`. Evaluating the stored
quadratic time polynomials first, the accepted atom at zero is
`(1-x)*AL + x*A`. Accepted density on full cell `r` is `K*G_r`.
It is restricted to positive starting loads below `1-s`. Thus cells
`0,...,K-j-1` contribute fully and the last cell contributes `(1-x)*G_(K-j)`.

The total acceptance mass is therefore:

```text
q(s) = (1-x)*AL + x*A
       + sum(G_r for r in 0,...,K-j-1) + (1-x)*G_(K-j)
```

The crossing envelope is affine in size between `L` and `U` and affine in load
between adjacent grid values. These are unconditional accepted masses and
crossing probabilities. Operational acceptance coins are obtained by dividing
accepted mass/density by the available target mass/density where positive.

### Lower functions

For each `(f,c,j,b,e)`, the two `minor` values define `beta(s)` by affine
interpolation in size. Together with `slope` and `curve`, they define:

```text
g(s,S) = beta(s) + v*S + xi*S*S,  with v <= 0 and xi <= 0
```

Here `S` is the total observed small-item size, including the distinguished
item if it is small. The lower function is shared across the compatible
small-mass labels. Its inequality is verified for all applicable labels and
continuous sizes using Bernstein controls.

## Reserve slopes and summary

`correction_slopes` is a separate list of **16 exact rational strings**. These
strings are parsed directly as fractions; **do not divide them by `Q` again**.
They specify a nonnegative reserve `R_b` on each time block. The checked
normalization is:

```text
H_star = (T/M)*sum(R_b) = 23407/5120000000000
alpha  = 1/(1+H_star)   = 5120000000000/5120000023407
```

The `summary` object records parameters, claimed exact constants, historical
counts, timings, and display decimals. Its stored `status: PASS` is not a
verification result for a new run. The checker reconstructs the numerical
inequalities and tests the claims against them. Display decimals and elapsed
times are never used to accept a mathematical bound.
