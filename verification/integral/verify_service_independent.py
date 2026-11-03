#!/usr/bin/env python3
"""Independent exact audit of the selected service and finite-size constants.

No modules from research/ are imported. Large companions are enumerated as
individual identities, not by the source verifier's binomial multiplicities.
Exponential integrals use a positive Taylor/Beta series, not the source
verifier's endpoint antiderivatives. All acceptance tests use Fraction.
Decimal is used only to format output. This is not an online simulation.
"""
from __future__ import annotations

import argparse
from decimal import Decimal, localcontext
from fractions import Fraction as F
from functools import lru_cache
import gzip
import hashlib
from itertools import product
import json
from math import comb
from pathlib import Path
import sys
import time

Interval = tuple[F, F]
GRID = 10**70
ROOT = Path(__file__).resolve().parent
DEFAULT = ROOT / 'research/runs/eplus1_search_expanded5_cert.json.gz'


def require(test: bool, message: str) -> None:
    if not test:
        raise ValueError(message)


def freeze(value):
    return tuple(freeze(x) for x in value) if isinstance(value, list) else value


def rounded(interval: Interval) -> Interval:
    lo, hi = interval
    require(lo <= hi, 'Reversed interval')
    return (F((lo.numerator * GRID) // lo.denominator, GRID),
            F(-((-hi.numerator * GRID) // hi.denominator), GRID))


@lru_cache(maxsize=None)
def exp_interval(x: F) -> Interval:
    """Rational Taylor enclosure, with a geometric bound on the omitted tail."""
    if x < 0:
        lo, hi = exp_interval(-x)
        return rounded((1 / hi, 1 / lo))
    terms = 220
    require(x < terms + 2, 'Taylor tail denominator would not be positive')
    term = total = F(1)
    for k in range(1, terms + 1):
        term *= x / k
        total += term
    remainder = term * x / (terms + 1) / (1 - x / (terms + 2))
    return rounded((total, total + remainder))


def signed_sum(terms: list[tuple[F, Interval]]) -> Interval:
    lo = hi = F(0)
    for coefficient, (a, b) in terms:
        if coefficient >= 0:
            lo += coefficient * a
            hi += coefficient * b
        else:
            lo += coefficient * b
            hi += coefficient * a
    return rounded((lo, hi))


@lru_cache(maxsize=None)
def positive_integral(m: int, block: int, e: int, T: F, M: int) -> Interval:
    """Integral p^(m-1) B_e^2(log(p/p_b)/Delta) dp on one block.

    Substitution p=p_b*exp(Delta*u) gives
      Delta*p_b^m * sum_{r>=0} (m*Delta)^r/r! *
                        2*(r+e)! / (e!*(r+3)!).
    Every term is nonnegative. The Beta factor is at most 1/3.
    """
    require(m >= 1 and 0 <= e <= 2, 'Unsupported integral')
    delta = T / M
    a = m * delta
    terms = 80
    require(a < terms + 2, 'Invalid integral tail bound')
    term = F(1)
    total = F(0)
    for r in range(terms + 1):
        if e == 0:
            beta = F(2, (r + 1) * (r + 2) * (r + 3))
        elif e == 1:
            beta = F(2, (r + 2) * (r + 3))
        else:
            beta = F(1, r + 3)
        total += term * beta
        if r < terms:
            term *= a / (r + 1)
    remainder = term * a / (terms + 1) / (1 - a / (terms + 2)) / 3
    p_lo, p_hi = exp_interval(m * (-T + block * delta))
    return rounded((delta * p_lo * total,
                    delta * p_hi * (total + remainder)))


@lru_cache(maxsize=None)
def weight(d: int, n: int, block: int, e: int, T: F, M: int) -> Interval:
    require(0 <= n <= d, 'Invalid binomial exponents')
    return signed_sum([
        (F((-1)**v * comb(d - n, v)),
         positive_integral(n + v + 1, block, e, T, M))
        for v in range(d - n + 1)
    ])


def decimal(x: F, precision: int = 55) -> str:
    with localcontext() as context:
        context.prec = precision
        return str(Decimal(x.numerator) / Decimal(x.denominator))


def audit(path: Path) -> dict:
    started = time.perf_counter()
    with gzip.open(path, 'rt', encoding='utf-8') as stream:
        data = json.load(stream)
    s = data['summary']
    K, h, M, Q = (s[k] for k in ('K', 'h', 'M', 'denominator'))
    T = F(s['T'])
    require((K, h, M, T, s['factor']) == (6, 2, 16, F(17, 4), 4),
            'This independent audit is for the selected certificate')
    require(s['variance_bound'] == 'exact' and s['tight_B'] == 1,
            'Sharp moments and actual-large-size bound are required')
    require(s['schema'] == 'finecap-moment-cross2-strict-v1', 'Unexpected schema')
    arrays = {name: {freeze(k): F(v, Q) for k, v in data[name]}
              for name in ('P', 'minor', 'slope', 'curve')}
    P, U, V, C = (arrays[name] for name in ('P', 'minor', 'slope', 'curve'))
    require(all(x <= 0 for x in V.values()), 'A slope is positive')
    require(all(x <= 0 for x in C.values()), 'A curvature is positive')
    reserves = tuple(F(x) for x in data['correction_slopes'])
    require(len(reserves) == M and all(x >= 0 for x in reserves), 'Bad reserves')
    hstar = T / M * sum(reserves)
    alpha = 1 / (1 + hstar)
    require(hstar == F(s['correction_Hstar']) and alpha == F(s['alpha']),
            'Reserve normalization mismatch')
    delta = F(h, K)

    def jmax(mass: F) -> F:
        require(mass >= 0, 'Negative remaining small mass')
        count = mass // delta
        return count * delta**2 + (mass - count * delta)**2

    classes = tuple(range(h + 1, K + 1))
    patterns = sorted(f for f in product(*[range(1 + (K - 1) // (j - 1))
                                          for j in classes])
                      if sum((j - 1) * x for j, x in zip(classes, f)) <= K - 1)
    results = []
    for full in patterns:
        lower = sum((j - 1) * x for j, x in zip(classes, full))
        B0 = F(K - lower, K)
        for j in range(1, K + 1):
            eligible = full[j - h - 1] > 0 if j > h else lower + j <= K
            if not eligible:
                continue
            companions = [cls for cls, count in zip(classes, full)
                          for _ in range(count)]
            if j > h:
                companions.remove(j)
            d = len(companions)
            for endpoint in (0, 1):
                size = F(j - 1 + endpoint, K)
                all_terms: list[tuple[F, Interval]] = []
                for flags in product((0, 1), repeat=d):
                    observed = [0] * len(classes)
                    if j > h:
                        observed[j - h - 1] += 1
                    for cls, present in zip(companions, flags):
                        observed[cls - h - 1] += present
                    obs = tuple(observed)
                    n = sum(flags)
                    for b in range(M):
                        for e in range(3):
                            key = (full, obs, j, b, e)
                            u, v, c = U[key + (endpoint,)], V[key], C[key]
                            if j > h:
                                mass = B0 - F(endpoint, K)
                                coefficients = (u, v * mass + c * jmax(mass),
                                                c * (mass**2 - jmax(mass)))
                            else:
                                mass = B0 - size
                                coefficients = (u + v * size + c * size**2,
                                                v * mass + c * (2 * size * mass + jmax(mass)),
                                                c * (mass**2 - jmax(mass)))
                            for shift, coeff in enumerate(coefficients):
                                all_terms.append((coeff, weight(d + shift, n + shift,
                                                               b, e, T, M)))
                lo, hi = signed_sum(all_terms)
                lo, hi = rounded((alpha * lo, alpha * hi))
                results.append({'pattern': full, 'class': j, 'endpoint': endpoint,
                                'lower': str(lo), 'upper': str(hi),
                                'lower_decimal': decimal(lo)})
    require(len(patterns) == 7 and len(results) == 38, 'Wrong service family size')
    claimed = F(s['limiting_coefficient'])
    minimum = min(F(row['lower']) for row in results)
    require(minimum >= claimed, 'Independent service lower bound is insufficient')
    require(claimed <= 1 - exp_interval(-T)[1], 'Zero-size service is insufficient')
    curvature = max(F(6 * M * M, 1) / T**2 *
                    abs(P[3*b+r+2, z] - 2*P[3*b+r+1, z] + P[3*b+r, z])
                    for b in range(M) for r in range(2) for z in range(K+1))
    reported_curvature = F(s['second_derivative_bound'])
    require(reported_curvature >= curvature, 'Curvature bound is insufficient')
    witness_norm = max([F(1)] + [abs(value) + abs(V[key[:-1]]) + abs(C[key[:-1]])
                                      for key, value in U.items()])
    exp_hi = exp_interval(T)[1]
    E = exp_hi.numerator // exp_hi.denominator + 1
    dmax = (K - 1) // h
    finite_required = witness_norm * (4*(dmax+2)**2+(20*E+20)*M+20) + F(E, 2)*T*reported_curvature
    finite_A = F(s['finite_loss_constant'])
    require(finite_A >= finite_required, 'Finite loss constant is insufficient')
    N = 10**14
    eps = F(1, N)
    target = F('0.25647580')
    finite = claimed - finite_A/N - eps
    require(finite > target, 'Requested finite guarantee fails')
    return {'status': 'PASS', 'certificate_sha256': hashlib.sha256(path.read_bytes()).hexdigest(),
            'method': 'Individual companion subsets; positive Taylor/Beta integrals; signed rational interval sums',
            'full_patterns': len(patterns), 'service_coefficients': len(results),
            'independent_minimum_lower': str(minimum),
            'independent_minimum_decimal': decimal(minimum),
            'claimed_limiting_lower': str(claimed), 'claimed_limiting_decimal': decimal(claimed),
            'witness_norm': str(witness_norm), 'curvature_bound': str(reported_curvature),
            'E': E, 'finite_loss_required': str(finite_required), 'finite_loss_A': str(finite_A),
            'N': N, 'sampling_failure': str(eps), 'target': str(target),
            'finite_coefficient': str(finite), 'finite_decimal': decimal(finite),
            'margin': str(finite-target), 'margin_decimal': decimal(finite-target),
            'coefficients': results, 'elapsed_seconds': round(time.perf_counter()-started, 3)}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('certificate', nargs='?', type=Path, default=DEFAULT)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    try:
        result = audit(args.certificate)
        if args.output:
            args.output.parent.mkdir(parents=True, exist_ok=True)
            args.output.write_text(json.dumps(result, indent=2)+'\n', encoding='utf-8')
        print(json.dumps({k:v for k,v in result.items() if k != 'coefficients'}, indent=2))
    except (ValueError, KeyError, OSError, ZeroDivisionError) as error:
        print(f'FAIL: {error}', file=sys.stderr)
        raise SystemExit(1) from error


if __name__ == '__main__':
    main()
