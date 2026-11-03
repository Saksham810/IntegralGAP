#!/usr/bin/env python3
"""Run all three integral verification checks, using Python >= 3.10 only.

Usage: python -S verify_all.py [--output NEW_DIRECTORY]
No third-party dependencies, optimizer, network access, or online simulation.
The supplied mathematical programs and rational table are preserved unchanged.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
from decimal import Decimal, localcontext
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import platform
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parent
INTEGRAL = Path('verification/integral')
RESEARCH = INTEGRAL / 'research'
TABLE = RESEARCH / 'runs/eplus1_search_expanded5_cert.json.gz'
TABLE_SHA256 = 'ed6fe9bcdcb3f24cacbc6b5a02cb823fb6f57bb124066c1e7d7e5f7d6ab4d264'
GAMMA = F(307770970018984529425873353861053278824139,
          1200000005486015625000000000000000000000000)
LOSS = F(4718849996128173, 17000000000)
INPUT_FILES = (
    (RESEARCH / 'eplus1_audit_true_spatial.py').as_posix(),
    (RESEARCH / 'exact_geometry.py').as_posix(),
    (RESEARCH / 'full_geometry.py').as_posix(),
    (RESEARCH / 'push027_cross_geometry.py').as_posix(),
    (RESEARCH / 'push027_cross_strict_replay.py').as_posix(),
    (RESEARCH / 'push027_strict_geometry.py').as_posix(),
    (RESEARCH / 'replay_exact.py').as_posix(),
    TABLE.as_posix(),
    (INTEGRAL / 'verify_service_independent.py').as_posix(),
    'verify_all.py',
)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def write_json(path: Path, value: dict) -> None:
    path.write_text(json.dumps(value, indent=2, ensure_ascii=True) + '\n',
                    encoding='utf-8')


def decimal(value: F) -> str:
    # Decimal strings are for display; every comparison below uses Fraction.
    with localcontext() as context:
        context.prec = 70
        return str(Decimal(value.numerator) / Decimal(value.denominator))


def integrity() -> int:
    manifest = json.loads((ROOT / 'SHA256SUMS.json').read_text(encoding='utf-8'))
    require(isinstance(manifest, dict), 'Invalid SHA256SUMS.json')
    require(set(manifest) == set(INPUT_FILES),
            'Hash manifest must list exactly the 10 required files')
    for name, expected in manifest.items():
        path = (ROOT / name).resolve()
        require(path.is_relative_to(ROOT), 'Unsafe manifest path: ' + name)
        require(path.is_file(), 'Missing required file: ' + name)
        actual = hashlib.sha256(path.read_bytes()).hexdigest()
        require(actual == expected, 'SHA-256 mismatch: ' + name)
    require(hashlib.sha256((ROOT / TABLE).read_bytes()).hexdigest() == TABLE_SHA256,
            'Unexpected rational table')
    return len(manifest)


def run_check(name: str, arguments: list[str], output: Path,
              *, result_from_stdout: bool = False) -> dict:
    print('[' + name + '] running ...', flush=True)
    command = [sys.executable, '-u', '-B', '-S', *arguments]
    # Child modules import local geometry from the directory of their script.
    # Use relative input paths so the saved mathematical output is relocatable.
    with (output / (name + '.stdout.log')).open('w', encoding='utf-8') as stdout, \
         (output / (name + '.stderr.log')).open('w', encoding='utf-8') as stderr:
        completed = subprocess.run(command, cwd=ROOT, stdout=stdout, stderr=stderr,
                                   check=False)
    require(completed.returncode == 0,
            name + ' failed (exit ' + str(completed.returncode) +
            '); inspect its stdout/stderr logs in ' + str(output))
    path = output / (name + '.json')
    if result_from_stdout:
        lines = (output / (name + '.stdout.log')).read_text(encoding='utf-8').splitlines()
        require(bool(lines), name + ' produced no result')
        result = json.loads(lines[-1])
        write_json(path, result)
    else:
        result = json.loads(path.read_text(encoding='utf-8'))
    require(isinstance(result, dict) and result.get('status') == 'PASS',
            name + ' did not report PASS')
    print('[' + name + '] PASS', flush=True)
    return result


def verify(output: Path) -> dict:
    started = time.perf_counter()
    count = integrity()
    print('Integrity: ' + str(count) + ' required files match SHA-256.', flush=True)
    full = run_check('integral_full',
                     [str(RESEARCH / 'push027_cross_strict_replay.py'), str(TABLE)],
                     output, result_from_stdout=True)
    spatial = run_check('integral_spatial',
                        [str(RESEARCH / 'eplus1_audit_true_spatial.py'), str(TABLE),
                         str(output / 'integral_spatial.json')], output)
    service = run_check('integral_service',
                        [str(INTEGRAL / 'verify_service_independent.py'), str(TABLE),
                         '--output', str(output / 'integral_service.json')], output)

    require(F(full['limiting_lower_bound']) == GAMMA,
            'Unexpected limiting comparison coefficient')
    require(F(full['finite_loss_constant']) == LOSS,
            'Unexpected finite-horizon loss constant')
    expected_counts = {
        'counts': 7, 'labels': 146, 'rules': 290,
        'availability_checks': 226688, 'crossing_checks': 723840,
        'packing_checks': 219072, 'minorant_checks': 47760,
        'service_coefficients': 38,
    }
    for key, expected in expected_counts.items():
        require(full.get(key) == expected, 'Unexpected full replay count: ' + key)
    require(spatial['certificate_sha256'] == TABLE_SHA256,
            'Spatial audit used a different rational table')
    require(spatial['exact_corner_and_quadratic_checks'] == 473280 and
            spatial['additional_interior_polynomial_identity_checks'] == 125280 and
            spatial['jump_endpoint_occurrences'] == 1914,
            'Unexpected spatial audit coverage')
    require(F(spatial['minimum_margin_unscaled']) >= 0,
            'Negative spatial audit margin')
    require(service['certificate_sha256'] == TABLE_SHA256,
            'Independent value audit used a different rational table')
    require(service['full_patterns'] == 7 and
            service['service_coefficients'] == 38 and
            len(service['coefficients']) == 38,
            'Unexpected independent value audit coverage')
    require(F(service['claimed_limiting_lower']) == GAMMA and
            F(service['independent_minimum_lower']) >= GAMMA,
            'Independent comparison bound is insufficient')
    require(F(service['finite_loss_A']) == LOSS and
            F(service['finite_loss_required']) <= LOSS,
            'Independent finite loss exceeds the claim')

    guarantees = []
    for horizon, failure, target in (
        (10**10, F(1, 10**8), F('0.2564')),
        (10**14, F(1, 10**14), F('0.25647580')),
    ):
        lower = GAMMA - LOSS / horizon - failure
        require(lower > target, 'Finite guarantee fails at N=' + str(horizon))
        guarantees.append({
            'N_min': horizon, 'sampling_failure': str(failure),
            'target': str(target), 'lower_coefficient': str(lower),
            'decimal': decimal(lower), 'strict_margin': str(lower - target),
        })
    require(F('0.2564') == F(641, 2500), 'Ratio conversion failed')
    require(F(2500, 641) < F('3.9002'), 'Competitive-factor rounding failed')
    # Detect a file edit during the run as well as edits made beforehand.
    integrity()
    return {
        'status': 'PASS', 'scope': 'integral verification only',
        'python': platform.python_version(), 'platform': platform.platform(),
        'manifest_files': count, 'integral_table_sha256': TABLE_SHA256,
        'integral_full': full, 'integral_spatial': spatial,
        'integral_service': {k: v for k, v in service.items() if k != 'coefficients'},
        'finite_integral_guarantees': guarantees,
        'displayed_reward_ratio': '0.2564',
        'competitive_factor_exact': '2500/641',
        'competitive_factor_decimal': decimal(F(2500, 641)),
        'elapsed_seconds': round(time.perf_counter() - started, 3),
        'proof_scope': 'Finite numerical hypotheses and finite-horizon arithmetic. '
                       'The manuscript supplies the reductions to continuous '
                       'sizes, online allocation, and integral GPAP.',
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=ROOT / 'results',
                        help='new or empty results directory (default: results beside this script)')
    args = parser.parse_args()
    output = None
    try:
        require(sys.version_info >= (3, 10), 'Use Python 3.10 or newer')
        require(sys.flags.optimize == 0,
                'Do not use -O/-OO or a nonzero PYTHONOPTIMIZE: '
                'the spatial audit uses assertions')
        destination = args.output.resolve()
        require(not destination.exists() or
                (destination.is_dir() and not any(destination.iterdir())),
                'Output directory must be new or empty. '
                'Choose a fresh directory with --output.')
        destination.mkdir(parents=True, exist_ok=True)
        output = destination
        started_utc = datetime.now(timezone.utc).isoformat()
        write_json(output / 'run_status.json',
                   {'status': 'RUNNING', 'started_utc': started_utc})
        result = verify(output)
        result['started_utc'] = started_utc
        result['completed_utc'] = datetime.now(timezone.utc).isoformat()
        write_json(output / 'all_checks.json', result)
        write_json(output / 'run_status.json',
                   {'status': 'PASS', 'started_utc': started_utc,
                    'completed_utc': result['completed_utc']})
        print('ALL INTEGRAL CHECKS PASS: finite ratio > 0.2564.', flush=True)
        print('Stronger finite check also passes: ratio > 0.25647580.', flush=True)
        print('Results: ' + str(output / 'all_checks.json'), flush=True)
        return 0
    except (OSError, ValueError, KeyError, TypeError, ZeroDivisionError,
            subprocess.SubprocessError, KeyboardInterrupt) as error:
        if output is not None:
            write_json(output / 'run_status.json',
                       {'status': 'FAIL', 'error': str(error) or type(error).__name__})
        print('FAIL: ' + (str(error) or type(error).__name__), file=sys.stderr)
        return 1


if __name__ == '__main__':
    raise SystemExit(main())
