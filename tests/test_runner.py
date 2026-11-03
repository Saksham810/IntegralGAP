"""Failure-path checks for the integral-only orchestrator.

These tests run no mathematical replay. A full verification is a separate step.
Each test uses a private copy and never modifies the published source files.
"""
from __future__ import annotations

import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
TABLE = Path('verification/integral/research/runs/eplus1_search_expanded5_cert.json.gz')
CHECKER = Path('verification/integral/research/push027_cross_strict_replay.py')


class RunnerFailureTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temp = tempfile.TemporaryDirectory(prefix='gpap-runner-test-')
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name) / 'copy'
        self.root.mkdir()
        shutil.copy2(ROOT / 'verify_all.py', self.root / 'verify_all.py')
        shutil.copy2(ROOT / 'SHA256SUMS.json', self.root / 'SHA256SUMS.json')
        shutil.copytree(ROOT / 'verification', self.root / 'verification')
        self.output = self.root / 'fresh-results'

    def invoke(self, flags=(), *, env=None, output=None) -> subprocess.CompletedProcess:
        clean_env = os.environ.copy()
        clean_env.pop('PYTHONOPTIMIZE', None)
        if env:
            clean_env.update(env)
        return subprocess.run(
            [sys.executable, '-B', '-S', *flags, str(self.root / 'verify_all.py'),
             '--output', str(output if output is not None else self.output)],
            cwd=self.temp.name, env=clean_env, capture_output=True, text=True,
            timeout=30, check=False,
        )

    def assertFailedBeforeReplay(self, completed, diagnostic: str) -> None:
        self.assertNotEqual(completed.returncode, 0)
        self.assertIn(diagnostic, completed.stderr)
        self.assertNotIn('ALL INTEGRAL CHECKS PASS', completed.stdout)
        self.assertNotIn('[integral_full] running', completed.stdout)
        self.assertFalse((self.output / 'all_checks.json').exists())

    def test_missing_checker_is_rejected(self) -> None:
        (self.root / CHECKER).unlink()
        completed = self.invoke()
        self.assertFailedBeforeReplay(completed, 'Missing required file')
        state = json.loads((self.output / 'run_status.json').read_text())
        self.assertEqual(state['status'], 'FAIL')

    def test_changed_checker_is_rejected(self) -> None:
        with (self.root / CHECKER).open('ab') as stream:
            stream.write(b'\n# changed in a private failure test\n')
        self.assertFailedBeforeReplay(self.invoke(), 'SHA-256 mismatch')

    def test_corrupted_table_is_rejected(self) -> None:
        path = self.root / TABLE
        value = bytearray(path.read_bytes())
        value[len(value) // 2] ^= 1
        path.write_bytes(value)
        self.assertFailedBeforeReplay(self.invoke(), 'SHA-256 mismatch')

    def test_manifest_omission_cannot_skip_a_file(self) -> None:
        path = self.root / 'SHA256SUMS.json'
        manifest = json.loads(path.read_text())
        del manifest[CHECKER.as_posix()]
        path.write_text(json.dumps(manifest))
        self.assertFailedBeforeReplay(self.invoke(), 'exactly the 10 required files')

    def test_optimized_flags_are_rejected(self) -> None:
        for flag in ('-O', '-OO'):
            with self.subTest(flag=flag):
                self.assertFailedBeforeReplay(self.invoke((flag,)), 'Do not use -O/-OO')

    def test_optimization_environment_is_rejected(self) -> None:
        self.assertFailedBeforeReplay(
            self.invoke(env={'PYTHONOPTIMIZE': '1'}), 'Do not use -O/-OO')

    def test_existing_results_are_not_overwritten(self) -> None:
        self.output.mkdir()
        old_summary = '{"status": "PASS", "old_run": true}\n'
        summary = self.output / 'all_checks.json'
        summary.write_text(old_summary)
        completed = self.invoke()
        self.assertNotEqual(completed.returncode, 0)
        self.assertIn('Output directory must be new or empty', completed.stderr)
        self.assertNotIn('ALL INTEGRAL CHECKS PASS', completed.stdout)
        self.assertEqual(summary.read_text(), old_summary)
        self.assertFalse((self.output / 'run_status.json').exists())

    def test_repository_root_cannot_be_used_for_results(self) -> None:
        completed = self.invoke(output=self.root)
        self.assertNotEqual(completed.returncode, 0)
        self.assertIn('Output directory must be new or empty', completed.stderr)
        self.assertFalse((self.root / 'run_status.json').exists())


if __name__ == '__main__':
    unittest.main()
