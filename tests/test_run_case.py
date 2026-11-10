"""Offline tests for the container entry point's tracer pre-check (docker/run_case.py).

The expected counts are the values SAM itself printed ("SAM-ASP tracer layout") for these decks.
Run:  python3 -m unittest discover -s tests -v
"""
import io
import os
import re
import shutil
import sys
import tempfile
import unittest
from contextlib import redirect_stderr, redirect_stdout

HERE = os.path.dirname(os.path.abspath(__file__))
TOP = os.path.dirname(HERE)
os.environ['SAMASP_HOME'] = TOP
sys.path.insert(0, os.path.join(TOP, 'docker'))
import run_case                                     # noqa: E402

WF = os.path.join(TOP, 'cases', 'williams_flats_2019', 'decks')
W09 = os.path.join(TOP, 'cases', 'williams_2009', 'decks')


def set_bins(deck_dir, nbins):
    """change the number of size bins in AerosolModes.in (keeps CRLF)"""
    p = os.path.join(deck_dir, 'AerosolModes.in')
    raw = open(p, 'rb').read().decode('latin-1')
    raw = re.sub(r'(?m)^\s*\d+(\s*!!\s*Number of Bins)', r'%d\1' % nbins, raw, count=1)
    open(p, 'wb').write(raw.encode('latin-1'))


class TracerCount(unittest.TestCase):
    def test_counts_match_sam(self):
        self.assertEqual(run_case.deck_tracers(WF)[0], 1488)
        self.assertEqual(run_case.deck_tracers(W09)[0], 1486)

    def test_forty_bins(self):
        d = tempfile.mkdtemp()
        try:
            dd = os.path.join(d, 'decks')
            shutil.copytree(W09, dd)
            set_bins(dd, 40)
            self.assertEqual(run_case.deck_tracers(dd)[0], 4096)
        finally:
            shutil.rmtree(d)

    def test_compiled_value(self):
        self.assertIsInstance(run_case.compiled_ntracers(), int)

    def test_mismatch_stops_before_work(self):
        have = run_case.compiled_ntracers()
        deck = WF if run_case.deck_tracers(WF)[0] != have else W09
        err = io.StringIO()
        with redirect_stderr(err), redirect_stdout(io.StringIO()):
            with self.assertRaises(SystemExit) as cm:
                run_case.main(['--check-decks', '--decks', deck])
        self.assertEqual(cm.exception.code, 2)
        self.assertIn('NTRACERS=%d' % run_case.deck_tracers(deck)[0], err.getvalue())

    def test_match_passes(self):
        have = run_case.compiled_ntracers()
        deck = WF if run_case.deck_tracers(WF)[0] == have else W09
        with redirect_stdout(io.StringIO()):
            self.assertEqual(run_case.main(['--check-decks', '--decks', deck]), 0)


if __name__ == '__main__':
    unittest.main()
