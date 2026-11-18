"""Offline tests for the container entry point (docker/run_case.py): the tracer pre-check
and the case settings.

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

    def check_decks(self, deck, ntracers):
        """run --check-decks with the model's ntracers fixed; -> (exit code, stdout, stderr)"""
        saved = run_case.compiled_ntracers
        run_case.compiled_ntracers = lambda: ntracers
        out, err = io.StringIO(), io.StringIO()
        try:
            with redirect_stderr(err), redirect_stdout(out):
                try:
                    code = run_case.main(['--check-decks', '--decks', deck])
                except SystemExit as e:
                    code = e.code
        finally:
            run_case.compiled_ntracers = saved
        return code, out.getvalue(), err.getvalue()

    def test_mismatch_stops_before_work(self):
        code, _, err = self.check_decks(W09, 1488)
        self.assertEqual(code, 2)
        self.assertIn('NTRACERS=1486', err)

    def test_match_passes(self):
        code, out, _ = self.check_decks(W09, 1486)
        self.assertEqual(code, 0)
        self.assertIn('tracers: 1486', out)

    def test_missing_deck_file(self):
        d = tempfile.mkdtemp()
        try:
            err = io.StringIO()
            with redirect_stderr(err), self.assertRaises(SystemExit):
                run_case.require_decks(d)
            self.assertIn('ASPInputDeck.in', err.getvalue())
        finally:
            shutil.rmtree(d)


TEMPLATE = os.path.join(TOP, 'RUNDATA', 'prm_template')


class CaseSettings(unittest.TestCase):
    def args(self, *argv):
        return run_case.parse_args(['--decks', W09] + list(argv))

    def test_template_fire_start_is_model_default(self):
        """the template sets fire_nstep_start to the value params.f90 would use anyway"""
        src = open(os.path.join(TOP, 'SRC_TUV', 'params.f90')).read()
        default = int(re.search(r'(?m)^integer:: *fire_nstep_start *= *(\d+)', src).group(1))
        self.assertEqual(run_case.prm_value(open(TEMPLATE).read(), 'fire_nstep_start'), default)

    def test_prm_value(self):
        prm = " dt = 8.,\n ! vg = 3.,\n vg = 2.5d0,  ! wind\n nstatfrq = 5,\n"
        self.assertEqual(run_case.prm_value(prm, 'vg'), 2.5)
        self.assertIsNone(run_case.prm_value(prm, 'nstat'))

    def test_short_run_lowers_nstat(self):
        nstat = run_case.prm_value(open(TEMPLATE).read(), 'nstat')
        sets = run_case.case_settings(self.args(), '/decks', '/out', int(nstat) - 1, TEMPLATE)
        self.assertIn('nstat=%d' % (nstat - 1), sets)
        sets = run_case.case_settings(self.args(), '/decks', '/out', int(nstat), TEMPLATE)
        self.assertFalse(any(s.startswith('nstat=') for s in sets))

    def test_settings_order(self):
        sets = run_case.case_settings(self.args('--fire-start-step', '7', '--set', 'perturb_type=2'),
                                      '/decks/', '/out', 900, TEMPLATE)
        self.assertEqual(sets, ["asp_deck_dir='/decks/'", "asp_output_dir='/out/ASP/'",
                                'fire_nstep_start=7', 'perturb_type=2'])

    def test_fire_window(self):
        prm = ' dy = 500.,\n vg = 2.,\n plume_y = 1000.,\n fire_nstep_start = 100,\n'
        with redirect_stdout(io.StringIO()) as out:
            self.assertEqual(run_case.report_fire_window(prm, 1000, 5.0, 4), (100, 200))
        self.assertNotIn('WARNING', out.getvalue())
        with redirect_stdout(io.StringIO()) as out:
            run_case.report_fire_window(prm, 150, 5.0, 4)
        self.assertIn('before the fire emissions end', out.getvalue())
        with redirect_stdout(io.StringIO()) as out:
            run_case.report_fire_window(prm, 100, 5.0, 4)
        self.assertIn('before the fire starts', out.getvalue())

    def test_fire_window_needs_fire_start(self):
        with redirect_stderr(io.StringIO()), self.assertRaises(SystemExit):
            run_case.report_fire_window(' dy = 500.,\n vg = 2.,\n plume_y = 1000.,\n', 1000, 5.0, 4)

if __name__ == '__main__':
    unittest.main()
