"""Fail-closed CLI controls; compile only inspected source, never bundled binaries."""
import unittest,subprocess,tempfile,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
class CLI(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp=tempfile.TemporaryDirectory();cls.binary=Path(cls.tmp.name)/'splitrel'
        subprocess.run(['cc','-O2','-Wall','-Wextra','-o',str(cls.binary),str(ROOT/'splitrel.c')],check=True,capture_output=True)
    @classmethod
    def tearDownClass(cls):cls.tmp.cleanup()
    def run_c(self,data,args=()):return subprocess.run([str(self.binary),*args],input=data,capture_output=True)
    def test_k3(self):
        r=self.run_c(b'Bw\n',('--conn',));self.assertEqual(r.returncode,0);self.assertEqual(r.stdout,b'Bw 3 0 0 3 1\n')
    def test_invalid(self):
        for data in (b'',b'z\n',b'K~~~~~~~~~~\n',b'B!\n',b'Bx\n',b'Bwextra\n',b'>>graph6<<Bw\n',b'Bw\x00junk\n',b'B'+b'w'*4096+b'\n'):
            with self.subTest(data=data):self.assertNotEqual(self.run_c(data,('--conn',)).returncode,0)
    def test_unknown_flag(self):self.assertNotEqual(self.run_c(b'Bw\n',('--unknown',)).returncode,0)
    def test_mixed_maxonly(self):self.assertNotEqual(self.run_c(b'Bw\nA_\n',('--maxonly',)).returncode,0)
    def test_injected_failures(self):
        for flag in ('--inject-failure','--inject-coverage-gap'):
            r=subprocess.run([sys.executable,'-O','-B',str(ROOT/'verification/replay_cones.py'),'--output-dir',self.tmp.name,'--max-N','7',flag],capture_output=True)
            self.assertNotEqual(r.returncode,0);self.assertIn(b'RuntimeError',r.stderr)
    def test_empty_cone_range(self):
        r=subprocess.run([sys.executable,'-O','-B',str(ROOT/'cone_star2.py'),'50','7'],capture_output=True);self.assertNotEqual(r.returncode,0)
if __name__=='__main__':unittest.main()
