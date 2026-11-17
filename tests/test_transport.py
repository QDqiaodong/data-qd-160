import copy,json,subprocess,tempfile,unittest,sys
from pathlib import Path
from reaction_batch.model import analyze
from reaction_batch.cli import run
ROOT=Path(__file__).resolve().parents[1]
def sample():return json.loads((ROOT/'fixtures/input.json').read_text())
class Transport(unittest.TestCase):
 def test_installed_fixture(self):self.assertEqual(run(ROOT/'fixtures/input.json')['schema_version'],1)
 def test_unknown_field(self):
  d=sample();d['extra']=1
  with self.assertRaises(ValueError):analyze(d)
 def test_cli_preserves_output(self):
  with tempfile.TemporaryDirectory() as tmp:
   p=Path(tmp);(p/'in.json').write_text('{}');(p/'out.json').write_text('old')
   c=subprocess.run([sys.executable,'-m','reaction_batch.cli','--source',str(p/'in.json'),'--output',str(p/'out.json')],capture_output=True)
   self.assertEqual(c.returncode,2);self.assertEqual((p/'out.json').read_text(),'old')
 def test_reject_duplicate_json(self):
  with tempfile.TemporaryDirectory() as tmp:
   p=Path(tmp)/'bad.json';p.write_text('{"a":1,"a":2}')
   with self.assertRaises(ValueError):run(p)
 def test_output_path_guards(self):
  with tempfile.TemporaryDirectory() as tmp:
   p=Path(tmp);src=p/'in.json';src.write_text(json.dumps(sample()));old=src.read_bytes();out=p/'link';out.symlink_to(src)
   for dest in [src,out]:
    c=subprocess.run([sys.executable,'-m','reaction_batch.cli','--source',str(src),'--output',str(dest)],capture_output=True)
    self.assertEqual(c.returncode,2);self.assertEqual(src.read_bytes(),old)
