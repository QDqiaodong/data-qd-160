import unittest,json
from fractions import Fraction
from pathlib import Path
from reaction_batch.model import analyze
from reaction_batch.chemistry import formula,balance
ROOT=Path(__file__).resolve().parents[1]
def fixture():return json.loads((ROOT/'fixtures/input.json').read_text())
class Reactions(unittest.TestCase):
 def test_water_batch(self):
  r=analyze(fixture());self.assertEqual(r['coefficients'],{'H2':2,'O2':1,'H2O':2});self.assertEqual(r['limiting_reactants'],['O2']);self.assertEqual(r['remaining_moles'],{'H2':'1','O2':'0'});self.assertEqual(r['produced_moles'],{'H2O':'4'})
 def test_parentheses(self):self.assertEqual(formula('Al2(SO4)3'),{'Al':2,'S':3,'O':12});self.assertEqual(formula('K4(ON(SO3)2)2'),{'K':4,'O':14,'N':2,'S':4})
 def test_iron_balance(self):self.assertEqual(balance(['Fe','O2'],['Fe2O3'])[0],[4,3,2])
 def test_combustion(self):self.assertEqual(balance(['C2H6','O2'],['CO2','H2O'])[0],[2,7,4,6])
 def test_underdetermined(self):
  with self.assertRaisesRegex(ValueError,'one-dimensional'):balance(['C','O2'],['CO','CO2'])
 def test_impossible(self):
  with self.assertRaises(ValueError):balance(['H2'],['CO2'])
 def test_invalid_formulas(self):
  for text in ['', 'H0','H01','H2O)','()','Mg(OH','CuSO4.5H2O','2H2','Na+']:
   with self.subTest(text=text),self.assertRaises(ValueError):formula(text)
 def test_colimiting_and_zero(self):
  d=fixture();d['available_moles']={'H2':'4','O2':'2'};self.assertEqual(analyze(d)['limiting_reactants'],['H2','O2']);d['available_moles']['H2']='0';self.assertEqual(analyze(d)['extent_moles'],'0')
 def test_exact_fraction(self):
  d=fixture();d['available_moles']={'H2':'0.3','O2':'0.1'};self.assertEqual(analyze(d)['produced_moles']['H2O'],'1/5');self.assertEqual(analyze(d)['remaining_moles']['H2'],'1/10')
 def test_atom_conservation_from_output(self):
  d={'reactants':['Al','HCl'],'products':['AlCl3','H2'],'available_moles':{'Al':'3','HCl':'7'}};r=analyze(d)
  for element in ['Al','H','Cl']:
   before=sum(formula(f).get(element,0)*Fraction(q) for f,q in r['consumed_moles'].items());after=sum(formula(f).get(element,0)*Fraction(q) for f,q in r['produced_moles'].items());self.assertEqual(before,after)
