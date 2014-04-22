#!/usr/bin/env python
"""
Basic unit tests for balancer functions. Not very complete.

Liz Wylie
February 2014
"""
# TODO: Refactor code to make better/more complete/faster unittests.

import unittest
from numpy import array
import json
from collections import Counter
from balancer import balance, check_balance, count, find_decomposition, candidate, get_atom_info, get_Hs, update
from rdkit import Chem
from extract_balance_rxns import check_entry, generate_data


class BalanceTest(unittest.TestCase):
    def test(self):
        self.assertEqual(balance('CCc1ccncc1>>CCC1CCNCC1.CCC1=CCNCC1', False), 'CCc1ccncc1.CCc1ccncc1>>CCC1CCNCC1.CCC1=CCNCC1')
        self.assertEqual(balance('CCc1ccncc1>>CCC1CCNCC1.CCC1=CCNCC1', True), None)
        self.assertEqual(balance('O=Cc1ccsc1>>O=C(c1ccsc1)C(O)c1ccsc1', False), 'O=Cc1ccsc1.O=Cc1ccsc1>>O=C(c1ccsc1)C(O)c1ccsc1')
        self.assertEqual(balance('O=Cc1ccsc1>>O=C(c1ccsc1)C(O)c1ccsc1', True), 'O=Cc1ccsc1.O=Cc1ccsc1>>O=C(c1ccsc1)C(O)c1ccsc1')
        self.assertEqual(balance('C1OCOCO1>>C=O', False), 'C1OCOCO1>>C=O.C=O.C=O')
        self.assertEqual(balance('C1OCOCO1>>C=O', True), 'C1OCOCO1>>C=O.C=O.C=O')
        self.assertEqual(balance("CC1CCCCC1CO.O=O>>CC1CCCCC1C(=O)O", False), None)
        self.assertEqual(balance("CC1CCCCC1CO.O=O>>CC1CCCCC1C(=O)O", True), None)


class CheckBalanceTest(unittest.TestCase):
    def test(self):
        self.assertEqual(check_balance('CCc1ccncc1>>CCC1CCNCC1.CCC1=CCNCC1', False), False)
        self.assertEqual(check_balance('CCc1ccncc1.CCc1ccncc1>>CCC1CCNCC1.CCC1=CCNCC1', False), True)
        self.assertEqual(check_balance('CCc1ccncc1.CCc1ccncc1>>CCC1CCNCC1.CCC1=CCNCC1', True), False)
        self.assertEqual(check_balance('O=Cc1ccsc1.O=Cc1ccsc1>>O=C(c1ccsc1)C(O)c1ccsc1', False), True)
        self.assertEqual(check_balance('O=Cc1ccsc1.O=Cc1ccsc1>>O=C(c1ccsc1)C(O)c1ccsc1', True), True)
        self.assertEqual(check_balance('C1OCOCO1>>C=O.C=O.C=O', True), True)
        self.assertEqual(check_balance('C1OCOCO1>>C=O', True), False)


class GetHsTest(unittest.TestCase):
    def setUp(self):
        with open("chemical_test.json") as f:
            self.data = json.load(f)

    def test(self):
        for k, v in self.data.items():
            self.assertEqual(get_Hs(Chem.MolFromSmiles(k.encode('ascii'))), v["get_Hs"])


class GetAtomInfoTest(unittest.TestCase):
    def setUp(self):
        with open("chemical_test.json") as f:
            self.data = json.load(f)

    def test(self):
        for k, v in self.data.items():
            self.assertEqual(get_atom_info(Chem.MolFromSmiles(k.encode('ascii'))), v["get_atom_info"])


class CountTest(unittest.TestCase):
    def test(self):
        self.assertEqual(count(['C', 'C', 'C', 'C', 'C', 'N', 'C', 'C'], 9), Counter({'H': 9, 'C': 7, 'N': 1}))
        self.assertEqual(count(['C', 'C', 'C', 'C', 'C', 'N', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'N', 'C', 'C'], 28), Counter({'H': 28, 'C': 14, 'N': 2}))
        self.assertEqual(count(['C', 'C', 'C', 'C', 'C', 'N', 'C', 'C'], 13), Counter({'H': 13, 'C': 7, 'N': 1}))
        self.assertEqual(count(['C', 'C', 'C', 'C', 'C', 'N', 'C', 'C']), Counter({'C': 7, 'N': 1}))
        self.assertEqual(count(['C', 'C', 'C', 'C', 'C', 'N', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'N', 'C', 'C']), Counter({'C': 14, 'N': 2}))


class CandidateTest(unittest.TestCase):
    def test(self):
        self.assertEqual(candidate(Counter({'C': 7, 'N': 1}), Counter({'C': 7, 'N': 1})), True)
        self.assertEqual(candidate(Counter({'H': 9, 'C': 7, 'N': 1}), Counter({'H': 19, 'C': 7, 'N': 1})), True)
        self.assertEqual(candidate(Counter({'H': 2, 'C': 1, 'O': 1}), Counter({'H': 4, 'C': 2, 'O': 2})), True)
        self.assertEqual(candidate(Counter({'H': 8, 'C': 5, 'O': 10}), Counter({'H': 4, 'C': 2, 'O': 2})), False)


class FindDecompTest(unittest.TestCase):
    def test(self):
        self.assertEqual(find_decomposition([Counter({'C': 7, 'N': 1})], Counter({'C': 7, 'N': 1})), [(Counter({'C': 7, 'N': 1}), array([ 1.]))])
        self.assertEqual(find_decomposition([Counter({'C': 1, 'O': 1})], Counter({'C': 2, 'O': 2})), [(Counter({'C': 1, 'O': 1}), array([ 2.]))])


class UpdateTest(unittest.TestCase):
    def test(self):
        self.assertEqual(update([(Counter({'C': 7, 'N': 1}), array([ 1.]))], [Counter({'C': 7, 'N': 1})], ['CCc1ccncc1']), ['CCc1ccncc1', 'CCc1ccncc1'])


class CheckEntryTest(unittest.TestCase):
    def setUp(self):
        with open("entry_test.json") as f:
            self.data = json.load(f)

    def test(self):
        for v in self.data.values():
            self.assertEqual(check_entry(v["entry"]), v["check_entry"])


class GenerateDataTest(unittest.TestCase):
    def setUp(self):
        with open("entry_test.json") as f:
            self.data = json.load(f)

    def test(self):
        for v in self.data.values():
            self.assertEqual(generate_data(v["entry"]), v["generate_data"])

if __name__ == '__main__':
    unittest.main()