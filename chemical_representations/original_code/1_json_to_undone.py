#!/usr/bin/env python
"""
Script that performs the first step of creating the desired chemical representation.
This was done in steps... this is primarily for my my own records.
January 2014
Liz Wylie
"""

from rdkit import Chem
import elasticsearch
import json
from os import listdir
from collections import defaultdict


def label(mol):
    for atom in mol.GetAtoms():
        hybridization = int(atom.GetHybridization())
        atom = atom.GetAtomicNum()
        aromaticity = int(atom.GetIsAromatic())
        hashed = int(str(hybridization) + str(atom) + str(aromaticity))
        atom.SetIsotope(hashed)


if __name__ == '__main__':
    es = elasticsearch.Elasticsearch()
    types = set()
    # must be directory containing .json files of list of chemical entries from database as
    # dictionaries. Each entry only needs the _id and SMILES from database.
    files = listdir('tmp')
    for i in files:
        with open('tmp/'+i) as f:
            data = json.load(f)
        for chem in data:
            smi = chem['smiles']
            _id = chem['_id']
            m = Chem.MolFromSmiles(smi.encode('ascii'))
            if not m:
                continue
            DM = Chem.GetDistanceMatrix(m)
            label(m)
            i = {}
            for a in m.GetAtoms():
                i[a.GetIdx()] = a.GetIsotope()
                types.add(a.GetIsotope())
            dat = {
                "status": "undone",
                "body": defaultdict(list)
            }
            for j in i.keys():
                t = defaultdict(list)
                for l, k in enumerate(DM[j]):
                    t[i[l]].append(k)
                dat["body"][i[j]].append(t)
            b = json.dumps(dat)
            es.index(index="dm_data", doc_type="chemical", id=int(_id), body=b)