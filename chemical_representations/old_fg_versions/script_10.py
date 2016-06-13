from rdkit import Chem
import itertools
import sys
import json
import pymongo
import bson
import cPickle as pickle
from os import listdir
import numpy as np
from scipy import sparse


def find_level(l):
    match_1, match_2 = l

    a_index_combos = itertools.product(match_1, match_2)
    levels = []
    for a_index_combo in a_index_combos:
        a_index_1, a_index_2 = a_index_combo
        levels.append(DM[a_index_1, a_index_2])
    return min(levels)


if __name__ == '__main__':
    db = pymongo.MongoClient().liz_data
    chem_count = 0
    with open('legend.json') as f:
        fgs = json.load(f)
    smarts = {Chem.MolFromSmarts(i.encode('ascii')):n for n,i in fgs.items()}
    dir_name = sys.argv[1]
    files = listdir(dir_name)
    for i in files[:1]:
        print i
        with open(dir_name+i) as f:
            data = json.load(f)
        data = data[:5]
        while data:
            chem = data.pop()
            base = [sparse.lil_matrix((361,361)) for i in range(10)]
            smi = chem['smiles'].encode('ascii')
            mol = Chem.MolFromSmiles(smi)
            if not mol:
                continue
            chem_count += 1
            potential = [j for j in smarts.keys() if mol.HasSubstructMatch(j)]
            matches = {j:mol.GetSubstructMatches(j) for j in potential}
            DM = Chem.GetDistanceMatrix(mol)
            
            for group in potential:
                g_matches = matches[group]
                number = len(g_matches)
                index = smarts[group]
                base[0][index, index] += number
                if number > 1:
                    combos = itertools.combinations(g_matches, 2)
                    levels = [find_level(c) for c in combos if find_level(c) < 9]
                    for level in levels:
                        base[level][index, index] += 1
            perms = itertools.permutations(potential,2)
            for p in perms:
                first, second = p
                d_1 = matches[first]
                d_2 = matches[second]
                n_1 = smarts[first]
                n_2 = smarts[second]
                combos = itertools.product(d_1, d_2)
                levels = [find_level(c) for c in combos if find_level(c) < 9]
                for level in levels:
                    base[level][n_1, n_2] += 1
            
            with open(str(chem['_id'])+'.dat', 'w') as f:
                pickle.dump(base,f, protocol = 2)
    print chem_count
