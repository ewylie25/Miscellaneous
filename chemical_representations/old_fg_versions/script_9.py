from rdkit import Chem
import cPickle as pickle
from rdkit.Chem import Descriptors
import itertools
import sys
import json
from os import listdir
import numpy as np


def find_level(pair):
    match_1, match_2 = pair

    # Get all combos of atom indicies
    a_index_combos = itertools.product(match_1, match_2)
    levels = []
    for a_index_combo in a_index_combos:
        a_index_1, a_index_2 = a_index_combo
        levels.append(DM[a_index_1, a_index_2])
    return min(levels)


def missing_substructures(m, matches):
    if Descriptors.NumRadicalElectrons(m) > 0:
        return False
    else:
        ats = [i.GetIdx() for i in m.GetAtoms()]
        for i in matches:
            for j in i:
                ats = set(ats) - set(j)
        for i in ats:
            a = m.GetAtomWithIdx(i)
            if a.GetAtomicNum() != 6:
                if a.GetAtomicNum() != 1:
                    return True
        return False


if __name__ == '__main__':
    # To count total processed chemicals
    chem_count = 0

    # To count total features mapped
    feat_count = 0

    # To help catch missing features
    check = []

    # Numbered list of features - arbitrary, but needs to be consistent
    with open('new_legend.json') as f:
        fgs = json.load(f)

    # Make SMARTS for features
    smarts = {Chem.MolFromSmarts(i.encode('ascii')): j for j, i in fgs.items()}

    # Directory of files in json format containing chemical data
    # formatted as list of dictionaries. Each dictionary represents one chemical entry
    # with a 'smiles' key containing SMILES.
    dir_name = sys.argv[1]
    files = listdir(dir_name)

    # 3-dimensional matrix to collect data - the first dimension is the number of bonds between
    #     features
    #       0 - occurance of the group on the diagonal and overlap of groups off-diagonal
    #       1 - one bond separated
    #       2 - two bonds (and one atom) separated
    #       ... though 9 bonds
    #
    # the next two indicies correspond to features as layed out in fgs/smarts
    #       matrix is symmetric for this reason, if algorithm adds to [1,1,2] it should
    #        also add to [1,2,1]
    base = np.zeros((10, 361, 361), np.float128)

    for i in files:
        with open(dir_name + i) as f:
            data = json.load(f)

        for chem in data:
            smi = chem['smiles'].encode('ascii')
            mol = Chem.MolFromSmiles(smi)

            # Move on if it's no good
            if not mol:
                continue

            # Only count chemicals that contribute to the data collected
            chem_count += 1

            # Only care about features that are in the chemical
            potential = [j for j in smarts.keys() if mol.HasSubstructMatch(j)]

            # Store matches by feature in dictionary
            matches = {j: mol.GetSubstructMatches(j) for j in potential}

            # Check if molecule has interesting features not in feature list
            # Hap-hazard algorithm but gives a rough idea. Could be greatly
            # improved.
            if missing_substructures(mol, matches.values()):
                check.append(smi)

            # Calculate the distance matrix for the chemical graph
            DM = Chem.GetDistanceMatrix(mol)

            # Collect feature frequency and intra-group data
            for group in potential:
                g_matches = matches[group]
                number = len(g_matches)

                # Count the total number of features
                feat_count += number

                index = smarts[group]

                # Fill in frequency of group in database
                base[0, index, index] += number

                # Deal with frequency of multiple instances of the same group
                #   relative to eachother
                if number > 1:
                    combos = itertools.combinations(g_matches, 2)
                    levels = [find_level(c)
                              for c in combos if find_level(c) < 9]
                    for level in levels:
                        base[level, index, index] += 1

            # Get all permutations of features
            perms = itertools.permutations(potential, 2)
            for p in perms:
                first, second = p
                first_matches = matches[first]
                first_index = smarts[first]
                second_matches = matches[second]
                second_index = smarts[second]

                # Get all combos of groups with each other
                combos = itertools.product(first_matches, second_matches)
                levels = [find_level(c) for c in combos if find_level(c) < 9]
                for level in levels:
                    base[level, first_index, second_index] += 1

    print "Features: ", feat_count
    print "Chemicals: ", chem_count
    with open('counts.pickle', 'w') as f:
        pickle.dump(base, f)
    with open('check.json', 'w') as f:
        json.dump(check, f)
