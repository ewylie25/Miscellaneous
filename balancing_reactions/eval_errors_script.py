#!/usr/bin/env python
"""
Script to evaluate what is missing in the balancing script by looking at properties of the 
reactions that were under-determined and cannot be balanced by balance function.
Produces a json with modified transforms.

Liz Wylie
March 2014
"""

import sys
from pymongo import MongoClient
from collections import defaultdict
from balancer import *
import json


def test_balance(smiles, explicit_h=False):
    """
    Function based off of balance function to return all attempts including incorrectly balanced reactions.
    """
    # Split up SMILES and get Mol objects
    reactant_smiles, product_smiles = [i.split('.') for i in smiles.encode('ascii').split('>>')]
    reactants = [Chem.MolFromSmiles(r) for r in reactant_smiles]
    products = [Chem.MolFromSmiles(p) for p in product_smiles]

    if None in reactants or None in products:
        return None

    # Get hydrogen atoms if necessary
    if explicit_h:
        reactant_hs = sum([get_Hs(i) for i in reactants])
        product_hs = sum([get_Hs(i) for i in products])

    # Get list of total atoms as atomic symbols
    reactant_symbols = get_atom_info(reactants)
    product_symbols = get_atom_info(products)

    # Find the in-balances
    if explicit_h:
        difference = count(reactant_symbols, reactant_hs)
        difference.subtract(count(product_symbols, product_hs))
    else:
        difference = count(reactant_symbols)
        difference.subtract(count(product_symbols))

    # Get list of atoms by molecule as atomic symbols
    reactants_atom_totals = [count(get_atom_info(i), get_Hs(i)) if explicit_h
                             else count(get_atom_info(i)) for i in reactants]
    products_atom_totals = [count(get_atom_info(i), get_Hs(i)) if explicit_h else
                            count(get_atom_info(i)) for i in products]

    # Keep positive values as atoms needed on product side
    atoms_needed_product = deepcopy(difference)
    atoms_needed_product += Counter()

    # Obtain negative values as atoms needed on reactant side.
    atoms_needed_reactant = Counter({k: abs(v) for k, v in difference.items() if v < 0})

    #Find candidate molecules for balancing reactant and product side.
    if atoms_needed_reactant:
        r_candidates = [i for i in reactants_atom_totals if candidate(i, atoms_needed_reactant)]
    else:
        r_candidates = None
    if atoms_needed_product:
        p_candidates = [i for i in products_atom_totals if candidate(i, atoms_needed_product)]
    else:
        p_candidates = None

    # Keep reference of what we had originally.
    new_reactant_smiles = deepcopy(reactant_smiles)
    new_product_smiles = deepcopy(product_smiles)

    # Try to balance and update lists as we go.
    if not p_candidates and not r_candidates:
        return None
    elif not p_candidates:
        prelim = find_decomposition(r_candidates, atoms_needed_reactant)
        new_reactant_smiles = update(prelim, reactants_atom_totals, new_reactant_smiles)
    elif not r_candidates:
        prelim = find_decomposition(p_candidates, atoms_needed_product)
        new_product_smiles = update(prelim, products_atom_totals, new_product_smiles)
    else:
        #TODO: This should be rewritten to consider both sides simultaneously.
        prelim = find_decomposition(r_candidates, atoms_needed_reactant)
        new_reactant_smiles = update(prelim, reactants_atom_totals, new_reactant_smiles)
        prelim = find_decomposition(p_candidates, atoms_needed_product)
        new_product_smiles = update(prelim, products_atom_totals, new_product_smiles)

    new_smiles = '>>'.join(['.'.join(new_reactant_smiles), '.'.join(new_product_smiles)])
    return new_smiles


def main(int_count):
    db = MongoClient(host="129.105.205.35").new_data
    results = []
    counts = defaultdict(int)
    for entry in db.reaction.find(snapshot=True).limit(int_count):
        new_reaction = test_balance(entry['smiles'])
        if not new_reaction:
            counts['no_change'] += 1
        else:
            if check_balance(new_reaction, False):
                counts['balanced'] += 1
            else:
                counts['unbalanced'] += 1
        results.append({'original': entry['smiles'],
                        'new': new_reaction})

    with open('balance_output.json', 'w') as f:
        json.dump(results, f)
    for k, v in counts.items():
        print k + '\t' + str(v)


if __name__ == '__main__':
    c = sys.argv[1]
    main(int(c))
