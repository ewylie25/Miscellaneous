#!/usr/bin/env python

import numpy as np
from copy import deepcopy
from collections import Counter
from scipy import linalg
from rdkit import Chem
from logger import logger


def balance(smiles, explicit_h=False):
    """
    Function to balance a reaction by modifying the stoichometric coefficients. Input
    is a reactions SMILES string and explicit_h flag. If explicit_h, function counts
    hydrogen atoms. Output is new reactions SMILES string or None.

    Algorithm only considers each side independently - not ideal, but a first pass.
    """
    logger.info("Entering balance function.")
    logger.debug("Starting Smiles:{0}".format(smiles))

    # Split up SMILES and get Mol objects
    reactant_smiles, product_smiles = [i.split('.') for i in smiles.encode('ascii').split('>>')]
    reactants = [Chem.MolFromSmiles(r) for r in reactant_smiles]
    products = [Chem.MolFromSmiles(p) for p in product_smiles]

    if None in reactants or None in products:
        return None

    # Get hydrogens if necessary
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

    logger.debug("Difference: {0}".format(str(difference)))

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

    logger.debug("Atoms needed on reactant side: {0}".format(str(atoms_needed_reactant)))
    logger.debug("Atoms needed on product side: {0}".format(str(atoms_needed_product)))

    #Find candidate molecules for balancing reactant and product side.
    if atoms_needed_reactant:
        r_candidates = [i for i in reactants_atom_totals if candidate(i, atoms_needed_reactant)]
    else:
        r_candidates = None
    logger.debug("Candidate reactants: {0}".format(str(r_candidates)))

    if atoms_needed_product:
        p_candidates = [i for i in products_atom_totals if candidate(i, atoms_needed_product)]
    else:
        p_candidates = None
    logger.debug("Candidate products: {0}".format(str(p_candidates)))

    # Keep reference of what we had originally.
    new_reactant_smiles = deepcopy(reactant_smiles)
    new_product_smiles = deepcopy(product_smiles)

    # Try to balance and update lists as we go.
    if not p_candidates and not r_candidates:
        logger.debug("No Candidates.")
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

    # If there's no change - continue on...
    if new_reactant_smiles == reactant_smiles and new_product_smiles == product_smiles:
        logger.debug("No change.")
        return None
    new_smiles = '>>'.join(['.'.join(new_reactant_smiles), '.'.join(new_product_smiles)])

    logger.info("Leaving balance function.")

    # Double check that we didn't mess up.
    if not check_balance(new_smiles, explicit_h):
        logger.debug("New but not balanced.")
        return None
    logger.debug("New and balanced.")
    return new_smiles


def check_balance(smiles, explicit_h):
    """
    Given reaction smiles, checks that the reaction is balanced by comparing the
    number of atoms of each element on the reactant side to the product side. Checks
    hydrogen atoms if explicit_H is set to True, otherwise it only checks heavy atoms.
    """
    logger.info("Entering check_balance function.")

    # Split up SMILES and get Mol objects
    reactant_smiles, product_smiles = [i.split('.') for i in smiles.encode('ascii').split('>>')]
    reactants = [Chem.MolFromSmiles(r) for r in reactant_smiles]
    products = [Chem.MolFromSmiles(p) for p in product_smiles]

    if None in reactants or None in products:
        return False

    # Get list of total atoms as atomic symbols
    reactant_symbols = get_atom_info(reactants)
    product_symbols = get_atom_info(products)

    # Deal with hydrogens if necessary
    if explicit_h:
        reactant_hs = sum(map(get_Hs, reactants))
        product_hs = sum(map(get_Hs, products))
        reactant_atom_totals = count(reactant_symbols, reactant_hs)
        product_atom_totals = count(product_symbols, product_hs)
    else:
        reactant_atom_totals = count(reactant_symbols)
        product_atom_totals = count(product_symbols)

    logger.debug("Reactant atom totals: {0}".format(str(reactant_atom_totals)))
    logger.debug("Product atom totals: {0}".format(str(product_atom_totals)))
    logger.info("Leaving check_balance function.")

    return reactant_atom_totals == product_atom_totals


def get_Hs(mol):
    """
    Helper function.
    Returns the number of hydrogen atoms in a molecule.
    """
    logger.info("Entering get_Hs")
    total_Hs = 0
    for atom in mol.GetAtoms():
        total_Hs += atom.GetNumImplicitHs() + atom.GetNumExplicitHs()
    logger.info("Leaving get_Hs, totaled {0} hydrogens.".format(str(total_Hs)))
    return total_Hs


def get_atom_info(m):
    """
    Helper function.
    Generates a list of atomic symbols.
    It can take either a mol object or a list of mol objects.
    """
    logger.info("Entering get_atom_info")
    total = []
    if isinstance(m, list):
        for mol in m:
            atomic_symbol = [atom.GetSymbol() for atom in mol.GetAtoms()]
            total += atomic_symbol
    elif isinstance(m, Chem.rdchem.Mol):
        total = [atom.GetSymbol() for atom in m.GetAtoms()]
    else:
        logger.error("Incorrect type passed to get_atom_info:{0}".format(type(m)))
    logger.info("Leaving get_atom_info, listed: {0}". format(str(total)))
    return total


def count(list_symbols, hs=None):
    """
    Helper function.
    Returns counter. Adds 'H' if
    specified in function call.
    """
    logger.info("Entering count")
    c = Counter(list_symbols)
    if hs:
        c['H'] = hs
    logger.info("Leaving count, counted: {0}".format(str(c)))
    return c


def candidate(sub_set, super_set):
    """
    Helper function.
    Returns boolean indicating whether the 1st counter
    is a subset of the 2nd counter.
    """
    if set(sub_set).issubset(set(super_set)):
        if all([v <= super_set[k] for k, v in sub_set.items()]):
            return True
    return False


def find_decomposition(potential, needs):
    """
    Uses matrices to solve a system of linear equations.
    Takes list of counters representing molecules and counter representing the total number of missing atoms to be
    satisfied. Returns list of tuples indicating how to modify side of reaction.
    """
    logger.info("Entering find_decomposition.")
    key = []
    val = []
    for k, v in needs.items():
        key.append(k)
        val.append([v])
    b = np.array(val)
    val = []
    for k in key:
        v = []
        for p in potential:
            v.append(p[k])
        val.append(v)
    a = np.array(val)
    solution = linalg.lstsq(a, b)[0]
    results = []
    for i, j in enumerate(potential):
        results.append((j, solution[i]))
    if all([sum([j*i[k] for i, j in results]) <= v for k, v in needs.items()]):
        logger.info("leaving find_decomposition - Successful")
        logger.debug("Decomp: {0}".format(str(results)))
        return results
    else:
        logger.info("leaving find_decompositions - Failure")
        return None


def update(res, d_list, smi_list):
    """
    Helper function.
    Takes results from the linear algebra problem, the list of atom counts and list of smiles.
    Returns the updated list of smiles.
    """
    logger.info("Entering update.")
    if not res:
        logger.info("Leaving update with no change.")
        return smi_list
    for i in res:
        d, num = i
        s = smi_list[d_list.index(d)]
        init = smi_list.count(s)
        while smi_list.count(s) < num+init:
            smi_list.append(s)
    logger.info("Leaving update with change.")
    logger.debug("Updated to: {0}".format(str(smi_list)))
    return smi_list
