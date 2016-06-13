#!/usr/bin/env python

"""
Full script to generate data from initial chemical entry in database.

Liz Wylie
August 2014

"""

from rdkit import Chem
import argparse
import elasticsearch
import json
from collections import defaultdict
from pymongo import MongoClient
from pymongo.errors import ConnectionFailure
from logger import logger
from scipy import sparse


def label(mol):
    for atom in mol.GetAtoms():
        hybridization = int(atom.GetHybridization())
        atom_num = atom.GetAtomicNum()
        aromaticity = int(atom.GetIsAromatic())
        hashed = int(str(hybridization) + str(atom_num) + str(aromaticity))
        atom.SetIsotope(hashed)


def build_m_of_m(data_set):
    base = [sparse.lil_matrix((len(legend), len(legend)))
            for number in range(10)]
    for j in data_set:
        first_index = legend[str(j)]
        for p in data_set[j]:
            for l in p:
                second_index = legend[str(l)]
                f_levels = [z for z in p[l] if z <= 9]
                for f_level in f_levels:
                    base[int(f_level)][first_index, second_index] += 1
    d = [b.nonzero() for b in base]
    temp_total = {}
    for h, b in enumerate(d):
        points = zip(list(b[0]), list(b[1]))
        new_data = {}
        for j in points:
            key_string = ".".join([str(z) for z in j])
            new_data[key_string] = base[h][j]
        temp_total[h] = new_data
    return temp_total


if __name__ == '__main__':
    # Set up the argument parsing
    logger.info("Setting up parser")
    parser = argparse.ArgumentParser(
        description="Script to convert chemical SMILES to topological matrix \
        representation used as an alternative to Fingerprints.", epilog="Contact Liz Wylie,\
        ewylie25@gmail.com for help.")
    parser.add_argument('-f', '--file', nargs='+',
                        type=str, help="specify file or files containing ids")

    # Connect to databases etc.
    try:
        logger.info("Connecting to MongoDB running on localhost.")
        db = MongoClient().liz_data
    except ConnectionFailure:
        logger.info("Connecting to localhost, failed.")
        ip = raw_input("Please enter ip address of computer running MongoDB: ")
        logger.info(
            "Attempting to connect to MongoDB running on {0}.".format(ip))
        db = MongoClient(host=ip).liz_data
        logger.info("Connected to MongoDB.")

    es = elasticsearch.Elasticsearch()

    # parsing arguments
    logger.info("Parsing arguments.")
    args = parser.parse_args()

    # Deal with input - if not files specified, use all of database
    if args.file:
        raw = []
        for file_name in args.file:
            logger.debug("Reading in _ids from file: {0}".format(file_name))
            with open(file_name, 'r') as f:
                raw += f.readlines()
        try:
            ids = [int(i.strip('\n')) for i in raw if i]
        except Exception as err:
            logger.error(
                "Input file(s) improperly formatted with error:{0}".format(err))
            raise Exception(
                "File(s) not formatted correctly. Should have one integer _id per line, no trailing whitespace.")
        cursor = db.chemical.find(
            {'$and': [{"_id": {'$in': ids}}, {"processed": False}]}, timeout=False)
    else:
        cursor = db.chemical.find(
            {'$and': [{'smiles': {'$exists': True}}, {'processed': False}]}, snapshot=True, timeout=False)

    with open("DM_data/legend.json") as f:
        legend = json.load(f)

    # data structure for keeping track of all of the existing and new hashes
    types = set(legend.keys())
    new_types = set()

    # ring info
    ring_smis = ["[*]1[*][*]1",
                 "[*]1[*][*][*]1",
                 "[*]1[*][*][*][*]1",
                 "[*]1[*][*][*][*][*]1",
                 "[*]1[*][*][*][*][*][*]1",
                 "[*]1[*][*][*][*][*][*][*]1",
                 "[*]1[*][*][*][*][*][*][*][*]1",
                 "[*]1[*][*][*][*][*][*][*][*][*]1",
                 "[*]1[*][*][*][*][*][*][*][*][*][*]1",
                 "[*]1[*][*][*][*][*][*][*][*][*][*][*]1",
                 "[*]1[*][*][*][*][*][*][*][*][*][*][*][*]1"]
    rings = [Chem.MolFromSmarts(i) for i in ring_smis]

    # iterate through db
    for chem in cursor:
        if not chem["smiles"]:
            db.chemical.update({'_id': chem['_id']}, {
                               '$set': {'status': 'fail', 'processed': True}})
            continue
        smi = chem['smiles']
        _id = chem['_id']
        m = Chem.MolFromSmiles(smi.encode('ascii'))
        if not m:
            db.chemical.update({'_id': chem['_id']}, {
                               '$set': {'status': 'fail', 'processed': True}})
            continue
        DM = Chem.GetDistanceMatrix(m)
        label(m)
        molecule_key = {}
        for a in m.GetAtoms():
            molecule_key[a.GetIdx()] = a.GetIsotope()
            if not str(a.GetIsotope()) in types:
                logger.warning("New atom hash type. Reconsider current legend")
                new_types.add(a.GetIsotope())
        if new_types:
            db.chemical.update({'_id': chem['_id']}, {
                               '$set': {'status': 'new_types', 'processed': True}})
            continue

        ring_matches = {ring_query_mol: m.GetSubstructMatches(
            ring_query_mol) for ring_query_mol in rings if m.HasSubstructMatch(ring_query_mol)}

        raw_ring_info = defaultdict(list)
        for key, value in ring_matches.items():
            temp = defaultdict(list)
            for match in value:
                ring_indices = set(match)
                non_ring_indices = set(molecule_key.keys()) - ring_indices
                for index_0 in ring_indices:
                    level = 0
                    temp[molecule_key[index_0]].append(level)
                for index_0 in non_ring_indices:
                    levels = []
                    for index_1 in ring_indices:
                        levels.append(DM[index_0, index_1])
                    level = min(levels)
                    temp[molecule_key[index_0]].append(level)
            raw_ring_info[ring_smis[rings.index(key)]].append(temp)

        data = defaultdict(list)
        for key in molecule_key.keys():
            temp = defaultdict(list)
            for index, value in enumerate(DM[key]):
                temp[molecule_key[index]].append(value)
            data[molecule_key[key]].append(temp)

        total_regular_atoms = build_m_of_m(data)
        total_rings = build_m_of_m(raw_ring_info)

        total = {}
        for k, v in total_regular_atoms.items():
            temp = dict(v.items() + total_rings[k].items())
            total[k] = temp
        info = {
            'smiles': smi,
            'body': total,
            'version': 'Ha_hybd_atom_aromat_1_Leg_hash_rings_1'
        }
        es.index(index='dm_final', doc_type='chemical', id=_id, body=info)
        db.chemical.update({'_id': chem['_id']}, {
                           '$set': {'status': 'success', 'processed': True}})
    cursor.close()
    logger.info("Completed {0}".format(file_name))
    logger.info("New Hashes: {0}".format(".".join([str(i) for i in new_types])))
