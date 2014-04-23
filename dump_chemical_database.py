#!/usr/bin/env python
"""
Script extracting data from chemical database

Liz Wylie
Feburary 2014
"""

import json
import pymongo


def dump_data(results, index):
    with open('tmp/data_'+str(index)+'.json', 'w') as f:
        json.dump(results,f, indent=4)


def main():
    cursor = db.chemical.find(snapshot=True)
    results = []
    for i, entry in enumerate(cursor):
        data = {}
        if not entry.get('smiles'):
            continue
        data['_id'] = entry['_id']
        data['smiles'] = entry['smiles'].encode('ascii')
        if not entry.get('functional_groups'):
            continue
        data['functional_groups'] = entry['functional_groups']
        results.append(data)
        if len(results) >100000:
            dump_data(results, i)
            results = []
    cursor.close()
    dump_data(results, i)
    return i

if __name__ == '__main__':
    db = pymongo.MongoClient().new_data
    entries = main()
    print "Processed Chemicals: ", entries
