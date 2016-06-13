#!/usr/bin/env python
"""
Script to categorize (unbalanced) reactions from pickled _id list by reaction 'topology'.
Takes a pickled list of ids.

Liz Wylie
January 25th 2014
"""

import pymongo
import sys
from collections import defaultdict
import cPickle as pickle


def main(unbal_list):
    data = defaultdict(list)
    for i in unbal_list:
        entry = db.reaction.find_one({'_id':i})
        smis = entry['smiles'].encode('ascii')
        reactant_str, product_str = smis.split('>>')

        rs = len(reactant_str.split('.'))
        ps = len(product_str.split('.'))
        t = str(rs)+'->'+str(ps)
        data[t].append(i)
    for k,v in data.items():
        print k, str(len(v))
        with open(k+'.txt', 'w') as f:
            pickle.dump(v, f)


if __name__ == '__main__':
    file_name = sys.argv[1]
    db = pymongo.MongoClient().new_data
    with open(file_name) as f:
        dat = pickle.load(f)
    main(dat)
