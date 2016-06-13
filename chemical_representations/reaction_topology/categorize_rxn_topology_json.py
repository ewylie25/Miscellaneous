#!/usr/bin/env python
"""
Script to categorize (unbalanced) reactions from pickled _id list by reaction 'topology'.
Takes a dumped json with modified reactions.

Liz Wylie
January 28th 2014
"""

import pymongo
import sys
from collections import defaultdict
import json


def main(unbal_json):
    data = defaultdict(list)
    for i in unbal_json:
        #entry = db.reaction.find_one({'_id':i})
        smis = i['new']
        reactant_str, product_str = smis.split('>>')

        rs = len(reactant_str.split('.'))
        ps = len(product_str.split('.'))
        t = str(rs)+'->'+str(ps)
        data[t].append(i)
    for k,v in data.items():
        print k, str(len(v))
        with open(k+'.json', 'w') as f:
            json.dump(v, f)


if __name__ == '__main__':
    file_name = sys.argv[1]
    db = pymongo.MongoClient().new_data
    with open(file_name) as f:
        dat = json.load(f)
    main(dat)
