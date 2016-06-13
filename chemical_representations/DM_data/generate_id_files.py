#!/usr/bin/env python

"""
Short script to generate files of _ids for distributing work load.

Liz Wylie
August 2014

"""

from pymongo import MongoClient
from pymongo.errors import ConnectionFailure
import os, errno

if __name__ == '__main__':

    # Create directory is it does not already exist, from http://stackoverflow.com/questions/600268/mkdir-p-functionality-in-python
    try:
        os.makedirs('./source_files')
    except OSError as exc: 
        if exc.errno == errno.EEXIST and os.path.isdir('./source_files'):
            pass
        else: raise

    # Connects to MongoDB on localhost unless altered or connection failure 
    try:
        db = MongoClient().new_data
    except ConnectionFailure:
        ip = raw_input("Please enter ip address of computer running MongoDB: ")
        db = MongoClient(host=ip).new_data

    # Creates cursor to iterate over the `chemical` collection in the mongoDB.
    # We only want entries to be analyzed once (snapshot=True) and 
    # we only care about entries with SMILES ({'smiles':{'$exists':True}}
    cursor = db.chemical.find({'smiles': {'$exists': True}}, snapshot=True)

    # Array to store _ids and a counter to keep track of file count 
    data = []
    counter = 0

    # count the number of entries we've processed as we go
    for count, entry in enumerate(cursor):
        # we need SMILES
        if not entry["smiles"]:
            continue
        data.append(entry['_id'])
        # dump the data every 500,000 and deal with the fact that enumerate() starts indexing at 0
        if (count + 1) % 500000 == 0:
            with open('_'.join(["source_files/chemical_ids", str(counter)])+".txt", 'w') as f:
                for value in data:
                    f.write(str(value)+'\n')
            data = []
            counter += 1
    # dump that last part of the database
    with open('_'.join(["source_files/chemical_ids", str(counter)])+".txt", 'w') as f:
        for value in data:
            f.write(str(value)+'\n')
