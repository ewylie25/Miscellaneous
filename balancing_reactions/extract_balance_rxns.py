#!/usr/bin/env python
"""
Script to balance reactions only needing an additional stoiciometric equivalent of a reactant.
Must provide a count or list of _ids to use

Liz Wylie
February 2014
"""

from pymongo import MongoClient
from pymongo.errors import ConnectionFailure
import argparse
from logger import logger
import time
from balancer import check_balance, balance


def check_entry(db_entry):
    logger.info("Entering check_entry.")
    logger.debug("Yield field: {0}".format(db_entry["yield"]))
    if db_entry["yield"] != 0.75:
        logger.info("Leaving check_entry under first condition.")
        return True
    if not db_entry["rxid"]:
        logger.info("Leaving check_entry under second condition.")
        return False
    logger.debug("RXID field of entry: {0}".format(str(db_entry["rxid"])))
    for rxid in db_entry["rxid"]:
        key = "rx" + str(rxid)
        logger.debug("Key: {0}".format(key))
        if key in db_entry:
            logger.debug("Key entry: {0}".format(db_entry[key]))
            yields = [i["yield"] for i in db_entry[key] if "yield" in i]
            logger.debug(str(yields))
            logger.debug(str([j != 0.75 for k in yields for j in k]))
            if any([j != 0.75 for k in yields for j in k]):
                logger.info("Leaving check_entry under third condition.")
                return True
    logger.info("Leaving check_entry under final condition.")
    return False


def generate_data(db_entry):
    logger.info("Entering generate_data")
    dat = []
    yields = [db_entry["yield"]]
    logger.debug("Yield field of entry: {0}".format(str(yields)))
    if db_entry["rxid"]:
        logger.debug("RXID field of entry: {0}".format(str(db_entry["rxid"])))
        for rxid in db_entry["rxid"]:
            key = "rx" + str(rxid)
            if key in db_entry:
                yields += [i["yield"] for i in db_entry[key] if "yield" in i]
    logger.debug("All yields: {0}".format(str(yields)))
    good_yields = [i for i in yields if i !=0.75]
    logger.debug("Good yields: {0}".format(str(good_yields)))
    for y in good_yields:
        for x in y:
            tmp = [
                db_entry['_id'],
                db_entry['smiles'],
                db_entry['solvent'] if 'solvent' in db_entry else '000',
                db_entry['min_temp'] if 'min_temp' in db_entry else '0.0',
                db_entry['max_temp'] if 'max_temp' in db_entry else '0.0',
                x
            ]
        logger.debug("Temp: {0}".format(str(tmp)))
        tmp = map(str, tmp)
        dat.append(tmp)
    logger.info("Leaving generate_data with data: {0}".format(str(dat)))
    return dat


if __name__ == '__main__':
    #Initialize structure for storing results before writing to file
    data = []

    #Set up the argument parsing
    logger.info("Setting up parser")
    parser = argparse.ArgumentParser(description="Script to balance reactions and retrieve relevant information for \
                                        training yield prediction algorithm.", epilog="Contact Liz Wylie, ewylie25@gmail.com for help.")
    parser.add_argument('-i', '--include_hydrogens', action='store_true', help="reactions will be balanced with hydrogens")
    parser.add_argument('-o', '--output_file', type=str, help="filename to write output data, otherwise defaults to 'bal_rxns_output_%Y%m%d_%H%M.txt' ")
    #Deal with different input types
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('-l', '--integer_ids', nargs='+', type=int, help="specify reaction ids")
    group.add_argument('-f', '--file', nargs='+', type=str, help="specify file or files containing ids")
    group.add_argument('-c', '--count', type=int, help="specify the number of reactions desired in the output file")

    #Connect to database
    try:
        logger.info("Connecting to MongoDB running on localhost.")
        db = MongoClient().new_data
    except ConnectionFailure:
        logger.info("Connecting to localhost, failed.")
        ip = raw_input("Please enter ip address of computer running MongoDB: ")
        logger.info("Attempting to connect to MongoDB running on {0}.".format(ip))
        db = MongoClient(host=ip).new_data
        logger.info("Connected to MongoDB.")

    #Collect git sCL arguments
    logger.info("Parsing arguments.")
    args = parser.parse_args()

    #Assign variables accordingly and create cursor.
    logger.info("Processing arguments and setting up cursor.")
    if args.file:
        raw = []
        for file_name in args.file:
            logger.debug("Reading in _ids from file: {0}".format(file_name))
            with open(file_name, 'r') as f:
                raw += f.readlines()
        try:
            ids = [int(i.strip('\n')) for i in raw if i]
        except Exception as err:
            logger.error("Input file(s) improperly formatted with error:{0}".format(err))
            raise Exception("File(s) not formatted correctly. Should have one integer _id per line, no trailing whitespace.")
        cursor = db.reaction.find({"_id": {'$in': ids}})
    elif args.count:
        logger.debug("Collecting {0} reactions.".format(str(args.count)))
        cursor = db.reaction.find(snaptshot=True)
    else:
        logger.debug("Reading in _ids from CL: {0}".format(str(args.integer_ids)))
        ids = args.integer_ids
        cursor = db.reaction.find({"_id": {'$in': ids}})

    logger.debug("Explicit hydrogen atoms flag set to: {0}".format(str(args.include_hydrogens)))

    #Iterate through reactions, making sure cursor is closed even if code errors.
    try:
        for entry in cursor:
            logger.info("New Entry")
            if not "smiles" in entry:
                continue
            if not entry["smiles"]:
                continue
            logger.debug("Reaction {0} with smiles: {1}". format((str(entry["_id"])), entry["smiles"].encode('ascii')))
            if not check_entry(entry):
                logger.debug("Entry failed yield requirements.")
                continue
            #Check if Reaction's SMILES are balanced.
            if not check_balance(entry['smiles'], args.include_hydrogens):
                reaction = balance(entry['smiles'], args.include_hydrogens)
                if reaction:
                    entry['smiles'] = reaction
                else:
                    continue
            # Collect info from entry
            desired_info = generate_data(entry)
            data += desired_info
            if args.count:
                if len(data) >= args.count:
                    break
    finally:
        cursor.close()

    output = args.output_file if args.output_file else 'bal_rxns_output_{0}.txt'.format(time.strftime('%Y%m%d_%H%M'))
    with open(output, 'w') as f:
        for i in data:
            f.write(' '.join(i)+'\n')