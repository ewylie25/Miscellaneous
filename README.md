Miscellaneous
=============

Miscellaneous scripts that are/have been/could be useful ideas/tools for solving challenges that come up in my research.

### Balancing Reactions:
Scripts to generate data files for yield prediction project.
#####Main code:
`python extract_balance_rxns.py [-h] [-i] [-o OUTPUT_FILE] (-l INTEGER_IDS [INTEGER_IDS ...] | -f FILE [FILE ...] | -c COUNT)`
#####Runs unit tests:
`python test.py`
#####Tool:
`python eval_error_script.py input`

Generates json file and counts of the results of the balancing code without taking into consideration entry data quality. Requires integer input denoted the number of entries to tes. Tool to explore how the code could be improved.

---
### Scripts:
Miscellaneous, single file scripts from various projects.
####dump_chemical_database.py
Used to generate .json files of chunks of the database. Useful for prep for running jobs on cluster or various other computers that may or may not have direct access to database.

