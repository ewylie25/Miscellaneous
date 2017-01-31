Miscellaneous
=============

Miscellaneous scripts and projects

### Balancing Reactions:
Scripts to generate data files for yield prediction project/paper.
#####Main code:
`python extract_balance_rxns.py [-h] [-i] [-o OUTPUT_FILE] (-l INTEGER_IDS [INTEGER_IDS ...] | -f FILE [FILE ...] | -c COUNT)`
#####Runs unit tests:
`python test.py`
#####Tool:
`python eval_error_script.py input`

Generates json file and counts of the results of the balancing code without taking into consideration entry data quality. Requires integer input denoted the number of entries to test. Tool to explore how the code could be improved.

---

### Chemical Representations:
Exploratory scripts from work around applying machine learning to chemical and reaction data for chematica stuff.

---
### Scripts:
Miscellaneous, single file scripts from various projects/tasks.
####dump_chemical_database.py
Used to generate .json files of chunks of the database. Useful for prep for running jobs on cluster or various other computers that may or may not have direct access to database.
####db_text_to_json.py
####generate_htmlpage_retrodb.py
Used to create html page of information from retro database. Beware - generates many-many png files. Very simplistic.
####generate_assets.sql
####trend_to_usable_csv.py

---
### SBBRG Project:
PHP/My SQL/HTML scripting task. See README.md in directory.
