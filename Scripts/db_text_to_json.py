 #!/usr/bin/env python
import sys
import json


def main(text_file):
    # list of lines from file
    with open(text_file) as f:
        raw = f.readlines()
    # remove pesky line returns and split on tabs
    data = [i.strip('\n').split('\t') for i in raw]

    # a little defensive coding
    data = [i for i in data if i != ['']]

    entries = []

    for i in data:
        #print i
        d = {}
        d['name'] = i[0].encode('ascii')
        d['reaction_smarts'] = i[1].encode('ascii')
        d['products'] = i[2].encode('ascii').split('.')
        d['product_smiles'] = i[3].encode('ascii').split('.')
        d['synthon_count'] = int(i[4])
        d['explicit_H'] = (i[5] == 'True')
        d['reference'] = i[8]
        d['condition'] = i[9]
        d['Corey3'] = (i[10] == 'True')
        d['protected_smarts'] = i[11].encode('ascii').split('.')
        d['green_conditions'] = i[12].encode('ascii').split('.')
        d['incompat_smarts'] = i[13].encode('ascii').split('.')
        d['condition_utility'] = int(i[14])
        #d['substitution_details'] = {'atom_index':int(i[15]), 'energy_type':i[16].encode('ascii'), 'reactivity':i[17].encode('ascii')}
        d['r_check'] = True
        d['atom_economy'] = -1
        #d['special_Heterocycles']=True
        entries.append(d)

    for e in entries:
        if e['reference'] == '':
            del e['reference']
        if e['product_smiles'] == ['']:
            del e['product_smiles']

    print "Number of entries: " + str(len(entries))
    print entries[0]
    is_okay = raw_input("Does this info look okay(enter 'y' or 'n'): ")
    if is_okay == 'y':
        #start = int(raw_input("Enter starting _id number: "))
        #for j, i in enumerate(entries):
        #    i['_id'] = start + j
        #    i['old_ID'] = start + j
        with open( text_file.split('.')[0] + '.json', 'w') as f:
            json.dump(entries, f, sort_keys=True, indent=4)
    elif is_okay == 'n':
        print "Sorry - edit files and try again"
        sys.exit()
    else:
        print "Stop using this code!"
        sys.exit()


if __name__ == '__main__':
    source_file = sys.argv[1]
    main(source_file)
