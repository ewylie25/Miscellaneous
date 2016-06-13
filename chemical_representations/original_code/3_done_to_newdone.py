import elasticsearch
from rdkit import Chem
from pymongo import MongoClient
from collections import defaultdict


def label(m):
    for a in m.GetAtoms():
        hybd = int(a.GetHybridization())
        atom = a.GetAtomicNum()
        aromat = int(a.GetIsAromatic())
        hashed = int(str(hybd) + str(atom) + str(aromat))
        a.SetIsotope(hashed)


if __name__ == '__main__':
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
    es = elasticsearch.Elasticsearch()
    db = MongoClient().liz_data
    while es.count(index=['dm_data'], doc_type=["chemical"], body={"query": {"match": {"status": "done"}}})['count'] > 0:
        data = es.search(index='dm_data', doc_type="chemical", body={
                         "query": {"match": {"status": "done"}}}, size=1000)
        data = data['hits']['hits']
        while data:
            entry = data.pop()
            info = entry['_source']
            _id = entry['_id']
            chem = db.chemical.find_one({"_id": int(_id)})
            smiles = chem['smiles'].encode('ascii')
            m = Chem.MolFromSmiles(smiles)
            if not m:
                continue
            DM = Chem.GetDistanceMatrix(m)
            label(m)
            key = {}
            for a in m.GetAtoms():
                key[a.GetIdx()] = a.GetIsotope()
            d = defaultdict(list)
            matches = {i: m.GetSubstructMatches(
                i) for i in rings if m.HasSubstructMatch(i)}
            for k, v in matches.items():
                t = defaultdict(list)
                for match in v:
                    ring_indicies = set(match)
                    non_ring_indicies = set(key.keys()) - ring_indicies
                    for i in ring_indicies:
                        level = 0
                        t[key[i]].append(level)
                    for i in non_ring_indicies:
                        levels = []
                        for j in ring_indicies:
                            levels.append(DM[i, j])
                        level = min(levels)
                        t[key[i]].append(level)
                d[ring_smis[rings.index(k)]].append(t)
            info['status'] = 'newdone'
            info['raw rings'] = d
            es.index(index=entry['_index'], doc_type=entry[
                     '_type'], id=entry['_id'], body=info)
