import elasticsearch
from scipy import sparse
import json

if __name__ == '__main__':
    es = elasticsearch.Elasticsearch()
    with open("DM_data/legend.json") as f:
        legend = json.load(f)
    while es.count(index=['dm_data'], doc_type=["chemical"], body={"query": {"match": {"status": "newdone"}}})['count'] > 0:
        data = es.search(index='dm_data', doc_type="chemical", body={
            "query": {"match": {"status": "newdone"}}}, size=500)
        data = data['hits']['hits']
        while data:
            entry = data.pop()
            base = [sparse.lil_matrix((len(legend), len(legend)))
                    for i in range(10)]
            info = entry['_source']
            info['status'] = 'final'
            for j in info['raw rings']:
                ring = j
                first_index = legend[j]
                for k in info['raw rings'][j]:
                    for l in k:
                        second = l
                        second_index = legend[l]
                        levels = [i for i in k[l] if i<=9]
                        for level in levels:
                            base[int(level)][first_index,second_index]+=1
            d = [b.nonzero() for b in base]
            new = {}
            for i, b in enumerate(d):
                points = zip(list(b[0]), list(b[1]))
                new_data = {}
                for j in points:
                    key = ".".join([str(k) for k in j])
                    new_data[key] = base[i][j]
                new[str(i)] = new_data
            total = info['processed']
            new_total = {}
            for k, v in total.items():
                t = dict(v.items() + new[k].items())
                new_total[k] = t
            new_info = {
                'body': new_total,
                'version': 'Ha_hybd_atom_aromat_1_Leg_hash_rings_1'
            }
            es.index(index='dm_final', doc_type='chemical', id=entry['_id'], body=new_info)
            del info['processed']
            es.index(index=entry['_index'], doc_type=entry['_type'], id=entry['_id'], body=info)