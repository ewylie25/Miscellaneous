import elasticsearch
import json
from elasticsearch import NotFoundError
from scipy import sparse


if __name__ == '__main__':
    es = elasticsearch.Elasticsearch()
    with open("DM_data/legend.json") as f:
        legend = json.load(f)
    while es.count(index=['dm_data'], doc_type=["chemical"], body={"query": {"match": {"status": "undone"}}})['count'] > 0:
        data = es.search(index='dm_data', doc_type="chemical", body={
            "query": {"match": {"status": "undone"}}}, size=500)
        data = data['hits']['hits']
        while data:
            entry = data.pop()
            base = [sparse.lil_matrix((len(legend), len(legend)))
                    for i in range(10)]
            info = entry['_source']
            info['status'] = 'done'
            if 'status' in info['body'] and info['body']['status'] == 'BAD SMILES':
                try:
                    es.delete(index=entry['_index'],
                              doc_type=entry['_type'], id=entry['_id'])
                    continue
                except NotFoundError:
                    continue
            for j in info['body']:
                first = j
                if first == 'status' or first == 'body':
                    continue
                first_index = legend[j]
                for k in info['body'][j]:
                    for l in k:
                        second = l
                        second_index = legend[l]
                        levels = [i for i in k[l] if i <= 9]
                        for level in levels:
                            base[int(level)][first_index, second_index] += 1
            d = [b.nonzero() for b in base]
            total = {}
            for i, b in enumerate(d):
                points = zip(list(b[0]), list(b[1]))
                new_data = {}
                for j in points:
                    key = ".".join([str(k) for k in j])
                    new_data[key] = base[i][j]
                total[i] = new_data
            info['processed'] = total
            es.index(index=entry['_index'], doc_type=entry['_type'], id=entry['_id'], body=info)