"""Read-only public snapshots and reproducible deletion-code preprocessing.

This client has no registration, POST, submission or credential functionality.
Downloaded verifiers are data: review their source before adding to verify.py.
"""
import argparse
import gzip
import hashlib
import json
from datetime import datetime, timezone
from urllib.request import urlopen
import numpy as np
from verify import ROOT, read_json

BASE='https://einsteinarena.com'

def snapshot(slug, limit=1):
    endpoints=[('problems_'+slug,f'/api/problems/{slug}')]
    problem=json.load(urlopen(BASE+endpoints[0][1],timeout=60))
    endpoints += [(slug+'-best',f'/api/solutions/best?problem_id={problem["id"]}&limit={limit}'),
                  (slug+'-leaderboard',f'/api/leaderboard?problem_id={problem["id"]}&limit=5')]
    manifest=[]
    for name,endpoint in endpoints:
        data=problem if name.startswith('problems_') else json.load(urlopen(BASE+endpoint,timeout=60))
        raw=json.dumps(data,ensure_ascii=False).encode()
        path=ROOT/'data'/(name+'.json')
        stored=path
        if len(raw)>200000:
            stored=path.with_suffix('.json.gz');stored.write_bytes(gzip.compress(raw,mtime=0))
            path.unlink(missing_ok=True)
        else:path.write_bytes(raw);path.with_suffix('.json.gz').unlink(missing_ok=True)
        manifest.append(dict(url=BASE+endpoint,file=stored.name,uncompressed_sha256=hashlib.sha256(raw).hexdigest(),fetched_at=datetime.now(timezone.utc).isoformat()))
    print(json.dumps(manifest,indent=2))

def deletion_precompute():
    source=read_json(ROOT/'data'/'two-deletion-code-16-best.json')[0]
    seed=np.array([int(w,2) for w in source['data']['words']])
    words=np.arange(65536,dtype=np.uint32);rows=[]
    for i in range(15):
        first=(words & ((1<<i)-1)) | ((words>>(i+1))<<i)
        for j in range(i+1,16):
            k=j-1
            rows.append((first & ((1<<k)-1)) | ((first>>(k+1))<<k))
    descendants=np.sort(np.array(rows,dtype=np.uint16).T,axis=1)
    # Separate exact string implementation checks bit indexing on fixed samples.
    for word in [0,1,3,21845,32768,65535]:
        s=format(word,'016b')
        strings={int(s[:i]+s[i+1:j]+s[j+1:],2) for i in range(15) for j in range(i+1,16)}
        assert set(descendants[word])==strings
    owner=np.full(16384,-1,dtype=np.int32)
    for index,w in enumerate(seed):
        assert np.all(owner[descendants[w]]==-1)
        owner[descendants[w]]=index
    blockers=np.sort(owner[descendants],axis=1)
    distinct=(blockers>=0)
    distinct[:,1:] &= blockers[:,1:]!=blockers[:,:-1]
    counts=distinct.sum(axis=1)
    np.savez_compressed(ROOT/'data'/'deletion-descendants.npz',descendants=descendants,conflict_counts=counts)
    print(json.dumps(dict(source_solution_id=source['id'],words=65536,positions_per_word=120,pool_at_most_four_blockers=int((counts<=4).sum()))))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--snapshot');p.add_argument('--limit',type=int,default=1)
    p.add_argument('--deletion-precompute',action='store_true');a=p.parse_args()
    if a.snapshot:snapshot(a.snapshot,a.limit)
    elif a.deletion_precompute:deletion_precompute()
    else:p.print_help()
