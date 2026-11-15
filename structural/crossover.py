"""Large cross-layer recombination, with exact bipartite union bounds.

The union of two feasible packings induces a bipartite conflict graph after
removing their common (isolated) words. A matching/vertex-cover certificate gives
the best possible recombination in that union. Three-or-more parents use CP-SAT;
optional graph clique cuts include conflicts not sharing a single received word.
"""
import argparse
import hashlib
import itertools
import json
import time
from pathlib import Path
import numpy as np
from fixed_verifier import verify

ROOT=Path(__file__).resolve().parent

def load_parents():
    rows=json.loads((ROOT/'public-parents.json').read_text())['parents'];parents=[]
    for row in rows:
        words=row['data']['words'];verify(words)
        base={int(w,2) for w in words}
        for transform in ['identity','complement','reverse','reverse-complement']:
            S={int(format(w,'016b')[::-1],2) if 'reverse' in transform else w for w in base}
            if 'complement' in transform:S={w^65535 for w in S}
            parents.append(dict(source_id=row['id'],source_agent=row['agentName'],transform=transform,words=S))
    return parents

def pair_optimum(A,B,d):
    common=A&B;left=sorted(A-B);right=sorted(B-A)
    owner=np.full(16384,-1,dtype=int)
    for i,w in enumerate(right):owner[d[w]]=i
    edges=[sorted(set(owner[d[w]])-{-1}) for w in left]
    mateR=[-1]*len(right)
    def dfs(u,seen):
        for v in edges[u]:
            if v in seen:continue
            seen.add(v)
            if mateR[v]<0 or dfs(mateR[v],seen):mateR[v]=u;return True
        return False
    matched=sum(dfs(u,set()) for u in range(len(left)))
    mateL=[-1]*len(left)
    for v,u in enumerate(mateR):
        if u>=0:mateL[u]=v
    zL={u for u in range(len(left)) if mateL[u]<0};zR=set();queue=list(zL)
    while queue:
        u=queue.pop()
        for v in edges[u]:
            if v not in zR:
                zR.add(v);u2=mateR[v]
                if u2>=0 and u2 not in zL:zL.add(u2);queue.append(u2)
    coverL=set(range(len(left)))-zL
    assert len(coverL)+len(zR)==matched
    assert all(u in coverL or v in zR for u,row in enumerate(edges) for v in row)
    result=common|{left[u] for u in zL}|{right[v] for v in range(len(right)) if v not in zR}
    assert len(result)==len(A|B)-matched
    return result,dict(union_size=len(A|B),matching=matched,optimum=len(result),common=len(common))

def pair_scan(args):
    parents=load_parents();d=np.load(ROOT.parent/'data'/'deletion-descendants.npz')['descendants']
    out=ROOT/'runs'/args.name;out.mkdir(parents=True,exist_ok=True);start=time.perf_counter();records=[]
    best=max(len(p['words']) for p in parents)
    for i,j in itertools.combinations(range(len(parents)),2):
        candidate,certificate=pair_optimum(parents[i]['words'],parents[j]['words'],d)
        verify([format(w,'016b') for w in candidate])
        row=dict(parents=[i,j],**certificate);records.append(row)
        if len(candidate)>best:
            best=len(candidate);(out/'best-valid.json').write_text(json.dumps({'words':[format(w,'016b') for w in sorted(candidate)]},indent=2)+'\n')
            print(json.dumps(dict(event='verified_improvement',size=best,parents=[i,j])),flush=True)
    metadata=[{k:v for k,v in p.items() if k!='words'} for p in parents]
    summary=dict(method='exact pair-union recombination',parents=metadata,pairs=records,best=best,elapsed_seconds=time.perf_counter()-start,
                 interpretation='Each optimum is restricted to the named pair union, not a global code bound.',publicly_submitted=False)
    (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps({k:v for k,v in summary.items() if k not in ['parents','pairs']}),flush=True)

def cp_union(args):
    from ortools.sat.python import cp_model
    parents=load_parents();chosen=list(map(int,args.parents.split(',')));P=[parents[i]['words'] for i in chosen]
    arrays=np.load(ROOT.parent/'data'/'deletion-descendants.npz');d=arrays['descendants']
    pool_set=set.union(*P)
    if args.max_blockers>=0:
        assert chosen[0]==0,'Stored blocker counts refer to the 226-word source parent'
        pool_set.update(map(int,np.flatnonzero(arrays['conflict_counts']<=args.max_blockers)))
    pool=sorted(pool_set)
    model=cp_model.CpModel();variables=[model.new_bool_var(f'w{w}') for w in pool]
    bins=[[] for _ in range(16384)]
    for i,w in enumerate(pool):
        for child in set(d[w]):bins[child].append(i)
    cliques={tuple(row) for row in bins if len(row)>1}
    adjacency=[set() for _ in pool]
    for row in cliques:
        for i in row:adjacency[i].update(set(row)-{i})
    original_count=len(cliques)
    if args.clique_cuts:
        for v in range(len(pool)):
            remaining=adjacency[v].copy();C={v}
            while remaining:
                w=max(remaining,key=lambda w:len(remaining&adjacency[w]));C.add(w);remaining &=adjacency[w]
            if len(C)>2:cliques.add(tuple(sorted(C)))
    for row in cliques:model.add_at_most_one([variables[i] for i in row])
    model.add(sum(variables)>=args.target)
    if args.maximize:model.maximize(sum(variables))
    for w,v in zip(pool,variables):model.add_hint(v,int(w in P[0]))
    out=ROOT/'runs'/args.name;out.mkdir(parents=True,exist_ok=True);handle=(out/'decisions.jsonl').open('w');start=time.monotonic()
    def log(row):
        row['seconds']=time.monotonic()-start;handle.write(json.dumps(row)+'\n');handle.flush();print(json.dumps(row),flush=True)
    log(dict(event='start',parents=chosen,pool_size=len(pool),max_blockers=args.max_blockers,slot_cliques=original_count,graph_cliques=len(cliques)-original_count,target=args.target,seconds_budget=args.seconds,clique_cuts=args.clique_cuts))
    solver=cp_model.CpSolver();solver.parameters.max_time_in_seconds=args.seconds;solver.parameters.num_search_workers=1;solver.parameters.random_seed=args.seed
    class Callback(cp_model.CpSolverSolutionCallback):
        def __init__(self):super().__init__();self.best=0
        def on_solution_callback(self):
            candidate=[format(w,'016b') for w,v in zip(pool,variables) if self.value(v)];result=verify(candidate);self.best=max(self.best,len(candidate))
            (out/'best-valid.json').write_text(json.dumps({'words':candidate},indent=2)+'\n');log(dict(event='independently_verified_candidate',**result))
            if not args.maximize:self.stop_search()
    cb=Callback();status=solver.solve(model,cb)
    summary=dict(status=solver.status_name(status),best=cb.best,parents=chosen,pool_size=len(pool),max_blockers=args.max_blockers,target=args.target,clique_cuts=args.clique_cuts,
                 seconds=time.monotonic()-start,solver_cpu_budget=args.seconds,publicly_submitted=False)
    (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');log(dict(event='finish',**summary))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--kind',choices=['pairs','cp'],default='pairs');p.add_argument('--parents',default='0,1,2,3');p.add_argument('--target',type=int,default=227)
    p.add_argument('--seconds',type=int,default=300);p.add_argument('--name',default='pair-scan');p.add_argument('--seed',type=int,default=42)
    p.add_argument('--max-blockers',type=int,default=-1)
    p.add_argument('--clique-cuts',action='store_true');p.add_argument('--maximize',action='store_true');a=p.parse_args()
    if a.kind=='pairs':pair_scan(a)
    else:cp_union(a)
