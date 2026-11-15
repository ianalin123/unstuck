"""Local, finite, auditable optimization. No paid API calls or arena writes.

The autocorrelation gradient is derived analytically from linear convolution.
The smoothing objective guides search; only the unsmoothed score accepts a result.
"""
import argparse
import json
import time
from pathlib import Path
import numpy as np
from scipy.fft import rfft, irfft, next_fast_len
from scipy.optimize import minimize
from scipy.special import logsumexp
from verify import ROOT, fast_score, read_json

class Objective:
    def __init__(self, f, kind, beta, keep, deadline, log):
        self.n=len(f);self.size=next_fast_len(2*self.n-1);self.kind=kind
        self.beta=beta;self.keep=keep;self.deadline=deadline;self.log=log
        self.calls=0;self.start=time.perf_counter()
    def __call__(self, x):
        if time.monotonic()>self.deadline: raise TimeoutError
        n=self.n
        if self.kind=='third': f=x + 1 - np.mean(x)
        else: f=x / np.mean(x)
        spectrum=rfft(f,self.size)
        c=irfft(spectrum*spectrum,self.size)[:2*n-1]/n
        target=np.abs(c) if self.kind=='third' else c
        maximum=float(target.max())
        weights=np.exp(self.beta*(target-maximum));weights/=weights.sum()
        soft=maximum+np.log(np.exp(self.beta*(target-maximum)).sum())/self.beta
        if self.kind=='third':weights*=np.sign(c)
        if self.kind=='second':
            q=(2*np.dot(c,c)+np.dot(c[:-1],c[1:]))/3
            # mean(f)=1 => sum(c)=n. Preserve denominator derivative below.
            total=c.sum(); value=-q/(total*soft)
            dq=4*c/3;dq[:-1]+=c[1:]/3;dq[1:]+=c[:-1]/3
            dc=-dq/(total*soft) + q/(total*total*soft) + q*weights/(total*soft*soft)
            score=float(q/(total*maximum))
        else:
            value=2*soft; dc=2*weights;score=2*maximum
        # d(f*f)/df is twice convolution with the reversed f, valid slice.
        grad=2*irfft(rfft(dc,self.size)*rfft(f[::-1],self.size),self.size)[n-1:2*n-1]/n
        if self.kind=='third':grad-=grad.mean()
        else:grad=(grad-np.dot(grad,f)/n)/np.mean(x)
        self.calls+=1
        # Do not turn roundoff at a copied incumbent into an improvement event.
        better=(score>self.keep['score']+1e-12 if self.kind=='second' else score<self.keep['score']-1e-12)
        if better:
            self.keep['score']=score;self.keep['f']=f.copy()
            self.log(dict(event='improvement',kind=self.kind,n=n,beta=self.beta,score=score,calls=self.calls))
        if self.calls%100==0:self.log(dict(event='progress',kind=self.kind,n=n,beta=self.beta,calls=self.calls,score=score,best=self.keep['score']))
        # Scale the objective so the quasi-Newton initial metric is sensible.
        return float(value*n),grad*n

def check_gradient():
    rng=np.random.default_rng(123)
    for kind in ['first','second','third']:
        x=rng.normal(1,3,17) if kind=='third' else rng.uniform(.2,2,17)
        keep=dict(score=-1 if kind=='second' else float('inf'))
        obj=Objective(x,kind,20,keep,time.monotonic()+20,lambda x:None)
        v,g=obj(x);eps=1e-6
        numeric=[]
        for j in range(len(x)):
            d=np.zeros(len(x));d[j]=eps
            numeric.append((obj(x+d)[0]-obj(x-d)[0])/(2*eps))
        err=float(np.max(np.abs(g-numeric)))
        assert err<1e-7,(kind,err)
        print(json.dumps(dict(gradient=kind,max_abs_error=err)),flush=True)
    h=rng.uniform(.1,.9,17);v,g,_=overlap_objective(h,20)
    numeric=[]
    for j in range(len(h)):
        d=np.zeros(len(h));d[j]=1e-6
        numeric.append((overlap_objective(h+d,20)[0]-overlap_objective(h-d,20)[0])/2e-6)
    err=float(np.max(np.abs(g-numeric)));assert err<1e-7,err
    print(json.dumps(dict(gradient='overlap',max_abs_error=err)),flush=True)

def overlap_objective(h,beta):
    n=len(h);size=next_fast_len(2*n-1)
    c=irfft(rfft(h,size)*rfft((1-h)[::-1],size),size)[:2*n-1]*2/n
    maximum=float(c.max());w=np.exp(beta*(c-maximum));total=w.sum();w/=total
    grad=(irfft(rfft(w,size)*rfft(1-h,size),size)[n-1:2*n-1]
          -irfft(rfft(w,size)*rfft(h[::-1],size),size)[n-1:2*n-1][::-1])*2/n
    return maximum+np.log(total)/beta,grad,maximum

def project_half(x):
    # Euclidean projection onto [0,1]^n with sum exactly n/2 to roundoff.
    low=float(x.min()-1);high=float(x.max())
    for _ in range(45):
        mid=(low+high)/2
        if np.clip(x-mid,0,1).sum()>len(x)/2:low=mid
        else:high=mid
    h=np.clip(x-(low+high)/2,0,1)
    active=(h>0)&(h<1)
    h[active]+=(len(h)/2-h.sum())/active.sum()
    return h

def run_overlap(args):
    source=read_json(ROOT/'data'/'erdos-min-overlap-best.json')[0]
    original=np.array(source['data']['values'],float)
    h=project_half(np.repeat(original,args.repeat))
    out=ROOT/'results'/args.name;out.mkdir(parents=True,exist_ok=True)
    start=time.monotonic();deadline=start+args.seconds
    handle=(out/'decisions.jsonl').open('a')
    def log(row):
        row['elapsed_seconds']=time.monotonic()-start
        handle.write(json.dumps(row)+'\n');handle.flush();print(json.dumps(row),flush=True)
    baseline=overlap_objective(h,1)[2];best=baseline;best_h=h.copy()
    log(dict(event='start',kind='overlap',source_solution_id=source['id'],source_agent=source['agentName'],baseline=baseline,args=vars(args)))
    if args.noise:h=project_half(h+np.random.default_rng(args.random_seed).normal(0,args.noise,len(h)))
    for beta in map(float,args.betas.split(',')):
        step=.01;previous_h=previous_g=None
        for iteration in range(args.iterations):
            if time.monotonic()>deadline:break
            value,g,score=overlap_objective(h,beta);g*=len(h)
            if score<best-1e-12:
                best=score;best_h=h.copy();log(dict(event='improvement',beta=beta,iteration=iteration,score=score))
            if previous_h is not None:
                s=h-previous_h;y=g-previous_g;sy=np.dot(s,y)
                if sy>1e-30:step=float(np.clip(np.dot(s,s)/sy,1e-9,10))
            previous_h=h.copy();previous_g=g.copy()
            accepted=False
            for trial in range(24):
                candidate=project_half(h-step*g);delta=candidate-h
                if np.max(np.abs(delta))<1e-14:break
                next_value,_,_=overlap_objective(candidate,beta)
                if next_value<=value+1e-4*np.dot(g,delta)/len(h):
                    accepted=True;h=candidate;break
                step*=.5
            if iteration%100==0:log(dict(event='progress',beta=beta,iteration=iteration,score=score,best=best,step=step))
            if not accepted:break
        log(dict(event='stage_end',beta=beta,iterations=iteration+1,best=best))
        if time.monotonic()>deadline:break
    np.save(out/'best.npy',best_h)
    (out/'candidate.json').write_text(json.dumps({'values':best_h.tolist()},separators=(',',':')))
    summary=dict(kind='overlap',baseline=baseline,best_fft_score=best,source_solution_id=source['id'],elapsed_seconds=time.monotonic()-start,officially_verified=False,publicly_submitted=False)
    (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');log(dict(event='finish',**summary))

def run_c42(args):
    """Public fixed-profile perturbation; finite-n diagnostics are NOT certificates."""
    import importlib.util
    from fractions import Fraction
    spec=importlib.util.spec_from_file_location('c42_finite',ROOT/'data'/'c42'/'finite_n.py')
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    source=read_json(ROOT/'data'/'c42'/'baseline-certificate.json')
    alpha=complex(*map(float,source['alpha']));eta=np.array([complex(*map(float,z)) for z in source['eta']])
    tau=Fraction(args.cutoff);C=args.bound;M=len(eta)
    assert Fraction(1,5)<tau<Fraction(1,2)
    scale=(C-1e-9)/max(abs(1-alpha),abs(eta).max())
    a=1-(1-alpha)*scale;e=eta*scale
    cert=dict(constant='C_42',C=str(C),alpha=[format(a.real,'.17f'),format(a.imag,'.17f')],
              breaks=[str(tau+Fraction(k,M)*(1-2*tau)) for k in range(M+1)],
              eta=[[format(z.real,'.17f'),format(z.imag,'.17f')] for z in e],
              provenance=dict(source='DennisKoshta/c42-turan-power-sum-certificate@8e1797e0d35454ad18cfebbbeb5858bf202b0748',
              method='lower cutoff; uniformly contract early and middle power sums',certified=False))
    out=ROOT/'results'/args.name;out.mkdir(parents=True,exist_ok=True)
    (out/'candidate.json').write_text(json.dumps(cert,indent=2)+'\n')
    start=time.monotonic();rows=[]
    for n in [65536,262144,1048576]:
        baseline=module.Model(n,Fraction(41,200),M).evaluate(alpha,eta)[2]
        value=module.Model(n,tau,M).evaluate(a,e)[2]
        row=dict(event='finite_n_diagnostic',n=n,baseline=baseline,candidate=value,difference=value-baseline,
                 heuristic_continuum_estimate=.6873358061865166+value-baseline,
                 elapsed_seconds=time.monotonic()-start,certified=False)
        rows.append(row);print(json.dumps(row),flush=True)
    (out/'decisions.jsonl').write_text(''.join(json.dumps(row)+'\n' for row in rows))
    (out/'summary.json').write_text(json.dumps(dict(kind='c42',claimed_bound=C,cutoff=str(tau),certified=False,publicly_submitted=False,
        source='Dennis Koshta, public MIT certificate',next_step='Arb checker plus independent verifier; do not call finite-n extrapolation a proof'),indent=2)+'\n')

def run(args):
    kind=args.kind;slug=kind+'-autocorrelation-inequality'
    seed=read_json(ROOT/'data'/f'{slug}-best.json')[args.seed_index]
    if args.seed_file:
        import gzip
        path=Path(args.seed_file)
        raw=gzip.decompress(path.read_bytes()) if path.suffix=='.gz' else path.read_bytes()
        numbers=list(map(int,raw.decode().split()))
        if numbers[0]==len(numbers)-1:numbers=numbers[1:]
        seed=dict(id=args.seed_file,agentName='Mats Heming Julner / intervalplace; L. S. Cumes construction',score=fast_score(slug,numbers),data=dict(values=numbers))
    f=np.array(seed['data']['values'],float);f/=f.mean()
    if args.repeat>1:f=np.repeat(f,args.repeat)
    assert len(f)<=2_000_000
    rng=np.random.default_rng(args.random_seed)
    if args.noise:
        f += args.noise*rng.normal(size=len(f))
        if kind!='third':f=np.maximum(f,0)
        f/=f.mean()
    out=ROOT/'results'/args.name;out.mkdir(parents=True,exist_ok=True)
    start=time.monotonic();deadline=start+args.seconds
    logfile=(out/'decisions.jsonl').open('a')
    def log(row):
        row['elapsed_seconds']=time.monotonic()-start
        logfile.write(json.dumps(row)+'\n');logfile.flush()
        print(json.dumps(row),flush=True)
    initial=fast_score(slug,f)
    if kind=='third':
        from scipy.signal import fftconvolve
        initial=float(2*len(f)*np.abs(fftconvolve(f,f)).max()/f.sum()**2)
    keep=dict(score=initial,f=f.copy())
    log(dict(event='start',kind=kind,source_solution_id=seed['id'],source_agent=seed['agentName'],
             source_advertised=seed['score'],baseline=initial,n=len(f),args=vars(args)))
    x=f.copy()
    if args.precondition:
        from scipy.ndimage import uniform_filter1d
        spectrum=np.abs(rfft(f))
        smooth=np.sqrt(uniform_filter1d(spectrum*spectrum,size=17,mode='wrap'))
        median=float(np.median(smooth));multiplier=np.clip(median/np.maximum(smooth,.01*median),.001,100)
        multiplier[0]=0
        def transform(u):return irfft(rfft(u)*multiplier,n=len(f))
        base=f.copy();x=np.zeros_like(f)
        log(dict(event='preconditioner',method='smoothed Fourier inverse magnitude; positive clipping; exact chain rule',window=17,clip=[.001,100]))
    for beta in map(float,args.betas.split(',')):
        obj=Objective(f,kind,beta,keep,deadline,log)
        if args.precondition:
            def objective(u):
                raw=base+transform(u);v,g=obj(np.maximum(raw,0))
                return v,transform(g*(raw>0))
        else:objective=obj
        try:
            res=minimize(objective,x,method='L-BFGS-B',jac=True,
                         bounds=None if kind=='third' or args.precondition else [(0,None)]*len(x),
                         options=dict(maxiter=args.iterations,maxcor=8,ftol=1e-16,gtol=1e-13,maxls=25))
            x=res.x
            log(dict(event='stage_end',beta=beta,iterations=int(res.nit),message=str(res.message),best=keep['score']))
        except TimeoutError:
            log(dict(event='time_budget_exhausted',beta=beta,best=keep['score']));break
        np.save(out/'best.npy',keep['f'])
    np.save(out/'best.npy',keep['f'])
    (out/'candidate.json').write_text(json.dumps({'values':keep['f'].tolist()},separators=(',',':')))
    summary=dict(kind=kind,baseline=initial,best_fft_score=keep['score'],source_score=seed['score'],
                 source_solution_id=seed['id'],elapsed_seconds=time.monotonic()-start,
                 officially_verified=False,publicly_submitted=False)
    (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');log(dict(event='finish',**summary))

def run_deletion(args):
    """Exact set packing in a bounded repair neighbourhood of a public code."""
    from ortools.sat.python import cp_model
    slug='two-deletion-code-16'
    source=read_json(ROOT/'data'/f'{slug}-best.json')[0]
    words=source['data']['words'];seed={int(w,2) for w in words}
    arrays=np.load(ROOT/'data'/'deletion-descendants.npz')
    d=arrays['descendants'];counts=arrays['conflict_counts']
    frozen=set()
    if args.window_start>=0:
        low=args.window_start;high=low+args.window_width
        frozen={w for w in seed if not low<=w.bit_count()<high}
        owner=np.zeros(16384,dtype=bool)
        for w in frozen:owner[d[w]]=True
        eligible=(~owner[d]).all(axis=1)
        weights=np.array([int(w).bit_count() for w in range(65536)])
        pool=np.flatnonzero(eligible & (weights>=low) & (weights<high))
    else:pool=np.flatnonzero(counts<=args.max_conflicts)
    assert (seed-frozen).issubset(set(pool.tolist()))
    out=ROOT/'results'/args.name;out.mkdir(parents=True,exist_ok=True)
    start=time.monotonic();handle=(out/'decisions.jsonl').open('a')
    def log(row):
        row['elapsed_seconds']=time.monotonic()-start
        handle.write(json.dumps(row)+'\n');handle.flush();print(json.dumps(row),flush=True)
    log(dict(event='start',kind='deletion',source_solution_id=source['id'],source_agent=source['agentName'],
             baseline=len(seed),max_conflicts=args.max_conflicts,drop_budget=args.drops,pool_size=len(pool),random_seed=args.random_seed,
             window_start=args.window_start,window_width=args.window_width,frozen=len(frozen)))
    model=cp_model.CpModel();variables=[model.new_bool_var(f'w{int(w)}') for w in pool]
    row=d[pool];valid=np.ones(row.shape,dtype=bool);valid[:,1:]=row[:,1:]!=row[:,:-1]
    parents=np.broadcast_to(np.arange(len(pool))[:,None],row.shape)[valid];children=row[valid]
    order=np.argsort(children,kind='stable');parents=parents[order];children=children[order]
    boundaries=np.r_[0,np.flatnonzero(children[1:]!=children[:-1])+1,len(children)]
    for left,right in zip(boundaries[:-1],boundaries[1:]):
        if right-left>1:model.add_at_most_one([variables[int(k)] for k in parents[left:right]])
    seedvars=[]
    for w,v in zip(pool,variables):
        in_seed=int(w) in seed;model.add_hint(v,int(in_seed))
        if in_seed:seedvars.append(v)
    if args.drops>=0:model.add(sum(seedvars)>=len(seed)-len(frozen)-args.drops)
    model.maximize(sum(variables))
    log(dict(event='model_ready',variables=len(variables),constraints=len(boundaries)-1))
    solver=cp_model.CpSolver();solver.parameters.max_time_in_seconds=args.seconds
    solver.parameters.num_search_workers=4;solver.parameters.random_seed=args.random_seed
    solver.parameters.log_search_progress=False
    class Callback(cp_model.CpSolverSolutionCallback):
        def __init__(self):super().__init__();self.best=len(seed)
        def on_solution_callback(self):
            score=int(round(self.objective_value))+len(frozen)
            log(dict(event='feasible_code',score=score,bound=self.best_objective_bound+len(frozen)))
            if score>self.best:
                selected=[format(w,'016b') for w in sorted(frozen)]+[format(int(w),'016b') for w,v in zip(pool,variables) if self.value(v)]
                # Independent exact string verifier, reviewed from the public source.
                seen=set()
                for word in selected:
                    descendants={word[:i]+word[i+1:j]+word[j+1:] for i in range(15) for j in range(i+1,16)}
                    assert seen.isdisjoint(descendants);seen.update(descendants)
                self.best=score
                (out/'candidate.json').write_text(json.dumps({'words':selected},indent=2)+'\n')
                log(dict(event='exactly_verified_improvement',score=score))
                if args.stop_after_improvement:self.stop_search()
    cb=Callback();status=solver.solve(model,cb)
    summary=dict(kind='deletion',status=solver.status_name(status),baseline=len(seed),best=cb.best,
                 objective_bound=solver.best_objective_bound+len(frozen),wall_seconds=solver.wall_time,
                 source_solution_id=source['id'],max_conflicts=args.max_conflicts,drop_budget=args.drops,
                 elapsed_seconds=time.monotonic()-start,publicly_submitted=False)
    (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');log(dict(event='finish',**summary))
    return cb.best>len(seed)

def run_snake(args):
    """Exact induced-path tail repair with a finite node/time budget."""
    import sys
    sys.setrecursionlimit(12000)
    slug='snake-in-the-box-13'
    source=read_json(ROOT/'data'/f'{slug}-best.json')[0]
    incumbent=source['data']['actions'];target=len(incumbent)
    out=ROOT/'results'/args.name;out.mkdir(parents=True,exist_ok=True)
    start=time.monotonic();deadline=start+args.seconds
    handle=(out/'decisions.jsonl').open('a');rng=np.random.default_rng(args.random_seed)
    def log(row):
        row['elapsed_seconds']=time.monotonic()-start
        handle.write(json.dumps(row)+'\n');handle.flush();print(json.dumps(row),flush=True)
    neighbours=[[v^(1<<k) for k in range(13)] for v in range(8192)]
    best=target;best_path=None
    log(dict(event='start',kind='snake',baseline=target,source_solution_id=source['id'],source_agent=source['agentName']))
    for reverse in [False,True]:
        seed=incumbent[::-1] if reverse else incumbent
        for cut in [0,2,4,8,16,32,64,128,256,512]:
            if time.monotonic()>deadline:break
            path=list(seed[:len(seed)-cut] if cut else seed)
            used=[False]*8192;adj=[0]*8192;used[0]=True
            for v in neighbours[0]:adj[v]+=1
            current=0
            for action in path:
                current ^=1<<action;used[current]=True
                for v in neighbours[current]:adj[v]+=1
            nodes=0;deepest=len(path);found=False
            def dfs(current):
                nonlocal nodes,deepest,found,best,best_path
                nodes+=1;deepest=max(deepest,len(path))
                if nodes>args.nodes or (nodes%1024==0 and time.monotonic()>deadline):return
                if len(path)>best:
                    best=len(path);best_path=path.copy();found=True
                    (out/'candidate.json').write_text(json.dumps({'actions':best_path})+'\n')
                    log(dict(event='exact_path_improvement',score=best,reverse=reverse,cut=cut,nodes=nodes));return
                choices=[(a,v) for a,v in enumerate(neighbours[current]) if not used[v] and adj[v]==1]
                # Try alternatives to the published tail before replaying it.
                original=seed[len(path)] if len(path)<len(seed) else -1
                choices.sort(key=lambda av: (av[0]==original,sum(not used[w] and adj[w]<=1 for w in neighbours[av[1]])))
                for action,v in choices:
                    used[v]=True;path.append(action)
                    for w in neighbours[v]:adj[w]+=1
                    dfs(v)
                    for w in neighbours[v]:adj[w]-=1
                    path.pop();used[v]=False
                    if found or nodes>args.nodes or time.monotonic()>deadline:return
            dfs(current)
            log(dict(event='tail_search_end',cut=cut,reverse=reverse,nodes=nodes,deepest=deepest,best=best))
            if found:break
        if best_path is not None:break
    summary=dict(kind='snake',baseline=target,best=best,source_solution_id=source['id'],
                 elapsed_seconds=time.monotonic()-start,publicly_submitted=False)
    (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');log(dict(event='finish',**summary))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--check-gradient',action='store_true')
    p.add_argument('--kind',choices=['first','second','third','deletion','snake','overlap','c42'],default='third');p.add_argument('--seconds',type=int,default=120)
    p.add_argument('--iterations',type=int,default=200);p.add_argument('--repeat',type=int,default=1)
    p.add_argument('--seed-index',type=int,default=0);p.add_argument('--random-seed',type=int,default=42)
    p.add_argument('--seed-file')
    p.add_argument('--precondition',action='store_true')
    p.add_argument('--noise',type=float,default=0);p.add_argument('--betas',default='1000,10000,100000,1000000,10000000')
    p.add_argument('--name',default='pilot')
    p.add_argument('--max-conflicts',type=int,default=4);p.add_argument('--drops',type=int,default=12)
    p.add_argument('--window-start',type=int,default=-1);p.add_argument('--window-width',type=int,default=2)
    p.add_argument('--window-all',action='store_true')
    p.add_argument('--nodes',type=int,default=500000)
    p.add_argument('--cutoff',default='51/250');p.add_argument('--bound',type=float,default=.68733)
    p.add_argument('--stop-after-improvement',action='store_true');args=p.parse_args()
    if args.check_gradient:check_gradient()
    elif args.kind=='snake':run_snake(args)
    elif args.kind=='overlap':run_overlap(args)
    elif args.kind=='c42':run_c42(args)
    elif args.kind=='deletion':
        if args.window_all:
            base_name=args.name
            for low in range(17-args.window_width+1):
                args.window_start=low;args.name=f'{base_name}-w{low}'
                improved=run_deletion(args)
                if improved and args.stop_after_improvement:break
        else:run_deletion(args)
    else:run(args)
