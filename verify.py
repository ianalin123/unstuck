"""Frozen public-verifier replay plus independently derived FFT scores."""
from pathlib import Path
import argparse
import hashlib
import json
import gzip
from fractions import Fraction
import time
import numpy as np
from scipy.signal import fftconvolve

ROOT = Path(__file__).resolve().parent

def read_json(path):
    path=Path(path)
    if path.suffix=='.gz':return json.loads(gzip.decompress(path.read_bytes()))
    if path.exists():return json.loads(path.read_text())
    return json.loads(gzip.decompress(path.with_suffix(path.suffix+'.gz').read_bytes()))

def problem(slug):
    return read_json(ROOT / 'data' / f'problems_{slug}.json')

def official(slug, data):
    # Only explicitly reviewed downloaded verifier sources can execute here.
    assert slug in {'first-autocorrelation-inequality', 'second-autocorrelation-inequality',
                    'third-autocorrelation-inequality', 'flat-polynomials', 'hadamard-det-51',
                    'two-deletion-code-16', 'snake-in-the-box-13', 'spencer-discrepancy', 'erdos-min-overlap'}
    scope = {}
    exec(compile(problem(slug)['verifier'], f'public_verifier:{slug}', 'exec'), scope)
    return float(scope['evaluate'](data))

def fast_score(slug, values):
    f = np.asarray(values, dtype=np.float64)
    assert f.ndim == 1 and np.isfinite(f).all() and abs(f.sum()) > 0
    c = fftconvolve(f, f)
    if slug.startswith('second'):
        assert (f >= 0).all()
        return float((2*np.dot(c,c) + np.dot(c[:-1],c[1:])) / (3*c.sum()*c.max()))
    if slug.startswith('first'):
        assert (f >= 0).all()
    return float(2*len(f)*abs(c.max()) / f.sum()**2)

def replay():
    out = []
    for slug in ['second-autocorrelation-inequality', 'third-autocorrelation-inequality',
                 'first-autocorrelation-inequality', 'flat-polynomials']:
        seed = read_json(ROOT/'data'/f'{slug}-best.json')[0]
        start = time.perf_counter()
        value = official(slug, seed['data'])
        row = dict(slug=slug, source_solution_id=seed['id'], source_agent=seed['agentName'],
                   advertised_score=seed['score'], official_score=value,
                   difference=value-seed['score'], elapsed_seconds=time.perf_counter()-start,
                   verifier_sha256=hashlib.sha256(problem(slug)['verifier'].encode()).hexdigest())
        if slug != 'flat-polynomials': row['independent_fft_score'] = fast_score(slug, seed['data']['values'])
        if slug.startswith('third'):
            f=np.asarray(seed['data']['values'],dtype=float)
            row['mathematical_max_abs_score']=float(2*len(f)*np.abs(fftconvolve(f,f)).max()/f.sum()**2)
            row['warning']='Official checker uses abs(max(c)), but the L-infinity norm requires max(abs(c)); do not treat its advertised leader as a mathematical upper bound.'
        out.append(row)
        print(json.dumps(row), flush=True)
        (ROOT/'results').mkdir(exist_ok=True)
        (ROOT/'results'/'baseline-replay.json').write_text(json.dumps(out, indent=2)+'\n')

def exact_first(path, save=False):
    """Every convolution coefficient, by carry-free integer polynomial squaring.

    Independent implementation from the defining polynomial identity. No float
    decides feasibility, the maximum or whether the rational bound improves.
    """
    path=Path(path);started=time.perf_counter()
    if path.suffix=='.npy':
        f=np.load(path);assert f.ndim==1 and np.isfinite(f).all() and (f>=0).all() and f.max()>0
        maximum=float(f.max())
        values=[int(round(float(v/maximum)*2**50)) for v in f]
    else:
        raw=gzip.decompress(path.read_bytes()) if path.suffix=='.gz' else path.read_bytes()
        values=list(map(int,raw.decode().split()))
    n=len(values);assert n>0 and min(values)>=0 and sum(values)>0
    width=(int(n*max(values)**2).bit_length()+7)//8
    assert 2**(8*width)>n*max(values)**2
    packed=b''.join(v.to_bytes(width,'little') for v in values)
    polynomial=int.from_bytes(packed,'little')
    product=(polynomial*polynomial).to_bytes((2*n-1)*width,'little')
    peak=-1;peak_index=-1;total=0
    for k in range(2*n-1):
        coefficient=int.from_bytes(product[k*width:(k+1)*width],'little');total+=coefficient
        if coefficient>peak:peak=coefficient;peak_index=k
    assert total==sum(values)**2
    rng=np.random.default_rng(20261115)
    sample={peak_index,0,n-1,2*n-2,*map(int,rng.integers(0,2*n-1,size=23))}
    for k in sample:
        direct=sum(values[i]*values[k-i] for i in range(max(0,k-n+1),min(n,k+1)))
        assert direct==int.from_bytes(product[k*width:(k+1)*width],'little')
    value=Fraction(2*n*peak,sum(values)**2)
    record=Fraction(735013322210029338350039415555883008,489139897323017668625465543540606041)
    raw=('\n'.join(map(str,values))+'\n').encode()
    if save:
        target=path.parent/'candidate-integers.txt.gz';target.write_bytes(gzip.compress(raw,mtime=0))
    row=dict(n=n,method='carry-free polynomial square; all coefficients; exact total and direct sample checks',
             score_fraction=str(value),score_float=float(value),peak_index=peak_index,peak_coefficient=str(peak),sum_heights=str(sum(values)),
             comparison_fraction=str(record),comparison_source='intervalplace/frontier C1_262144_1.5026648332.txt, checked 2026-11-15',
             beats_comparison_exactly=value<record,improvement_fraction=str(record-value),
             integer_witness_sha256=hashlib.sha256(raw).hexdigest(),elapsed_seconds=time.perf_counter()-started,
             publicly_submitted=False)
    if save:(path.parent/'exact-verification.json').write_text(json.dumps(row,indent=2)+'\n')
    print(json.dumps(row),flush=True)

def plot_progress():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    reference=float(Fraction(735013322210029338350039415555883008,489139897323017668625465543540606041))
    fig,ax=plt.subplots(figsize=(10,5.5));comparison=[]
    for name,label,color in [('first-current-sharp','Original grid, direct refinement','#3678b5'),
                             ('first-fourier','Original grid, Fourier preconditioner','#d78923'),
                             ('first-current-finer','Doubled grid, direct refinement','#219979')]:
        folder=ROOT/'results'/name
        rows=[json.loads(s) for s in (folder/'decisions.jsonl').read_text().splitlines()]
        start=rows[0];events=[r for r in rows if r['event']=='improvement']
        xs=[0]+[r['elapsed_seconds'] for r in events];ys=[start['baseline']]+[r['score'] for r in events]
        xs.append(rows[-1]['elapsed_seconds']);ys.append(ys[-1])
        ax.step(xs,(np.array(ys)-reference)*1e10,where='post',color=color,label=label,lw=2)
        exact_path=folder/'exact-verification.json'
        exact=read_json(exact_path) if exact_path.exists() else None
        if exact:ax.scatter(xs[-1],(exact['score_float']-reference)*1e10,color=color,edgecolor='white',s=65,zorder=4)
        comparison.append(dict(run=name,source_baseline=start['baseline'],first_accepted_improvement_seconds=events[0]['elapsed_seconds'] if events else None,
                               last_logged_seconds=xs[-1],best_numerical_score=ys[-1],exact_endpoint=exact['score_fraction'] if exact else None))
    ax.axhline(0,color='#687383',ls='--',lw=1,label='Newest public witness checked')
    ax.set(xlabel='Wall time within each run (seconds)',ylabel='Bound minus public reference (× 10⁻¹⁰)\nLower is better',title='Useful search gains can arrive after several minutes')
    ax.spines[['top','right']].set_visible(False);ax.grid(axis='y',alpha=.2)
    ax.legend(frameon=False,fontsize=9,loc='upper right')
    fig.text(.08,.025,'Cumes / Julner warm starts credited. Runs overlapped on one laptop; these settings differ beyond grid size.\nCircles show exact integer endpoint checks. Local evidence; no public record acceptance.',fontsize=8,color='#596373')
    fig.tight_layout(rect=(0,.09,1,1));fig.savefig(ROOT/'results'/'search-progress.png',dpi=180)
    svg=ROOT/'results'/'search-progress.svg'
    fig.savefig(svg);plt.close(fig)
    svg.write_text('\n'.join(line.rstrip() for line in svg.read_text().splitlines())+'\n')
    (ROOT/'results'/'search-comparison.json').write_text(json.dumps(comparison,indent=2)+'\n')
    print(json.dumps(comparison,indent=2))

if __name__ == '__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--candidate');parser.add_argument('--slug');parser.add_argument('--exact-first');parser.add_argument('--save-exact',action='store_true');parser.add_argument('--plot-progress',action='store_true')
    args=parser.parse_args()
    if args.plot_progress:plot_progress()
    elif args.exact_first:exact_first(args.exact_first,args.save_exact)
    elif args.candidate:
        data=read_json(Path(args.candidate))
        print(json.dumps({'score':official(args.slug,data),'slug':args.slug}))
    else: replay()
