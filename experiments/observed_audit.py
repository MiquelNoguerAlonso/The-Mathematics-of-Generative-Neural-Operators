"""Observed-sample audit of bounded spectral conditional generators.

The audit function receives observed calibration fields, independent aggregation
contexts, model parameters, and stated support/regularity constants. It never
calls target_mean. Truth is used only by the sampler and a separate diagnostic.
"""
from pathlib import Path
import json, platform
import numpy as np
from scipy.special import zeta
from numpy.polynomial.legendre import leggauss
import torch
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1]
SEEDS=[481,482,483]
BUDGETS=[4096,16384,65536]
BINS=[8,12,16]
M=32
SIGMA=.1
DELTA=.01

def target_mean(x):
    return .7*np.tanh(.8*x)+.1*np.sin(3*x)

def observed_fields(rng,n,m):
    x=rng.uniform(-1,1,n)
    out=np.empty((n,m));out[:,0]=x
    out[:,1:]=(target_mean(x)[:,None]+SIGMA*rng.uniform(-1,1,(n,m-1)))/np.arange(2,m+1)
    return out

def predict(x,model):
    return np.tanh(np.asarray(x)[...,None]*np.asarray(model['slopes']))@np.asarray(model['weights'])

def train(seed):
    rng=np.random.default_rng(seed)
    data=observed_fields(rng,16384,8)
    torch.manual_seed(seed)
    slopes=torch.nn.Parameter(torch.randn(12)*.4)
    logits=torch.nn.Parameter(torch.zeros(12))
    opt=torch.optim.Adam([slopes,logits],lr=.025)
    for step in range(1200):
        i=rng.integers(len(data),size=1024);j=rng.integers(1,8,size=1024)
        x=torch.from_numpy(data[i,0]);y=torch.from_numpy(data[i,j]*(j+1))
        pred=torch.tanh(x[:,None]*slopes)@logits.softmax(0)
        loss=((pred-y)**2).mean();opt.zero_grad();loss.backward();opt.step()
        with torch.no_grad(): slopes.clamp_(-1.3,1.3)
    a=slopes.detach().numpy();w=logits.detach().softmax(0).numpy()
    return {'seed':seed,'slopes':a.tolist(),'weights':w.tolist(),'context_lipschitz':float(w@abs(a)),
            'sigma':SIGMA,'parameters':24,'training_noisy_mse':float(loss.detach())}

def empirical_uniform_w2(samples,center,halfwidth):
    z=np.sort(samples);n=len(z);q=(np.arange(n)+.5)/n
    return np.sqrt(np.mean((z-center-halfwidth*(2*q-1))**2)+halfwidth**2/(3*n*n))

def audit(calibration,aggregation_x,models,bins,delta):
    """Target-independent audit; fixed constants are model-class assumptions."""
    n,m=calibration.shape;J=m-1;K=len(models);diameter=2.2;target_lip=.86
    edges=np.linspace(-1,1,bins+1);centers=(edges[1:]+edges[:-1])/2
    ci=np.minimum(np.searchsorted(edges,calibration[:,0],side='right')-1,bins-1)
    ai=np.minimum(np.searchsorted(edges,aggregation_x,side='right')-1,bins-1)
    assert np.all((ci>=0)&(ai.min()>=0))
    counts=np.bincount(ci,minlength=bins)
    kappa=np.minimum(1,np.sqrt(np.log(4*J*bins/delta)/(2*np.maximum(counts,1))))
    rows=[]
    for model in models:
        radius=np.zeros((J,bins));error=np.zeros(m)
        for j in range(1,m):
            for c in range(bins):
                if counts[c]==0: radius[j-1,c]=diameter/(j+1);continue
                # Normalize details before exact empirical-to-uniform transport.
                w=empirical_uniform_w2(calibration[ci==c,j]*(j+1),predict(centers[c],model),SIGMA)
                radius[j-1,c]=min(diameter,(target_lip+model['context_lipschitz'])*(2/bins)+w+diameter*np.sqrt(kappa[c]))/(j+1)
            d=diameter/(j+1)
            error[j]=min(d*d,np.mean(radius[j-1,ai]**2)+d*d*np.sqrt(np.log(2*J*K/delta)/(2*len(ai))))
        radii=np.zeros(m)
        for j in range(1,m):
            radii[j]=np.hypot(radii[j-1],np.sqrt(error[j])+model['context_lipschitz']/(j+1)*radii[j-1])
        cuts=[]
        for cut in [4,8,16,32]:
            tail=(.8**2+SIGMA**2/3)*zeta(2,cut+1)
            dependency_squared=float(np.sum(error[:cut]))
            cuts.append({'m':cut,'projected_certificate':float(np.sqrt(dependency_squared)),
                         'generic_prefix_projected_certificate':float(radii[cut-1]),
                         'generic_prefix_full_law_certificate':float(np.sqrt(tail+radii[cut-1]**2)),
                         'tail_bound':float(np.sqrt(tail)),
                         'full_law_certificate':float(np.sqrt(tail+dependency_squared))})
        rows.append({'model':model['seed'],'cuts':cuts,'conditional_error_squared_upper':error.tolist(),
                     'min_count':int(counts.min()),'empty_cells':int(np.sum(counts==0))})
    return rows

def main():
    torch.set_default_dtype(torch.float64);torch.set_num_threads(1);torch.use_deterministic_algorithms(True)
    models=[train(s) for s in SEEDS]
    models.append({'seed':'no_context','slopes':[0.],'weights':[1.],'context_lipschitz':0.,'sigma':SIGMA,'parameters':0})
    cal=observed_fields(np.random.default_rng(7401),max(BUDGETS),M)
    agg=observed_fields(np.random.default_rng(7402),max(BUDGETS),1)[:,0]
    results=[]
    for n,bins in zip(BUDGETS,BINS):
        rows=audit(cal[:n],agg[:n],models,bins,DELTA/len(BUDGETS))
        results.append({'n_calibration':n,'n_aggregation':n,'bins':bins,'rows':rows})
        print(n,[round(r['cuts'][-1]['full_law_certificate'],6) for r in rows],flush=True)
    q,w=leggauss(512);truth=target_mean(q)
    diagnostics=[]
    for model in models:
        mse=float(w@((truth-predict(q,model))**2)/2)
        cuts=[]
        for cut in [4,8,16,32]:
            projected=mse*sum(1/j**2 for j in range(2,cut+1))
            tail=(float(w@(truth**2)/2)+SIGMA**2/3)*zeta(2,cut+1)
            cuts.append({'m':cut,'shared_noise_projected_rms':float(np.sqrt(projected)),
                         'shared_noise_full_coupling_upper':float(np.sqrt(projected+tail))})
        diagnostics.append({'model':model['seed'],'conditional_mean_mse_quadrature':mse,'cuts':cuts})
    # Meaningful boundary checks of the empirical quantile integral.
    assert abs(empirical_uniform_w2(np.array([0.]),0.,1.)**2-1/3)<1e-14
    assert abs(empirical_uniform_w2(np.array([-1.,1.]),0.,0.)-1)<1e-14
    for b in results:
        for row,diag in zip(b['rows'],diagnostics):
            assert all(a['shared_noise_full_coupling_upper']<=c['full_law_certificate'] for a,c in zip(diag['cuts'],row['cuts']))
    # A second quadrature order checks the separate truth diagnostic.
    q2,w2=leggauss(1024)
    quadrature_diff=max(abs(float(w2@((target_mean(q2)-predict(q2,mo))**2)/2)-d['conditional_mean_mse_quadrature']) for mo,d in zip(models,diagnostics))
    for model in models[:-1]:
        (ROOT/f"models/observed_generator_{model['seed']}.json").write_text(json.dumps(model,indent=2)+'\n')
    report={'description':'Training on noisy observed fields; conditional distribution audit without target means',
            'confidence':1-DELTA,'joint_scope':'All four fixed candidates, all audited bands, all cutoffs, and all three stated sample budgets',
            'training':{'fields':16384,'modes':8,'steps':1200,'batch_pairs':1024,'width':12,'seeds':SEEDS},
            'audit_uses_target_mean':False,'propagation':'Dependency certificate: the only parent is the exactly sampled first mode; generic prefix comparison also recorded',
            'target_regularization_assumptions':{'feature':'first coefficient','target_feature_lipschitz':.86,'normalized_support':[-1.1,1.1]},
            'budgets':results,'truth_diagnostics':diagnostics,'truth_quadrature_order_difference':quadrature_diff,
            'environment':{'python':platform.python_version(),'numpy':np.__version__,'torch':torch.__version__},
            'arithmetic_scope':'Float64 theorem-bound evaluation; not an interval-arithmetic certificate'}
    (ROOT/'results/observed_audit.json').write_text(json.dumps(report,indent=2)+'\n')
    table=[r'\begin{tabular}{rrrrr}',r'\toprule',r'Fields per split & Cells & Learned bound & No-context bound & Coupling upper \\',r'\midrule']
    for b in results:
        table.append(f"{b['n_calibration']:,} & {b['bins']} & {b['rows'][0]['cuts'][-1]['full_law_certificate']:.4f} & {b['rows'][-1]['cuts'][-1]['full_law_certificate']:.4f} & {diagnostics[0]['cuts'][-1]['shared_noise_full_coupling_upper']:.4f} "+r'\\')
    table += [r'\bottomrule',r'\end{tabular}'];(ROOT/'tables/observed_audit.tex').write_text('\n'.join(table)+'\n')
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    fig,ax=plt.subplots(1,2,figsize=(10.4,3.8),layout='constrained')
    for i in range(4):
        ax[0].plot(BUDGETS,[b['rows'][i]['cuts'][-1]['full_law_certificate'] for b in results],'o-',label=f'Seed {i+1}' if i<3 else 'No context')
    ax[0].set(xscale='log',xlabel='Observed fields per audit split',ylabel='99% full-law upper bound',title='No conditional-mean oracle');ax[0].legend(fontsize=8)
    for i in [0,3]:
        ax[1].plot([4,8,16,32],[r['full_law_certificate'] for r in results[-1]['rows'][i]['cuts']],'o-',label='Learned certificate' if i==0 else 'No-context certificate')
        ax[1].plot([4,8,16,32],[r['shared_noise_full_coupling_upper'] for r in diagnostics[i]['cuts']],':',label='Learned coupling upper' if i==0 else 'No-context coupling upper')
    ax[1].set(xlabel='Retained modes',ylabel='Law / coupling upper bound',title='Resolution and statistical uncertainty');ax[1].legend(fontsize=8)
    fig.savefig(ROOT/'figures/observed_audit.png',dpi=220,facecolor='white');plt.close(fig)
    c=np.linspace(0,.99,100)
    fig,ax=plt.subplots(1,2,figsize=(10,3.5),layout='constrained')
    ax[0].plot(c,1+c,label='Exact aligned-error ratio');ax[0].axhline(1,color='grey',ls='--',label='Ignoring the Gram matrix')
    ax[0].set(xlabel='Overlap c',ylabel='Physical / coefficient squared error',title='Nonorthogonal observation geometry');ax[0].legend(fontsize=8)
    ax[1].semilogy(c,np.sqrt((1+c)/(1-c)));ax[1].set(xlabel='Overlap c',ylabel='Condition number of the encoder',title='Forward accuracy and inverse stability')
    fig.savefig(ROOT/'figures/observation_geometry.png',dpi=220,facecolor='white');plt.close(fig)
    print('truth quadrature difference',quadrature_diff,flush=True)
if __name__=='__main__':main()
