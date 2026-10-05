"""Regenerate P2D figures from committed compact campaign data."""
import csv
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
DATA=ROOT/'results/phase2/p2d'
OUT=ROOT/'report/figures/phase2'
OUT.mkdir(parents=True,exist_ok=True)
summary=json.loads((DATA/'summary.json').read_text())
results={(r['case'],r['mode'],r['policy']):r['metrics'] for r in summary['results']}
policies=['cyclic','information','p2b_uniform']
labels={'cyclic':'Cyclic','information':'Information','p2b_uniform':'All links, 1 Hz'}
colors={'cyclic':'#2474a6','information':'#d56b25','p2b_uniform':'#6e6e6e'}
plt.rcParams.update({'font.size':9,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42})

def save(fig,name):
    temporary = OUT/f'p2d_{name}.tmp.pdf'
    fig.savefig(temporary,bbox_inches='tight')
    if not temporary.read_bytes().rstrip().endswith(b'%%EOF'):
        raise RuntimeError(f'Incomplete PDF: {temporary}')
    temporary.replace(OUT/f'p2d_{name}.pdf')
    fig.savefig(OUT/f'p2d_{name}.png',dpi=160,bbox_inches='tight')
    plt.close(fig)

cases=['clean','loss30','loss50','network_outage','node_outage','outliers','nlos']
names=['Clean','30% loss','50% loss','Network\noutage','Node 1\noutage','5% positive\noutliers','Link 1-2\nbias']
fig,axs=plt.subplots(1,2,figsize=(10.4,3.9))
for ax,key,title in zip(axs,['relative_shape_rmse_m','worst_node_rmse_m'],['(a) Centered relative shape','(b) Worst node']):
    x=np.arange(len(cases))
    for index,policy in enumerate(policies):
        vals=[results[c,'naive',policy][key] for c in cases]
        ax.errorbar(x+(index-1)*.16,[v['mean'] for v in vals],yerr=[v['std'] for v in vals],
                    fmt='o',ms=4,capsize=2,color=colors[policy],label=labels[policy])
    ax.set_xticks(x,names,rotation=35,ha='right');ax.set_ylabel('Whole-run RMSE [m]');ax.set_title(title)
    ax.grid(axis='y',alpha=.2)
axs[0].legend(fontsize=8);fig.tight_layout();save(fig,'reliability')

fig,axs=plt.subplots(1,2,figsize=(10.4,3.6))
for ax,case,title in zip(axs,['outliers','nlos'],['(a) 5% positive outliers, +1 to +3 m','(b) Link 1-2 bias, +0.75 m from 20 to 40 s']):
    for pi,p in enumerate(policies):
        for mi,mode in enumerate(['naive','gated']):
            v=results[case,mode,p]['relative_shape_rmse_m']
            ax.errorbar(pi+(mi-.5)*.18,v['mean'],yerr=v['std'],fmt='o' if mi==0 else 's',
                        color=colors[p],capsize=3,label=mode.capitalize() if pi==0 else None)
    ax.set_xticks(range(3),[labels[p] for p in policies]);ax.set_ylabel('Shape RMSE [m]');ax.set_title(title)
    ax.grid(axis='y',alpha=.2);ax.legend(fontsize=8)
fig.tight_layout();save(fig,'gating')

profiles=list(csv.DictReader((DATA/'time_profile.csv').open()))
fig,axs=plt.subplots(1,3,figsize=(10.4,3.2),sharey=True)
for ax,case,mode,title in zip(axs,['network_outage','nlos','nlos'],['naive','naive','gated'],
                            ['(a) Network outage, naive','(b) Biased link, naive','(c) Biased link, gated']):
    for p in ['cyclic','information']:
        group=[r for r in profiles if (r['case'],r['mode'],r['policy'])==(case,mode,p)]
        t=np.array([float(r['time_s']) for r in group]);mean=np.array([float(r['mean_shape_rms_error_m']) for r in group]);std=np.array([float(r['std_shape_rms_error_m']) for r in group])
        ax.plot(t,mean,color=colors[p],label=labels[p]);ax.fill_between(t,np.maximum(0,mean-std),mean+std,color=colors[p],alpha=.14)
    ax.axvspan(20,30 if case=='network_outage' else 40,color='black',alpha=.07)
    ax.set_title(title);ax.set_xlabel('Time [s]');ax.grid(alpha=.2)
axs[0].set_ylabel('Instantaneous shape RMS error [m]');axs[0].legend(fontsize=8)
fig.tight_layout();save(fig,'time_profiles')

fig,axs=plt.subplots(1,2,figsize=(10.4,3.7))
comps={(r['case'],r['mode'],r['policy'],r['reference']):r for r in summary['comparisons']}
for ax,mode,cs in zip(axs,['naive','gated'],[cases,['clean','network_outage','outliers','nlos']]):
    x=np.arange(len(cs));v=[comps[c,mode,'information','same_case_cyclic']['metrics']['relative_shape_rmse_m'] for c in cs]
    m=np.array([r['mean'] for r in v]);lo=np.array([r['ci95'][0] for r in v]);hi=np.array([r['ci95'][1] for r in v])
    ax.errorbar(x,m,yerr=[m-lo,hi-m],fmt='o',capsize=3,color=colors['information'])
    ax.axhline(0,color='black',lw=.8);ax.set_xticks(x,[names[cases.index(c)] for c in cs],rotation=35,ha='right')
    ax.set_ylabel('Information minus cyclic shape RMSE [m]');ax.set_title(mode.capitalize()+' EKF: paired mean and bootstrap 95% interval');ax.grid(axis='y',alpha=.2)
fig.tight_layout();save(fig,'paired_differences')
print('Generated four P2D PDF figures')
