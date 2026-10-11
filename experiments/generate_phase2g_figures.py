"""P2G scientific figures: all headline errors score the same three nodes."""
import csv
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT/'results/phase2/p2g'
OUT = ROOT/'report/figures/phase2'
SUMMARY = json.loads((DATA/'summary.json').read_text())
PROFILE = list(csv.DictReader((DATA/'time_profile.csv').open()))
CASES = ['unanchored', 'exact_100hz', 'noisy_1hz', 'noisy_1hz_outage']
LABELS = ['no continuing reference', 'exact position, 100 Hz', 'noisy position, 1 Hz', 'noisy 1 Hz, outage']
COLORS = ['#555555', '#0072b2', '#009e73', '#d55e00']
plt.rcParams.update({'font.size':9,'pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False})


def save(fig,name):
    OUT.mkdir(parents=True,exist_ok=True)
    fig.savefig(OUT/f'{name}.pdf',bbox_inches='tight')
    fig.savefig(OUT/f'{name}.png',dpi=170,bbox_inches='tight')
    plt.close(fig)


fig, axes = plt.subplots(2,2,figsize=(10,7))
for ax,key,title in zip(axes.flat,['fleet_rmse_m','centroid_rmse_m','relative_shape_rmse_m','worst_node_rmse_m'],
                       ['Absolute fleet (nodes 2–4)','Centroid (nodes 2–4)','Relative shape (nodes 2–4)','Worst node among nodes 2–4']):
    for case,label,color in zip(CASES,LABELS,COLORS):
        rows=[r for r in SUMMARY['results'] if r['case']==case and r['budget']>0]
        x=[r['budget'] for r in rows];y=[r['metrics'][key]['mean'] for r in rows]
        ax.plot(x,y,'o-',color=color,label=label,ms=4)
    if key=='fleet_rmse_m':ax.axhline(.1,color='black',ls=':',label='10 cm absolute target')
    ax.set(xscale='log',xlabel='Attempted UWB exchanges / 60 s',ylabel='Mean RMSE (m)',title=title)
    ax.set_xticks([180,360,720,3600],labels=['180','360','720','3600']);ax.grid(alpha=.2)
axes[0,0].legend(fontsize=7)
fig.tight_layout();save(fig,'p2g_budget_curves')

fig,axes=plt.subplots(1,3,figsize=(12,3.8))
for ax,key,title in zip(axes,['fleet_rms_error_m','shape_rms_error_m','reference_error_m'],
                       ['Absolute fleet error, nodes 2–4','Shape error, nodes 2–4','Reference-node error (reference cases)']):
    for case,label,color in zip(CASES,LABELS,COLORS):
        if key == 'reference_error_m' and case == 'unanchored':
            continue
        rows=[r for r in PROFILE if r['case']==case and int(r['budget'])==360]
        t=np.array([float(r['time_s']) for r in rows]);mean=np.array([float(r['mean_'+key]) for r in rows])
        std=np.array([float(r['std_'+key]) for r in rows]);ax.plot(t,mean,color=color,label=label,lw=1.3)
        ax.fill_between(t,np.maximum(0,mean-std),mean+std,color=color,alpha=.12)
    ax.axvspan(20,30,color='#d55e00',alpha=.08)
    ax.set(xlabel='Time (s)',ylabel='Error (m)',title=title);ax.grid(alpha=.2)
axes[0].legend(fontsize=7)
fig.tight_layout();save(fig,'p2g_time_profiles')

fig,axes=plt.subplots(1,2,figsize=(10,4))
for ax,key,title in zip(axes,['fleet_rmse_m','relative_shape_rmse_m'],
                       ['Paired absolute fleet difference','Paired relative-shape difference']):
    for i,(case,label,color) in enumerate(zip(CASES[1:],LABELS[1:],COLORS[1:])):
        points=[r for r in SUMMARY['comparisons'] if r['case']==case and r['reference']=='unanchored' and r['budget']>0]
        x=np.array([r['budget'] for r in points])*np.exp((i-1)*.055)
        y=np.array([r['metrics'][key]['mean'] for r in points])
        lo=np.array([r['metrics'][key]['ci95'][0] for r in points]);hi=np.array([r['metrics'][key]['ci95'][1] for r in points])
        ax.errorbar(x,y,yerr=[y-lo,hi-y],fmt='o',capsize=3,color=color,label=label)
    ax.axhline(0,color='grey',ls=':');ax.set(xscale='log',xlabel='Matched UWB attempt budget',
        ylabel='Reference case minus unanchored (m)',title=title)
    ax.set_xticks([180,360,720,3600],labels=['180','360','720','3600']);ax.grid(alpha=.2)
axes[0].legend(fontsize=7)
fig.tight_layout();save(fig,'p2g_paired_differences')
