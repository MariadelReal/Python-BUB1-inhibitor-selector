"""Independent verification using archived workbook labels and Maestro subsets.
Run from the workspace containing upload/. Requires openpyxl and both supplied scripts.
Does not recreate human selection judgments or run docking/MMGBSA.
"""
from pathlib import Path
from collections import Counter,defaultdict
import csv,json,math
from openpyxl import load_workbook
from bub1_selection_to_mmgbsa import combine_selector,extract_archive_members,load_csv
from reconstruir_BUB1_61384_a_183 import values_only_ct_headers
U=Path('upload'); O=Path('outputs/bub1_reproduction'); O.mkdir(parents=True,exist_ok=True)
w=load_workbook(U/'Structures_183_to_30(2).xlsx',read_only=True,data_only=True)
def sheet(name):
 rows=iter(w[name].values); h=next(rows);return [dict(zip(h,r)) for r in rows]
master,first=combine_selector(load_csv(U/'BUB1_SELECTOR_183x3_with_compound_ID(2).csv'),sheet('All183_clasificados'),[])
assert [r['drugbank_id'] for r in first]==[r['drugbank_id'] for r in sheet('Top30_presentable')]
by={r['drugbank_id']:r for r in master}; stages={}
for label,fn,n in [('revised30','NewTop30(2).maegz',30),('top10','NewTop10(2).maegz',10),('top3','BUB1_top3_Atorvastatin_Vatalanib_Montelukast_with_models_references(6).maegz',3)]:
 a,_=extract_archive_members(U/fn);assert len(a)==n
 for db,poses in a.items():
  for pdb,r in poses.items():assert math.isclose(r['gscore'],by[db]['score_'+pdb],abs_tol=1e-8,rel_tol=0)
 stages[label]=set(a)
assert stages['top3']<=stages['top10']<=stages['revised30']<=set(by)
fids={r['drugbank_id'] for r in first}
summary={'initial30_exact_order_match_to_workbook':True,'archived_scores_match_selector':True,'nested_subsets':True,'categories':dict(Counter(r['historical_category'] for r in master)),'revised30_categories':dict(Counter(by[x]['historical_category'] for x in stages['revised30'])),'removed':sorted(fids-stages['revised30']),'added':sorted(stages['revised30']-fids),'top10':[by[x]['compound_name'] for x in sorted(stages['top10'])],'top3':[by[x]['compound_name'] for x in sorted(stages['top3'])]}
# MM/GBSA receptor assignment uses per-record source/grid metadata, never row order.
energies={}; receptor=None; duplicates=0
for ct,p in values_only_ct_headers(U/'Structures_MMGBSAMARIJOSE_Nora-Nidia (2)(1).maegz'):
 title=p.get('s_m_title','');db=p.get('s_sd_DRUGBANK_ID','');e=p.get('r_psp_MMGBSA_dG_Bind','')
 metadata=p.get('s_m_Source_File','')+' '+p.get('s_i_glide_gridfile','')
 hits={x for x in ('4QPM','5DMZ','6F7B') if x in metadata.upper()}
 receptor=next(iter(hits)) if len(hits)==1 else None
 if db and e and e!='<>':
  assert receptor, (ct,title)
  key=(db,receptor);value=float(e)
  if key in energies:
   assert math.isclose(energies[key],value,abs_tol=1e-8);duplicates+=1
  energies[key]=value
ids={k[0] for k in energies};assert ids==stages['revised30'];assert len(energies)==90
rank=[]
for db in ids:
 vals=[energies[(db,p)] for p in ('4QPM','5DMZ','6F7B')]
 rank.append({'drugbank_id':db,'name':by[db]['compound_name'],'mean_mmgbsa':sum(vals)/3,**dict(zip(('4QPM','5DMZ','6F7B'),vals))})
rank.sort(key=lambda r:r['mean_mmgbsa'])
for i,r in enumerate(rank,1):r['rank']=i
summary['mmgbsa']={'unique_pairs':len(energies),'compounds':len(ids),'duplicate_pairs_removed':duplicates,'ranking':rank}
with (O/'mmgbsa_recomputed.csv').open('w',newline='') as f:
 wr=csv.DictWriter(f,fieldnames=list(rank[0]));wr.writeheader();wr.writerows(rank)
(O/'selection_verification.json').write_text(json.dumps(summary,indent=2,ensure_ascii=False))
print(json.dumps(summary,indent=2,ensure_ascii=False))
