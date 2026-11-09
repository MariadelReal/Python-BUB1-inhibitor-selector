"""Build supplementary CSV tables from archived exports. Python 3.9+."""
from pathlib import Path
import csv, hashlib, json, math
from statistics import mean

base = Path(__file__).resolve().parent
inputs = base / 'inputs'
output = base / 'tables'
output.mkdir(exist_ok=True)
pdbs = ('4QPM', '5DMZ', '6F7B')
def read(name):
    with (inputs / name).open(encoding='utf-8-sig', newline='') as f:
        return list(csv.DictReader(f))
def save(name, rows, columns=None):
    with (output / name).open('w', encoding='utf-8-sig', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=columns or list(rows[0]), extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)
def check(condition, message):
    if not condition:
        raise ValueError(message)
for name, digest in json.loads((base/'input_hashes.json').read_text()).items():
    check(hashlib.sha256((inputs/name).read_bytes()).hexdigest()==digest, name+' checksum')
trace = read('trace183.csv')
poses = read('poses549.csv')
energy = read('mmgbsa30.csv')
annotations = json.loads((inputs/'annotations.json').read_text(encoding='utf-8'))
verification = json.loads((inputs/'selection_verification.json').read_text(encoding='utf-8'))
ids = {r['drugbank_id'] for r in trace}
check(len(ids)==len(trace)==183, '183 distinct IDs')
pairs = {(r['drugbank_id'],r['pdb']):float(r['glide_gscore']) for r in poses}
check(len(poses)==len(pairs)==549, '549 unique poses')
stage = lambda key: {r['drugbank_id'] for r in trace if r[key]=='1'}
first, revised, ten, three = map(stage,['first30','revised30','top10','top3_MMGBSA_input'])
check(list(map(len,(first,revised,ten,three)))==[30,30,10,3] and three<=ten<=revised, 'historical subsets')
check(first-revised==set(verification['removed']) and revised-first==set(verification['added']), 'substitutions')
for r in trace:
    scores = [pairs[r['drugbank_id'],p] for p in pdbs]
    expected = scores+[max(scores), mean(scores),max(scores)-min(scores)]
    keys = ['score_'+p for p in pdbs]+['worst_gscore','mean_gscore','spread']
    check(all(abs(v-float(r[k]))<1e-8 for v,k in zip(expected,keys)),r['drugbank_id'])
columns = ['drugbank_id','compound_name','historical_category','robustness_global_rank']+keys+['first30','revised30','top10','top3_MMGBSA_input']
save('S1_consensus183.csv',trace,columns)
save('S1b_poses549.csv',poses,['drugbank_id','generic_name','pdb','variant','glide_gscore','entry_id','ct_index'])
save('S2_shortlist30.csv',[r for r in trace if r['drugbank_id'] in revised],columns)
changes=[dict(drugbank_id=r['drugbank_id'],compound_name=r['compound_name'],change='added' if r['drugbank_id'] in revised else 'removed') for r in trace if r['drugbank_id'] in first^revised]
save('S2b_changes.csv',changes)
check(set(annotations)==revised,'annotation IDs')
save('S3_contacts6F7B.csv',[dict(drugbank_id=r['drugbank_id'],compound_name=r['compound_name'],**annotations[r['drugbank_id']]) for r in trace if r['drugbank_id'] in revised])
check(len(energy)==30 and {r['drugbank_id'] for r in energy}==revised,'MM/GBSA membership')
energies=[]
for rank,r in enumerate(sorted(energy,key=lambda x:(mean(float(x[p]) for p in pdbs),x['drugbank_id'])),1):
    avg=mean(float(r[p]) for p in pdbs)
    check(abs(avg-float(r['mean_mmgbsa']))<1e-8 and rank==int(r['rank']),'MM/GBSA rank or mean')
    for p in pdbs:
        energies.append(dict(drugbank_id=r['drugbank_id'],compound_name=r['name'],pdb=p,MMGBSA_kcal_mol=r[p],three_receptor_mean_kcal_mol=avg,rank_of_mean=rank))
save('S4_MMGBSA90.csv',energies)
print('PASS: 183 compounds, 549 poses, 30/10/3 subsets, 6 out + 6 in, 30 diagrams and 90 energies.')
