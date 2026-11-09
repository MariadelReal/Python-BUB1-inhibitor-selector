#!/usr/bin/env python3
"""BUB1: audit existing docking and trace documented selection to MM/GBSA inputs.

Retrospective, NOT the original docking or a learned pharmacological ranker.
Must supply archived 549-pose selector, historical A/B/C/D curation, and the
three archived MAEGZ selection files. Computes original A-first Top30 and
verifies the recorded revised30 -> top10 -> top3 (does NOT invent their ranking).

Python >=3.9, standard library only. The three chosen inputs were documented
for MM/GBSA; this script does not compute MM/GBSA energies.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import sys
from collections import Counter, defaultdict
from pathlib import Path

from reconstruir_BUB1_61384_a_183 import values_only_ct_headers

PDBS = ("4QPM", "5DMZ", "6F7B")
CATS = {"A_prioritario":25, "B_revisar":82, "C_reserva_ATP_like":75, "D_excluir":1}
FINAL = {"DB01076": "Atorvastatin", "DB04879": "Vatalanib", "DB00471": "Montelukast"}
RU = {"DB01908", "DB02336", "DB03591"}


def load_csv(path):
    with open(path, newline="", encoding="utf-8-sig") as f:
        return list(csv.DictReader(f))


def save_csv(path, rows, cols):
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def rank_key(x):
    return (x["worst_gscore"], x["mean_gscore"], x["spread"], x["drugbank_id"])


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def extract_archive_members(path):
    """Read archived stages by DrugBank ID, not by positional row count.

    Some Maestro subset exports lose s_lp_Variant and s_i_glide_gridfile.
    Their entries are ordered in PDB groups introduced by protein records.
    Thus infer PDB from each receptor title, falling back to gridfile when present.
    Reference compounds have no DrugBank ID; kept separate and NOT ranked.
    """
    stage = defaultdict(dict)
    references = Counter()
    receptor = None
    receptor_titles = []
    all_count = 0
    for ct_index, props in values_only_ct_headers(path):
        all_count += 1
        title = props.get("s_m_title", "")
        dbid = props.get("s_sd_DRUGBANK_ID", "").strip()
        gscore = props.get("r_i_glide_gscore", "")
        if not dbid and not gscore and any(p in title.upper() for p in PDBS):
            receptor = next(p for p in PDBS if p in title.upper())
            receptor_titles.append((ct_index, receptor, title))
            continue
        grid = props.get("s_i_glide_gridfile", "")
        pdb_grid = next((p for p in PDBS if p in grid.upper()), None)
        pdb = pdb_grid or receptor
        if dbid:
            if not pdb:
                raise ValueError(f"Candidate before receptor/group in {path}: CT {ct_index}")
            if not gscore:
                raise ValueError(f"Missing GlideScore for {dbid} in {path}")
            if pdb in stage[dbid]:
                raise ValueError(f"Duplicate candidate/PDB for {dbid}/{pdb} in {path}")
            stage[dbid][pdb] = {
                "gscore": float(gscore),
                "name": props.get("s_sd_GENERIC_NAME", ""),
                "variant": props.get("s_lp_Variant", ""),
                "ct_index": ct_index,
            }
        elif gscore:
            references[pdb] += 1
    if {x[1] for x in receptor_titles} != set(PDBS):
        raise ValueError(f"Expected one receptor group per PDB in {path}; found {receptor_titles}")
    for dbid, entries in stage.items():
        if set(entries) != set(PDBS):
            raise ValueError(f"{dbid} missing PDB(s) in {path}: {set(PDBS)-set(entries)}")
    return dict(stage), {"total_ct": all_count, "reference_pose_counts": dict(references),
                         "receptors": receptor_titles}


def combine_selector(selector, category_rows, audit_rows):
    if len(selector) != 549:
        raise ValueError(f"Expected 549 poses; got {len(selector)}")
    best = {}
    for r in selector:
        dbid, pdb = r["drugbank_id"].strip(), r["pdb"].strip()
        if pdb not in PDBS:
            raise ValueError(f"Unexpected PDB {pdb}")
        key = (dbid, pdb)
        score = float(r["glide_gscore"])
        if not math.isfinite(score):
            raise ValueError(f"Non-finite score for {key}")
        if key in best:
            raise ValueError(f"Duplicate receptor–compound pair {key}; input must be collapsed")
        best[key] = r
    ids = {dbid for dbid, _ in best}
    if len(ids) != 183 or len(best) != 549:
        raise ValueError(f"Wrong consensus dimensions: {len(ids)} compounds, {len(best)} pairs")
    if any((id_, p) not in best for id_ in ids for p in PDBS):
        raise ValueError("Missing receptor for one or more IDs")

    categories = {r["drugbank_id"].strip():r for r in category_rows}
    if len(categories) != len(category_rows) or set(categories) != ids:
        raise ValueError("Historical categories must include the same 183 distinct DrugBank IDs")
    distribution = Counter(r["category"] for r in category_rows)
    if distribution != CATS:
        raise ValueError(f"Historical categories altered: {distribution}")
    chemistry = {r["DrugBank ID"].strip():r for r in audit_rows} if audit_rows else {}
    if audit_rows and (len(chemistry) != 183 or set(chemistry) != ids):
        raise ValueError("Chemical audit must contain the same 183 DrugBank IDs")

    master = []
    for dbid in sorted(ids):
        ref = best[(dbid, PDBS[0])]
        scores = [float(best[(dbid,p)]["glide_gscore"]) for p in PDBS]
        chem = chemistry.get(dbid, {})
        r = {
            "drugbank_id":dbid, "compound_id":ref["compound_id"],
            "compound_name":ref["generic_name"],
            "historical_category":categories[dbid]["category"],
            "historical_category_reason":categories[dbid]["reason"],
            "revised_category":("C_reserva_ATP_like" if dbid in RU else categories[dbid]["category"]),
            "chemical_motifs_audited":chem.get("Motivos químicos detectados (SMILES)", ""),
            "physicochemical_audit_note":chem.get("Interpretación auditada (NUEVA)", ""),
            "chemical_audit_evidence_status":chem.get("Estatus de la explicación", ""),
            **{f"score_{p}":s for p,s in zip(PDBS,scores)},
            "worst_gscore":max(scores), "mean_gscore":sum(scores)/3,
            "spread":max(scores)-min(scores),
        }
        master.append(r)
    ordered = sorted(master, key=rank_key)
    for i, r in enumerate(ordered, 1): r["robustness_global_rank"] = i
    A = sorted((r for r in master if r["historical_category"]=="A_prioritario"),key=rank_key)
    B = sorted((r for r in master if r["historical_category"]=="B_revisar"),key=rank_key)
    initial30 = A + B[:5]
    if len(initial30)!=30: raise ValueError("Could not regenerate the original A+B Top 30")
    return ordered, initial30


def main():
    ap = argparse.ArgumentParser(description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--selector",required=True,type=Path,help="Historical 549 pose CSV")
    ap.add_argument("--categories",required=True,type=Path,help="183 historic labels A/B/C/D")
    ap.add_argument("--curation",type=Path,help="Optional private chemistry audit CSV")
    ap.add_argument("--historical30",type=Path,help="Optional first historical Top30 CSV for exact order verification")
    ap.add_argument("--revised30",required=True,type=Path,help="Archived NewTop30.maegz")
    ap.add_argument("--top10",required=True,type=Path,help="Archived NewTop10.maegz")
    ap.add_argument("--top3",required=True,type=Path,help="Archived final three ligands + reference structures .maegz")
    ap.add_argument("--references",required=True,type=Path,help="Publication references.json")
    ap.add_argument("--notes",required=True,type=Path,help="Historical decision_notes.json")
    ap.add_argument("--outdir",required=True,type=Path)
    args=ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    sources = [args.selector,args.categories,args.revised30,args.top10,args.top3,args.references,args.notes]
    if args.curation: sources.append(args.curation)
    if args.historical30: sources.append(args.historical30)
    master,initial30 = combine_selector(load_csv(args.selector),load_csv(args.categories),
                                       load_csv(args.curation) if args.curation else [])
    by_id={r["drugbank_id"]:r for r in master}
    history_order = [r["drugbank_id"] for r in initial30]
    if args.historical30:
        actual_order = [r["drugbank_id"] for r in load_csv(args.historical30)]
        if history_order != actual_order:
            raise AssertionError("Original A-first Top30 order differs from original worksheet")
    archived={}
    metainfo={}
    for label,path,n in (("revised30",args.revised30,30), ("top10",args.top10,10), ("top3",args.top3,3)):
        archive,summary=extract_archive_members(path)
        if len(archive)!=n or not set(archive).issubset(by_id):
            raise AssertionError(f"Archive {label} size or consensus inclusion mismatch")
        for dbid, entries in archive.items():
            for pdb,data in entries.items():
                expected=by_id[dbid][f"score_{pdb}"]
                if not math.isclose(data["gscore"],expected,rel_tol=0,abs_tol=1e-8):
                    raise AssertionError(f"Score mismatch: {label} {dbid} {pdb} {data['gscore']} != {expected}")
        archived[label]=set(archive)
        metainfo[label]=summary
    if not (archived["top3"] <= archived["top10"] <= archived["revised30"]):
        raise AssertionError("Not nested: Top 3 / Top 10 / revised Top 30")
    if archived["top3"] != set(FINAL):
        raise AssertionError(f"MM/GBSA candidate archive has unexpected names/IDs: {archived['top3']}")
    dropped = set(history_order) - archived["revised30"]
    added = archived["revised30"] - set(history_order)
    documented_dropped={"DB05552","DB05416","DB01136"}|RU
    documented_added={"DB06959","DB04760","DB06907","DB02936","DB01988","DB03471"}
    if dropped!=documented_dropped or added!=documented_added:
        raise AssertionError(f"Six substitutions disagree: drop {dropped}, add {added}")

    bibliography = json.loads(args.references.read_text(encoding="utf-8"))
    refs = {r["id"]: r for r in bibliography}
    notes=json.loads(args.notes.read_text(encoding="utf-8"))
    revised=archived["revised30"]; top10=archived["top10"]; top3=archived["top3"]
    out=[]
    for r in master:
        dbid = r["drugbank_id"]
        note=notes.get(dbid,{})
        refids=note.get("references",[])
        missing=[v for v in refids if v not in refs]
        if missing:raise AssertionError(f"Unknown bibliography IDs for {dbid}: {missing}")
        stage = ("top3_MMGBSA_input" if dbid in top3 else
                 "top10_archived" if dbid in top10 else
                 "revised30_archived" if dbid in revised else
                 "first30_only" if dbid in history_order else "consensus_only")
        if dbid in dropped: change = "removed_after_first_top30"
        elif dbid in added: change="added_to_revised_top30"
        else: change="unchanged_or_not_shortlisted"
        out.append({**r,
            "first30":int(dbid in history_order),
            "revised30":int(dbid in revised),
            "top10":int(dbid in top10),
            "top3_MMGBSA_input":int(dbid in top3),
            "latest_stage":stage,"revision_type":change,
            "documented_selection_reason":note.get("historical_reason",""),
            "retrospective_scientific_context":note.get("scientific_context",""),
            "evidence_level":note.get("evidence_level","archive membership only; individual manual reason not documented"),
            "references":";".join(refids),
            "citation_urls":";".join(refs[x]["url"] for x in refids),
            "method_context_references":"R01;R02;R03" + (";R08" if "fosf" in r["chemical_motifs_audited"].lower() else ""),
            "method_context_note":"General methodological support for repurposing/chemical curation, not an individual selection proof",
            "requires_confirming_original_judgement":note.get("needs_confirmation", "yes"),
        })
    columns=list(out[0])
    save_csv(args.outdir/"TRACE_183_to_MMGBSA_three.csv",out,columns)
    stages=[]
    for stage,items in (("183_consensus",set(by_id)),("30_first",set(history_order)),
                        ("30_revised",revised),("10",top10),("3_for_MMGBSA",top3)):
        stages.append({"stage":stage,"unique_compounds":len(items),"representative_ligand_poses":3*len(items)})
    save_csv(args.outdir/"stage_counts.csv",stages,list(stages[0]))
    decisions=[]
    for r in out:
        if r["revision_type"]!="unchanged_or_not_shortlisted" or r["top10"]:
            decisions.append({k:r[k] for k in ("drugbank_id","compound_name","first30","revised30","top10","top3_MMGBSA_input","revision_type","documented_selection_reason","retrospective_scientific_context","evidence_level","references","citation_urls","requires_confirming_original_judgement")})
    save_csv(args.outdir/"DECISIONS_and_candidates.csv", decisions, list(decisions[0]))
    selected3=sorted([r for r in out if r["top3_MMGBSA_input"]],key=lambda x:["DB01076","DB04879","DB00471"].index(x["drugbank_id"]))
    save_csv(args.outdir/"THREE_submitted_for_MMGBSA.csv",selected3,columns)
    audit={
        "status":"verified retrospective reconstruction; not prospective ranking or redocking",
        "counts":stages,
        "category_distribution":CATS,
        "first30_membership_exactly_matches_historical_csv": bool(args.historical30),
        "revised30_overlap_with_first30":len(set(history_order)&revised),
        "six_removed":sorted(dropped),"six_added":sorted(added),
        "top10_in_revised30":True,"top3_in_top10":True,
        "three_selected_for_MMGBSA":sorted(top3),
        "scores_match_original_selector_abs_tolerance":1e-8,
        "archive_metadata":metainfo,
        "input_sha256":{p.name:sha256(p) for p in sources},
        "caveats":[
          "Postscreening 30->10->3 is documented by archive membership. Exact historic numeric decision function was not found.",
          "Literature supports pharmacological rationale and assay relevance, NOT direct BUB1 inhibition.",
          "The source MAEGZ export does not fully preserve pose variants in all intermediate lists.",
          "The selection for MM/GBSA predates MM/GBSA results; never use those results to retrospectively justify the input list.",
          "Classification A/B/C/D was manually curated; retrospective chemistry audit does not recreate human judgments.",
          "PAINS absence was noted in historical email, not verified with an exhaustive computational filter.",
          "Early manual ATP-pocket pose checking cannot be independently reproduced from this script.",
        ]
    }
    (args.outdir/"VALIDATION.json").write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding="utf-8")
    print("OK: 183 -> 30 historical -> 30 revised (24 same, 6 replaced) -> 10 -> 3")
    print("Final 3 (DrugBank IDs):", ", ".join(f"{x} {FINAL[x]}" for x in FINAL))
    print("Wrote", args.outdir)


if __name__ == "__main__":
    main()
