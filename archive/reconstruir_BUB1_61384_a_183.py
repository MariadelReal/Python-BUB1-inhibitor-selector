#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Auditoría y reconstrucción reproducible de BUB1 desde MAEGZ bruto hasta 183.

Método COMPUTACIONAL reconstruido (no el script histórico original):
    MAEGZ multiestructura -> excluir 2LAH/3ESL -> quitar copias duplicadas
    -> filtrar por GlideScore -> intersectar DrugBank ID en 4QPM/5DMZ/6F7B
    -> elegir la variante de GlideScore mínimo por DrugBank ID/PDB.

El archivo original contiene referencias de inhibidores sin DrugBank ID. Éstas
se auditan por separado; NO se integran como fármacos candidatos.

Limitación: La inspección manual del bolsillo ATP que antecedió los cortes
no se puede verificar geométricamente usando únicamente los metadatos MAE.
Este script usa la identidad del grid de docking correspondiente a cada PDB.

Dependencias: ÚNICAMENTE biblioteca estándar Python 3.9+ (csv, gzip, etc).

Ejemplo:
  python reconstruir_BUB1_61384_a_183.py \\
      --input Structures_DrugBankReferences_BUB1_raw.maegz \\
      --historical BUB1_SELECTOR_183x3_with_compound_ID.csv \\
      --outdir resultados_BUB1
"""
from __future__ import annotations

import argparse
import csv
import gzip
import json
import re
from collections import Counter, defaultdict
from pathlib import Path

PDBS = ("4QPM", "5DMZ", "6F7B")
EXCLUDED = ("2LAH", "3ESL")
COLUMNS = [
    "ct_index", "pdb", "record_type", "compound_id", "drugbank_id",
    "generic_name", "variant", "glide_gscore", "docking_score",
    "glide_emodel", "gridfile", "title", "source_file",
    "entry_id", "entry_name", "source_file_index", "drug_groups",
    "inchi_key", "smiles", "molecular_weight", "logp", "psa",
]


def values_only_ct_headers(path: Path):
    """Lee solamente la cabecera de cada f_m_ct, sin cargar átomos ni poses en RAM.

    En ESTE MAEGZ, las propiedades CT son claves de línea individuales antes
    de la primera línea ':::', seguidas de los mismos valores (1 por línea).
    Los bloques m_atom[] y m_bond[] no se interpretan para la criba por score.
    """
    names = []
    vals = []
    phase = 0
    ct_index = 0
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith("f_m_ct {"):
                if ct_index:
                    if len(vals) != len(names):
                        raise ValueError(f"Cabecera truncada en CT {ct_index}")
                    yield ct_index, dict(zip(names, vals))
                ct_index += 1
                names, vals, phase = [], [], 1
                continue
            if phase == 1:
                stripped = line.strip()
                if stripped == ":::":
                    phase = 2
                elif re.match(r"^[bisr]_\S+", stripped):
                    names.append(stripped.replace("\\_", "_"))
            elif phase == 2:
                vals.append(line.strip().strip('"'))
                if len(vals) == len(names):
                    phase = 3
        if ct_index:
            if len(vals) != len(names):
                raise ValueError(f"Cabecera final truncada en CT {ct_index}")
            yield ct_index, dict(zip(names, vals))


def pdb_from_props(p):
    combined = " ".join((p.get("s_i_glide_gridfile", ""),
                         p.get("s_m_Source_File", ""),
                         p.get("s_pdb_PDB_ID", ""))).upper()
    for pdb in PDBS + EXCLUDED:
        if pdb in combined:
            return pdb
    if "ORTOLOGO" in combined:
        return "3ESL"
    return "UNKNOWN"


def make_record(ct, p):
    pdb = pdb_from_props(p)
    source = p.get("s_m_Source_File", "")
    # Clasificar receptores ANTES de inspeccionar el nombre del archivo.
    # Los receptores también conservan un s_m_Source_File con «Drugbank_BUB1».
    if p.get("b_glide_receptor", "").strip() == "1":
        kind = "receptor"
    elif "REFERENCES_BUB1" in source.upper():
        kind = "reference"
    elif "DRUGBANK_BUB1" in source.upper():
        kind = "drugbank"
    else:
        kind = "other"
    return {
        "ct_index": ct,
        "pdb": pdb,
        "record_type": kind,
        "compound_id": p.get("s_lp_Variant", "").split("-")[0] or p.get("s_m_title", ""),
        "drugbank_id": p.get("s_sd_DRUGBANK_ID", ""),
        "generic_name": p.get("s_sd_GENERIC_NAME", ""),
        "variant": p.get("s_lp_Variant", ""),
        "glide_gscore": p.get("r_i_glide_gscore", ""),
        "docking_score": p.get("r_i_docking_score", ""),
        "glide_emodel": p.get("r_i_glide_emodel", ""),
        "gridfile": p.get("s_i_glide_gridfile", ""),
        "title": p.get("s_m_title", ""),
        "source_file": source,
        "entry_id": p.get("s_m_entry_id", ""),
        "entry_name": p.get("s_m_entry_name", ""),
        "source_file_index": p.get("i_m_Source_File_Index", ""),
        "drug_groups": p.get("s_sd_DRUG_GROUPS", ""),
        "inchi_key": p.get("s_sd_INCHI_KEY", ""),
        "smiles": p.get("s_sd_SMILES", ""),
        "molecular_weight": p.get("r_sd_MOLECULAR_WEIGHT", ""),
        "logp": p.get("r_sd_JCHEM_LOGP", ""),
        "psa": p.get("r_sd_JCHEM_POLAR_SURFACE_AREA", ""),
    }


def export_csv(path, records, fields):
    with path.open("w", newline="", encoding="utf-8-sig") as f:
        writer = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(records)


def retrieve_selected_structures(maegz: Path, dest: Path, indices):
    """Copia los CT de receptor y las poses elegidas (sin editar coordenadas)."""
    indices = set(indices)
    n = 0
    copied = 0
    with gzip.open(maegz, "rb") as source, gzip.open(dest, "wb", compresslevel=6) as sink:
        sink.write(b"{\n s_m_m2io_version\n :::\n 2.0.0\n}\n\n")
        include = False
        for line in source:
            if line.startswith(b"f_m_ct {"):
                n += 1
                include = n in indices
                if include:
                    copied += 1
            if include:
                sink.write(line)
        if copied != len(indices):
            raise AssertionError(f"CT exportados {copied}; solicitados {len(indices)}")
    return copied


def validate_historical(path, selected_rows):
    with path.open(encoding="utf-8-sig", newline="") as f:
        historical = list(csv.DictReader(f))
    expected = {(r["drugbank_id"], r["pdb"]): r for r in historical}
    observed = {(r["drugbank_id"], r["pdb"]): r for r in selected_rows}
    issues = []
    for key in sorted(set(expected) | set(observed)):
        a, b = expected.get(key), observed.get(key)
        if a is None or b is None:
            issues.append({"compound_pdb": key, "reason": "missing_composite_key"})
        elif (abs(float(a["glide_gscore"]) - float(b["glide_gscore"])) > 1e-8
              or a.get("variant", "") != b.get("variant", "")):
            issues.append({"compound_pdb": key, "reason": "score_or_variant_mismatch",
                           "history": [a["glide_gscore"], a.get("variant", "")],
                           "reconstructed": [b["glide_gscore"], b.get("variant", "")]})
    return {"historical_entries": len(historical),
            "matched_pairs": len(expected) - len(issues),
            "discrepancies": issues,
            "exact_match": len(issues) == 0 and len(expected) == 549}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--input", required=True, type=Path, help="MAEGZ bruto original")
    ap.add_argument("--historical", type=Path, help="CSV histórico de 549 entradas, opcional")
    ap.add_argument("--outdir", type=Path, default=Path("resultados_BUB1"))
    ap.add_argument("--cutoff", type=float, default=-7.4)
    ap.add_argument("--precutoff", type=float, default=-6.0)
    ap.add_argument("--no-maegz", action="store_true", help="No escribir MAEGZ reconstruido")
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    inventory = []
    n_all = 0
    n_excluded = Counter()
    dup = Counter()
    seen_drugbank = set()
    seen_reference = set()
    selected_receptors = {}
    drugbank = []
    reference = []
    for ct, p in values_only_ct_headers(args.input):
        n_all += 1
        r = make_record(ct, p)
        pdb, kind = r["pdb"], r["record_type"]
        if kind == "receptor":
            if pdb in PDBS and pdb not in selected_receptors:
                selected_receptors[pdb] = ct
            continue
        if pdb not in PDBS:
            n_excluded[pdb] += 1
            continue
        if kind == "drugbank":
            if not r["drugbank_id"] or not r["glide_gscore"]:
                continue
            # 5DMZ aparece repetido (el mismo pose/variante/score dos veces).
            # La llave incluye PDB y DBID para no fusionar compuestos distintos.
            key = (pdb, r["drugbank_id"], r["variant"], r["glide_gscore"])
            if key in seen_drugbank:
                dup[pdb] += 1
                continue
            seen_drugbank.add(key)
            drugbank.append(r)
        elif kind == "reference" and r["glide_gscore"]:
            # Los controles carecen de DBID; conservamos las variantes únicas
            # para el inventario, NO como candidatos.
            key = (pdb, r["variant"], r["glide_gscore"], r["title"])
            if key not in seen_reference:
                seen_reference.add(key)
                reference.append(r)

    pass_pre = [r for r in drugbank if float(r["glide_gscore"]) <= args.precutoff]
    pass_cut = [r for r in drugbank if float(r["glide_gscore"]) <= args.cutoff]
    references_pre = [r for r in reference if float(r["glide_gscore"]) <= args.precutoff]
    references_cut = [r for r in reference if float(r["glide_gscore"]) <= args.cutoff]
    coverage = defaultdict(set)
    for r in pass_cut:
        coverage[r["drugbank_id"]].add(r["pdb"])
    common_ids = {dbid for dbid, cov in coverage.items() if all(p in cov for p in PDBS)}

    # Mejor pose por PDB/DrugBank (en empates, conservar primera en MAEGZ).
    best = {}
    for r in pass_cut:
        if r["drugbank_id"] not in common_ids:
            continue
        key = (r["drugbank_id"], r["pdb"])
        if key not in best or float(r["glide_gscore"]) < float(best[key]["glide_gscore"]):
            best[key] = r

    master = []
    selected = []
    for dbid in common_ids:
        trio = [best[(dbid, p)] for p in PDBS]
        scores = [float(r["glide_gscore"]) for r in trio]
        first = trio[0]
        worst, mean, spread = max(scores), sum(scores) / 3, max(scores) - min(scores)
        master.append({"drugbank_id": dbid, "compound_id": first["compound_id"],
                       "generic_name": first["generic_name"], "drug_groups": first["drug_groups"],
                       "score_4QPM": scores[0], "score_5DMZ": scores[1],
                       "score_6F7B": scores[2], "worst_gscore": worst,
                       "mean_gscore": mean, "spread": spread})
        selected.extend(trio)
    master.sort(key=lambda r: (r["worst_gscore"], r["mean_gscore"],r["spread"], r["drugbank_id"]))
    rank = {r["drugbank_id"]: i+1 for i, r in enumerate(master)}
    for r in master: r["rank_robustez"] = rank[r["drugbank_id"]]
    selected.sort(key=lambda r: (rank[r["drugbank_id"]], PDBS.index(r["pdb"])))
    for r in selected: r["rank_robustez"] = rank[r["drugbank_id"]]
    mastercols = ["rank_robustez", "drugbank_id", "compound_id", "generic_name", "drug_groups",
                  "score_4QPM", "score_5DMZ", "score_6F7B", "worst_gscore", "mean_gscore", "spread"]
    export_csv(args.outdir / "183_compuestos_robustez.csv", master, mastercols)
    export_csv(args.outdir / "549_mejores_poses.csv", selected, ["rank_robustez"] + COLUMNS)
    export_csv(args.outdir / "referencias_separadas.csv", reference, COLUMNS)
    export_csv(args.outdir / "poses_drugbank_GlideScore_leq_m6.csv", pass_pre, COLUMNS)
    export_csv(args.outdir / "poses_drugbank_GlideScore_leq_m7p4.csv", pass_cut, COLUMNS)

    # Intersecciones Venn exclusivas (NO simples intersecciones inclusivas).
    pdb_sets = {p: {r["drugbank_id"] for r in pass_cut if r["pdb"] == p} for p in PDBS}
    union_ids = set.union(*pdb_sets.values())
    venn_counts = {}
    for p0 in (False, True):
        for p1 in (False, True):
            for p2 in (False, True):
                if not (p0 or p1 or p2):
                    continue
                members = (p0, p1, p2)
                name = "+".join(p for p, in_set in zip(PDBS, members) if in_set)
                venn_counts[name] = sum(all(((db in pdb_sets[p]) == in_set)
                                              for p, in_set in zip(PDBS, members))
                                        for db in union_ids)
    export_csv(args.outdir / "venn_3_PDB_regiones_exclusivas.csv",
               [{"region": key, "number_of_drugbank_ids": value}
                for key, value in sorted(venn_counts.items())],
               ["region", "number_of_drugbank_ids"])

    pdb_totals = Counter(r["pdb"] for r in drugbank)
    pdb_pre = Counter(r["pdb"] for r in pass_pre)
    pdb_cut = Counter(r["pdb"] for r in pass_cut)
    ref_pdb_pre = Counter(r["pdb"] for r in references_pre)
    ref_pdb_cut = Counter(r["pdb"] for r in references_cut)
    unique_cut = {p: len({r["drugbank_id"] for r in pass_cut if r["pdb"] == p}) for p in PDBS}
    audit = {
        "source_maegz": str(args.input),
        "total_ct_records_raw": n_all,
        "excluded_nonselected_pdb_entries": dict(n_excluded),
        "duplicated_drugbank_poses_removed": dict(dup),
        "drugbank_poses_in_3_pdbs_after_dedup": len(drugbank),
        "poses_per_pdb_before_score_cut": {p:pdb_totals[p] for p in PDBS},
        "drugbank_cutoff_pre": args.precutoff,
        "drugbank_poses_at_pre_cut": len(pass_pre),
        "reference_poses_at_pre_cut": len(references_pre),
        "combined_drugbank_plus_references_at_pre_cut": len(pass_pre)+len(references_pre),
        "poses_at_pre_cut_per_pdb": {p:pdb_pre[p] for p in PDBS},
        "ref_at_pre_cut_per_pdb": {p:ref_pdb_pre[p] for p in PDBS},
        "drugbank_cutoff_final": args.cutoff,
        "drugbank_poses_at_final_cut": len(pass_cut),
        "reference_poses_at_final_cut": len(references_cut),
        "combined_drugbank_plus_references_at_final_cut": len(pass_cut)+len(references_cut),
        "poses_at_final_cut_per_pdb": {p:pdb_cut[p] for p in PDBS},
        "reference_at_final_cut_per_pdb": {p:ref_pdb_cut[p] for p in PDBS},
        "unique_drugbank_ids_at_final_cut_per_pdb": unique_cut,
        "venn_exclusive_memberships": venn_counts,
        "common_drugbank_ids_3_of_3": len(common_ids),
        "best_drugbank_pdb_pairs": len(best),
        "receptor_ct_kept": selected_receptors,
        "historical_slide_numbers": {"initial": 61384, "pre_cavity": 22588,
                                     "pre_cutoff_minus6": 11691, "final_cutoff": 5049,
                                     "common_ids": 183, "precollapse_entries": 1063,
                                     "after_collapse_entries": 549},
        "notes": [
            "El original tiene dos copias idénticas de las poses de 5DMZ; se deduplican por PDB+DrugBank+variante+score.",
            "El conteo de 5,049 de la diapositiva no coincide exactamente con las 5,047 poses DrugBank que pasan <=-7.4; se mantienen referencias fuera del conjunto de candidatos.",
            "22,588 y 1,063 requieren definir exactamente qué entradas/controles se contaban tras selección manual; no se han forzado los totales.",
            "El script no reconstruye inspección geométrica humana de la cavidad ni reejecuta docking Glide.",
            "La robustez calculada aquí sólo ordena los 183: NO clasifica A/B/C/D ni produce directamente el Top 30 curado.",
        ],
    }
    if args.historical:
        audit["historical_comparison"] = validate_historical(args.historical, selected)

    if not args.no_maegz:
        selected_ct = [r["ct_index"] for r in selected]
        selected_ct += [selected_receptors[p] for p in PDBS if p in selected_receptors]
        maegz_out = args.outdir / "183_common_reconstruido_desde_raw.maegz"
        audit["ct_exported_to_maegz"] = retrieve_selected_structures(args.input, maegz_out, selected_ct)
        audit["output_maegz"] = str(maegz_out)
    audit_path = args.outdir / "auditoria_raw_a_183.json"
    audit_path.write_text(json.dumps(audit, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps({
        "raw": n_all, "drugbank_after_structure_and_dedup": len(drugbank),
        "pre_cut": len(pass_pre), "pre_cut_including_reference": len(pass_pre)+len(references_pre),
        "cut": len(pass_cut), "common": len(common_ids), "best_pairs": len(best),
        "history_exact_match":audit.get("historical_comparison",{}).get("exact_match",None),
        "outdir": str(args.outdir)
    }, indent=2, ensure_ascii=False))
    if len(common_ids) != 183 or len(best) != 549:
        raise AssertionError("La cifra de consenso no coincide con el conjunto histórico esperado.")
    if args.historical and not audit["historical_comparison"]["exact_match"]:
        raise AssertionError("No reproduce todas las variantes y scores históricos.")

if __name__ == "__main__":
    main()
