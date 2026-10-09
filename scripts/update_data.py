#!/usr/bin/env python3
# Copyright (c) 2026 Florian LAUBSCHER
"""Build static Nucleye JSON files from existing reference-based alignments.

Python 3.9+, standard library only. Does not perform sequence alignment.
"""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import re
import sys
import tempfile

ROLES = ("F", "P", "R")
DNA = set("ACGTRYSWKMBDHVN")
COMPLEMENT = str.maketrans("ACGTRYSWKMBDHVN-", "TGCAYRSWMKVHDBN-")
SLUG = re.compile(r"[A-Za-z0-9_-]+\Z")


class DataError(ValueError):
    pass


def require(condition, message):
    if not condition:
        raise DataError(message)


def reverse_complement(sequence):
    return sequence.translate(COMPLEMENT)[::-1]


def read_json(path):
    def unique_pairs(pairs):
        result = {}
        for key, value in pairs:
            require(key not in result, f"{path}: clé JSON dupliquée {key!r}.")
            result[key] = value
        return result
    try:
        return json.loads(path.read_text(encoding="utf-8"), object_pairs_hook=unique_pairs)
    except json.JSONDecodeError as exc:
        raise DataError(f"{path}: JSON invalide ({exc}).") from exc


def read_fasta(path):
    records, current = {}, None
    for line_number, raw in enumerate(path.read_text(encoding="utf-8-sig").splitlines(), 1):
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            header = line[1:].split()
            require(header, f"{path}:{line_number}: identifiant vide.")
            current = header[0]
            require(current not in records, f"{path}: identifiant dupliqué {current}.")
            records[current] = []
        else:
            require(current is not None, f"{path}:{line_number}: séquence avant le premier en-tête.")
            sequence = "".join(line.split()).upper()
            require(set(sequence) <= DNA | {"-"}, f"{path}:{line_number}: alphabet ADN IUPAC ou '-' attendu.")
            records[current].append(sequence)
    records = {key: "".join(chunks) for key, chunks in records.items()}
    require(records and all(records.values()), f"{path}: FASTA vide ou séquence vide.")
    return records


def checked_alignment(path, reference_id, reference):
    alignment = read_fasta(path)
    require(reference_id in alignment, f"{path}: référence {reference_id!r} absente.")
    require(len({len(s) for s in alignment.values()}) == 1, f"{path}: les séquences alignées doivent avoir la même longueur.")
    require(alignment[reference_id].replace("-", "") == reference,
            f"{path}: la référence diffère de reference.fasta (hors gaps).")
    return alignment


def reference_map(aligned_reference):
    """Map alignment columns to 1-based ungapped reference coordinates."""
    position = 0
    coordinates = []
    for base in aligned_reference:
        if base != "-":
            position += 1
            coordinates.append(position)
        else:
            coordinates.append(None)
    return coordinates


def unverified_status(row):
    definition = (row.get("genbank") or {}).get("definition") or ""
    fallback = bool(re.match(r"^\s*UNVERIFIED\s*:", definition, re.I))
    flag = row.get("ncbi_unverified", fallback)
    require(type(flag) is bool, f"{row.get('id', '?')}: ncbi_unverified doit être un booléen.")
    return flag


def read_metadata(path, sample_ids):
    data = read_json(path)
    require(isinstance(data, dict) and data.get("schema_version") == 1
            and isinstance(data.get("records"), list), f"{path}: schema_version: 1 et records attendus.")
    metadata = {}
    for row in data["records"]:
        require(isinstance(row, dict) and isinstance(row.get("id"), str) and row["id"], f"{path}: identifiant absent.")
        key = row["id"]
        require(key not in metadata, f"{path}: identifiant dupliqué {key}.")
        require(row.get("collection_date") is None or isinstance(row["collection_date"], str), f"{path}: date invalide pour {key}.")
        loc = row.get("location") or {}
        require(isinstance(loc, dict), f"{path}: lieu invalide pour {key}.")
        for field in ("raw", "country_code", "country", "region", "locality"):
            require(loc.get(field) is None or isinstance(loc[field], str), f"{path}: {key}/{field} doit être un texte ou null.")
        require(loc.get("country_code") is None or re.fullmatch(r"[A-Za-z]{2}", loc["country_code"]), f"{path}: code pays invalide pour {key}.")
        metadata[key] = {"id": key, "collection_date": row.get("collection_date"), "location": loc,
                         "ncbi_unverified": unverified_status(row)}
    missing, extra = set(sample_ids) - metadata.keys(), metadata.keys() - set(sample_ids)
    require(not missing and not extra, f"{path}: identifiants différents du FASTA ; absents={sorted(missing)[:10]}, supplémentaires={sorted(extra)[:10]}.")
    return metadata


def read_set(folder, reference_id, reference):
    config_path = folder / "set.json"
    config = read_json(config_path)
    require(isinstance(config, dict) and config.get("schema_version") == 1, f"{config_path}: schema_version: 1 attendu.")
    require(config.get("id") == folder.name and SLUG.fullmatch(folder.name), f"{config_path}: id doit correspondre au nom du dossier.")
    require(isinstance(config.get("name"), str) and config["name"].strip(), f"{config_path}: nom absent.")
    require(isinstance(config.get("oligos"), dict) and set(config["oligos"]) == set(ROLES), f"{config_path}: rôles F, P et R attendus.")
    alignment_path = folder / "primers-aligned.fasta"
    alignment = checked_alignment(alignment_path, reference_id, reference)
    coordinates = reference_map(alignment[reference_id])
    result = {"id": config["id"], "name": config["name"], "oligos": {}}
    citation = config.get("citation")
    if citation is not None:
        require(isinstance(citation, dict) and all(isinstance(citation.get(key), str) and citation[key].strip()
                for key in ("label", "url")), f"{config_path}: citation invalide (label et url requis).")
        result["citation"] = citation
    used_ids = {reference_id}
    for role in ROLES:
        entry = config["oligos"][role]
        require(isinstance(entry, dict), f"{config_path}: entrée {role} invalide.")
        alignment_id, strand = entry.get("alignment_id"), entry.get("strand")
        require(isinstance(alignment_id, str) and alignment_id in alignment and alignment_id not in used_ids,
                f"{config_path}: alignment_id absent, dupliqué ou inconnu pour {role}.")
        require(strand in ("+", "-"), f"{config_path}: strand '+' ou '-' requis pour {role}.")
        used_ids.add(alignment_id)
        columns = [i for i, b in enumerate(alignment[alignment_id]) if b != "-"]
        require(columns, f"{alignment_path}: oligo {alignment_id} vide.")
        positions = [coordinates[i] for i in columns]
        require(None not in positions, f"{alignment_path}: {role} contient une base dans un gap de référence ; découpage fixe impossible.")
        require(positions == list(range(positions[0], positions[-1] + 1)),
                f"{alignment_path}: {role} ne couvre pas une plage de référence continue.")
        sequence = "".join(alignment[alignment_id][i] for i in columns)
        start, end = positions[0], positions[-1]
        if strand == "-":
            sequence, positions = reverse_complement(sequence), positions[::-1]
        result["oligos"][role] = {
            "name": entry.get("name") or role, "sequence": sequence, "strand": strand,
            "reference_start": start, "reference_end": end, "reference_positions": positions,
        }
        require(isinstance(result["oligos"][role]["name"], str), f"{config_path}: nom invalide pour {role}.")
    require(set(alignment) == used_ids, f"{alignment_path}: lignes supplémentaires non déclarées dans set.json.")
    return result


def extract_target(aligned_sequence, aligned_reference, oligo, coordinate_columns=None):
    if coordinate_columns is None:
        coordinate_columns = {pos: i for i, pos in enumerate(reference_map(aligned_reference)) if pos is not None}
    start, end = oligo["reference_start"], oligo["reference_end"]
    columns = [coordinate_columns[p] for p in range(start, end + 1)]
    # An insertion between covered reference bases cannot be shown faithfully
    # as a fixed-length target: mark that zone unavailable rather than drop it.
    inserted = [i for i in range(columns[0], columns[-1] + 1)
                if aligned_reference[i] == "-" and aligned_sequence[i] != "-"]
    if inserted:
        return None, {"reason": "insertion_in_target", "inserted_base_count": len(inserted)}
    target = "".join(aligned_sequence[i] for i in columns)
    if oligo["strand"] == "-":
        target = reverse_complement(target)
    return target, None


def sequence_bounds(sequence):
    """Observed extent in alignment columns, independent of any primer set."""
    return (len(sequence) - len(sequence.lstrip('-')), len(sequence.rstrip('-')) - 1)


def target_status(sequence, oligo, columns, bounds, unavailable=False):
    """Distinguish missing ends from internal alignment gaps; follow oligo 5'->3'."""
    first, last = bounds
    statuses = []
    for position in range(oligo['reference_start'], oligo['reference_end'] + 1):
        column = columns[position]
        base = sequence[column]
        status = ('uncovered' if column < first or column > last else
                  'unavailable' if unavailable else
                  'deletion' if base == '-' else
                  'base' if base in 'ACGT' else 'ambiguous')
        statuses.append(status)
    return statuses[::-1] if oligo['strand'] == '-' else statuses


def build_virus(folder):
    config_path = folder / "virus.json"
    config = read_json(config_path)
    require(isinstance(config, dict) and config.get("schema_version") == 1, f"{config_path}: schema_version: 1 attendu.")
    require(config.get("id") == folder.name and SLUG.fullmatch(folder.name), f"{config_path}: id doit correspondre au nom du dossier.")
    require(isinstance(config.get("name"), str) and config["name"].strip(), f"{config_path}: nom absent.")
    reference_id = config.get("reference_id")
    require(isinstance(reference_id, str) and reference_id, f"{config_path}: reference_id absent.")
    reference_path, alignment_path, metadata_path = (folder / name for name in ("reference.fasta", "alignment.fasta", "metadata.json"))
    references = read_fasta(reference_path)
    require(set(references) == {reference_id}, f"{reference_path}: exactement une référence, {reference_id!r}, attendue.")
    reference = references[reference_id]
    require("-" not in reference, f"{reference_path}: référence non alignée sans gaps attendue.")
    alignment = checked_alignment(alignment_path, reference_id, reference)
    samples = {key: seq for key, seq in alignment.items() if key != reference_id}
    require(samples, f"{alignment_path}: aucune séquence hors référence.")
    metadata = read_metadata(metadata_path, samples)
    folders = sorted(p.parent for p in folder.glob("*/set.json"))
    require(folders, f"{folder}: aucun dossier de set contenant set.json.")
    reference_info = {"id": reference_id, "length": len(reference)}
    catalog = {"id": config["id"], "name": config["name"], "reference": reference_info, "sets": []}
    catalog["last_updated"] = config.get("last_updated")
    columns = {pos: i for i, pos in enumerate(reference_map(alignment[reference_id])) if pos is not None}
    bounds = {key: sequence_bounds(seq) for key, seq in samples.items()}
    outputs = []
    for set_folder in folders:
        primer_set = read_set(set_folder, reference_id, reference)
        records = []
        for sample_id, sequence in samples.items():
            # Membership belongs to this set only. Never alter the source alignment.
            first, last = bounds[sample_id]
            if any(first > columns[o['reference_start']] or last < columns[o['reference_end']]
                   for o in primer_set['oligos'].values()):
                continue
            targets, notes, statuses = {}, {}, {}
            for role, oligo in primer_set["oligos"].items():
                targets[role], note = extract_target(sequence, alignment[reference_id], oligo, columns)
                statuses[role] = target_status(sequence, oligo, columns, bounds[sample_id], targets[role] is None)
                if note:
                    notes[role] = note
            # The full source annotations stay in metadata.json, shared by all sets.
            source = metadata[sample_id]
            location = source.get('location') or {}
            row = {"id": sample_id, "collection_date": source.get('collection_date'),
                   "ncbi_unverified": source['ncbi_unverified'],
                   "location": {key: location[key] for key in ('raw', 'country', 'country_code', 'region', 'locality') if key in location},
                   "targets": {primer_set["id"]: targets}, "target_status": statuses}
            if notes:
                row["target_notes"] = notes
            records.append(row)
        data = {
            "schema_version": 1, "virus_id": config["id"], "reference": reference_info,
            "coordinate_system": "1-based-inclusive", "set": primer_set,
            "records": records,
        }
        outputs.append((set_folder / "data.json", data))
        catalog["sets"].append({"id": primer_set["id"], "name": primer_set["name"],
                                "file": f"{folder.name}/{set_folder.name}/data.json"})
    return catalog, outputs


def atomic_write_json(path, value):
    encoded = (json.dumps(value, ensure_ascii=False, indent=2) + "\n").encode("utf-8")
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(dir=path.parent, prefix=".nucleye-", delete=False) as handle:
            temporary = Path(handle.name)
            handle.write(encoded)
        os.replace(temporary, path)
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)


def run(data_dir, check=False):
    folders = sorted(p.parent for p in data_dir.glob("*/virus.json"))
    require(folders, f"{data_dir}: aucun dossier de virus contenant virus.json.")
    viruses, outputs = [], []
    # Validate every input before replacing any generated file. Catalog last.
    for folder in folders:
        virus, virus_outputs = build_virus(folder)
        viruses.append(virus)
        outputs.extend(virus_outputs)
    outputs.append((data_dir / "viruses.json", {"schema_version": 1, "viruses": viruses}))
    if not check:
        for path, data in outputs:
            atomic_write_json(path, data)
    for path, data in outputs[:-1]:
        notes = sum(bool(row.get("target_notes")) for row in data["records"])
        print(f"{path.relative_to(data_dir)} : {len(data['records'])} séquences, {notes} avec une zone contenant une insertion.")
    print(f"{'Validation terminée' if check else 'JSON mis à jour'} : {len(viruses)} virus, {len(outputs)-1} sets.")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=Path(__file__).resolve().parents[1] / "data")
    parser.add_argument("--check", action="store_true", help="Valider les sources sans écrire les JSON.")
    args = parser.parse_args()
    try:
        run(args.data_dir.resolve(), args.check)
    except (DataError, OSError) as exc:
        print(f"Erreur : {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
