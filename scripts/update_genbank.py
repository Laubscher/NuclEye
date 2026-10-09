#!/usr/bin/env python3
# Copyright (c) 2026 Florian LAUBSCHER
"""Update one Nucleye virus from GenBank: alignment, metadata and blacklist.

Python 3.9+ (standard library), MAFFT 7 in PATH. Run manually, or from a scheduler.
"""
from __future__ import annotations

import argparse
from contextlib import contextmanager
from datetime import date, datetime, timezone
import fcntl
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile
import time
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode
from urllib.request import Request, urlopen
import uuid
import xml.etree.ElementTree as ET

from update_data import DNA, DataError, read_fasta, read_json, read_metadata, require, unverified_status

BASE_URL = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'
ACCESSION = re.compile(r'[A-Za-z][A-Za-z0-9_]*[0-9](?:\.[0-9]+)?\Z')
DEFAULTS = {
    'taxid': None, 'email': None, 'since': None, 'lookback_months': 12,
    'include_descendants': True, 'source': 'insdc',
    'min_length': 150, 'max_length': None, 'trim_ends': 25,
    'max_n_fraction': 0.02, 'max_n_run': 20, 'max_ambiguous_fraction': None,
    'exclude_unverified': False, 'exclude_lab_host': True,
    'excluded_divisions': ['PAT', 'SYN'],
    'excluded_terms': ['patent', 'synthetic construct', 'synthetic sequence', 'synthetic DNA',
        'synthetic RNA', 'artificial sequence', 'artificial construct', 'engineered construct',
        'cloning vector', 'expression vector', 'recombinant plasmid', 'infectious clone',
        'culture', 'cultured', 'cell culture', 'culture-derived', 'passaged', 'serial passage',
        'laboratory strain', 'lab strain', 'laboratory-adapted', 'cell-adapted', 'lab stock',
        'synthetic control', 'artificial control'],
    'mafft': 'mafft', 'alignment_mode': 'addfragments', 'threads': 2,
    'country_aliases': {},
}
MISSING = ('missing', 'not collected', 'not provided', 'restricted access', 'unknown', 'not applicable')
MONTHS = {name: i for i, name in enumerate(('jan','feb','mar','apr','may','jun','jul','aug','sep','oct','nov','dec'), 1)}
MUTABLE = {'alignment.fasta', 'metadata.json', 'blacklist.json', 'genbank_raw.fasta', 'virus.json'}


def now():
    return datetime.now(timezone.utc).isoformat(timespec='seconds')


def base_accession(value):
    return re.sub(r'\.[0-9]+$', '', value.strip().upper())


def encoded_json(value):
    return (json.dumps(value, ensure_ascii=False, indent=2) + '\n').encode('utf-8')


def digest(data):
    return hashlib.sha256(data).hexdigest() if data is not None else None


def atomic_bytes(path, data):
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(dir=path.parent, prefix='.nucleye-', delete=False) as handle:
            temporary = Path(handle.name)
            handle.write(data)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)


@contextmanager
def lock(folder):
    with (folder / '.genbank-update.lock').open('a') as handle:
        try:
            fcntl.flock(handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as exc:
            raise DataError('Une mise à jour est déjà en cours dans ce dossier.') from exc
        try:
            yield
        finally:
            fcntl.flock(handle, fcntl.LOCK_UN)


def archive_transaction(folder, pending, label):
    archive = folder / '.genbank-backups'
    archive.mkdir(exist_ok=True)
    pending.rename(archive / (datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S') + '-' + label + '-' + uuid.uuid4().hex[:8]))


def recover(folder):
    pending = folder / '.genbank-update-pending'
    if not pending.exists():
        return
    manifest_path = pending / 'transaction.json'
    if not manifest_path.exists():
        # Preparation had not reached the point where a target may be changed.
        shutil.rmtree(pending)
        return
    manifest = read_json(manifest_path)
    if manifest['phase'] == 'committed':
        archive_transaction(folder, pending, 'completed')
        return
    require(set(manifest['files']) <= MUTABLE, 'Journal de restauration invalide.')
    for name, item in manifest['files'].items():
        path = folder / name
        current = digest(path.read_bytes() if path.exists() else None)
        require(current in (item['before'], item['after']), f'{name} a été modifié depuis l’interruption ; restauration manuelle nécessaire.')
    for name, item in manifest['files'].items():
        original = pending / ('before-' + name)
        if item['before'] is None:
            (folder / name).unlink(missing_ok=True)
        else:
            content = original.read_bytes()
            require(digest(content) == item['before'], 'Sauvegarde de restauration altérée.')
            atomic_bytes(folder / name, content)
    archive_transaction(folder, pending, 'restored')
    print('Mise à jour interrompue restaurée avant de continuer.', file=sys.stderr)


def commit(folder, updates, report, expected):
    require(set(updates) <= MUTABLE, 'Fichier de sortie non autorisé.')
    # Refuse to overwrite concurrent user edits made during network/alignment work.
    for name, before in expected.items():
        path = folder / name
        require((path.read_bytes() if path.exists() else None) == before, f'{name} a changé pendant le traitement ; aucune mise à jour appliquée.')
    pending = folder / '.genbank-update-pending'
    pending.mkdir()
    files = {}
    for name, content in updates.items():
        before = expected[name]
        if before is not None:
            atomic_bytes(pending / ('before-' + name), before)
        atomic_bytes(pending / ('after-' + name), content)
        files[name] = {'before': digest(before), 'after': digest(content)}
    manifest = {'phase': 'prepared', 'files': files, 'run': report}
    atomic_bytes(pending / 'transaction.json', encoded_json(manifest))
    try:
        for name, content in updates.items():
            atomic_bytes(folder / name, content)
        manifest['phase'] = 'committed'
        atomic_bytes(pending / 'transaction.json', encoded_json(manifest))
    except BaseException:
        recover(folder)
        raise
    archive_transaction(folder, pending, 'completed')


def date_window(since, months, until=None):
    today = date.today()
    if since:
        require(re.fullmatch(r'\d{4}-\d{2}', since), '--since doit être AAAA-MM.')
        start = date.fromisoformat(since + '-01')
    else:
        index = today.year * 12 + today.month - 1 - months
        year, month = divmod(index, 12)
        start = date(year, month + 1, 1)
    if until:
        require(re.fullmatch(r'\d{4}-\d{2}-\d{2}', until), '--until doit être AAAA-MM-JJ.')
    end = date.fromisoformat(until) if until else today
    require(start <= end, 'Le début de la période est après sa fin.')
    return start, end


def collection_month(raw):
    if not raw:
        return None
    value = raw.strip()
    if re.fullmatch(r'\d{4}', value):
        return value if 1 <= int(value) <= 9999 else None
    if re.fullmatch(r'\d{4}-\d{2}', value):
        try:
            date.fromisoformat(value + '-01')
            return value
        except ValueError:
            return None
    if re.fullmatch(r'\d{4}-\d{2}-\d{2}(?:T.*)?', value):
        try:
            date.fromisoformat(value[:10])
            return value[:7]
        except ValueError:
            return None
    match = re.fullmatch(r'(?:(\d{1,2})-)?([A-Za-z]{3})-(\d{4})', value)
    if match and match[2].lower() in MONTHS:
        try:
            parsed = date(int(match[3]), MONTHS[match[2].lower()], int(match[1] or 1))
            return parsed.strftime('%Y-%m')
        except ValueError:
            return None
    return None  # ranges/uncertain dates are kept in the raw source field


class EntrezClient:
    def __init__(self, email, api_key=None):
        self.email, self.api_key, self.last_request = email, api_key, 0.0

    def xml(self, utility, **params):
        params.update(tool='nucleye_update', email=self.email)
        if self.api_key:
            params['api_key'] = self.api_key
        request = Request(BASE_URL + utility, data=urlencode(params).encode('utf-8'),
                          headers={'User-Agent': 'Nucleye/1.0', 'Content-Type': 'application/x-www-form-urlencoded'})
        # A conservative rate also works without an API key; no secret URLs in logs.
        for attempt in range(4):
            time.sleep(max(0, 0.4 - (time.monotonic() - self.last_request)))
            self.last_request = time.monotonic()
            try:
                with urlopen(request, timeout=60) as response:
                    payload = response.read()
                require(b'<!ENTITY' not in payload.upper(), 'Réponse XML avec entités non prise en charge.')
                root = ET.fromstring(payload)
                error = root.find('.//ERROR')
                require(error is None and root.tag != 'ERROR', 'NCBI a renvoyé une erreur ; aucune exclusion ajoutée.')
                return root
            except HTTPError as exc:
                if exc.code != 429 and not 500 <= exc.code < 600:
                    raise DataError(f'NCBI {utility}: HTTP {exc.code}.') from None
            except (URLError, TimeoutError, ET.ParseError):
                pass
            if attempt < 3:
                time.sleep(min(20, 2 ** (attempt + 1)))
        raise DataError(f'NCBI {utility}: requête échouée après quatre tentatives ; fichiers conservés.')

    def accessions(self, config, start, end):
        field = 'Organism:exp' if config['include_descendants'] else 'Organism:noexp'
        source = 'srcdb_ddbj/embl/genbank[PROP]' if config['source'] == 'insdc' else 'srcdb_genbank[PROP]'
        query = f'txid{config["taxid"]}[{field}] AND {source}'
        first = self.xml('esearch.fcgi', db='nuccore', term=query, usehistory='y', idtype='acc',
                         datetype='pdat', mindate=start.strftime('%Y/%m/%d'), maxdate=end.strftime('%Y/%m/%d'),
                         retmax=10000, retmode='xml')
        require(first.find('.//ErrorList') is None, 'NCBI ne reconnaît pas un champ de la recherche.')
        count = int(first.findtext('Count', '0'))
        result = [node.text for node in first.findall('IdList/Id')]
        key, webenv = first.findtext('QueryKey'), first.findtext('WebEnv')
        while len(result) < count:
            require(key and webenv, 'Recherche tronquée : historique NCBI absent.')
            page = self.xml('esearch.fcgi', db='nuccore', term='#' + key, WebEnv=webenv, usehistory='y',
                            idtype='acc', retstart=len(result), retmax=10000, retmode='xml')
            require(int(page.findtext('Count', '-1')) == count, 'Historique NCBI incohérent.')
            batch = [node.text for node in page.findall('IdList/Id')]
            require(batch, 'Recherche NCBI tronquée pendant la pagination.')
            result.extend(batch)
        require(len(result) == count and len(set(result)) == count, 'Résultats NCBI incomplets ou dupliqués.')
        require(all(isinstance(acc, str) and ACCESSION.fullmatch(acc) for acc in result), 'NCBI n’a pas renvoyé des accessions valides.')
        return result, query

    def records(self, accessions):
        for offset in range(0, len(accessions), 100):
            batch = accessions[offset:offset + 100]
            root = self.xml('efetch.fcgi', db='nuccore', id=','.join(batch), rettype='gb', retmode='xml')
            require(root.tag == 'GBSet', 'Réponse GenBank inattendue.')
            records = [parse_record(node) for node in root.findall('GBSeq')]
            covered = set()
            for record in records:
                aliases = {base_accession(a) for a in record['aliases']}
                matched = {a for a in batch if base_accession(a) in aliases}
                require(matched, 'GenBank a renvoyé une accession non demandée.')
                covered.update(matched)
            require(covered == set(batch), 'Lot GenBank incomplet ; aucune modification appliquée. Relancer la commande.')
            yield from records


def parse_record(node):
    accession = node.findtext('GBSeq_accession-version')
    require(accession is not None and ACCESSION.fullmatch(accession), 'Enregistrement GenBank sans accession.version valide.')
    qualifiers = {}
    for feature in node.findall('GBSeq_feature-table/GBFeature'):
        if feature.findtext('GBFeature_key') != 'source':
            continue
        for qualifier in feature.findall('GBFeature_quals/GBQualifier'):
            key, value = qualifier.findtext('GBQualifier_name'), qualifier.findtext('GBQualifier_value')
            if key and value and value not in qualifiers.setdefault(key, []):
                qualifiers[key].append(value)
    aliases = [accession, node.findtext('GBSeq_primary-accession') or base_accession(accession)]
    aliases += [n.text for n in node.findall('GBSeq_secondary-accessions/GBSecondary-accn') if n.text]
    return {'accession': accession.upper(), 'aliases': aliases,
            'sequence': ''.join((node.findtext('GBSeq_sequence') or '').split()).upper().replace('U', 'T'),
            'declared_length': int(node.findtext('GBSeq_length', '0')),
            'definition': node.findtext('GBSeq_definition'), 'organism': node.findtext('GBSeq_organism'),
            'division': node.findtext('GBSeq_division'),
            'create_date': node.findtext('GBSeq_create-date'), 'update_date': node.findtext('GBSeq_update-date'),
            'qualifiers': qualifiers}


def quality_reason(record, config):
    seq = record['sequence']
    if not seq:
        return 'sequence_absente'
    if set(seq) - DNA:
        return 'alphabet_non_iupac'
    input_length = record.get('original_length', len(seq))
    if input_length < config['min_length']:
        return 'sequence_trop_courte'
    if config['max_length'] is not None and input_length > config['max_length']:
        return 'sequence_trop_longue'
    # The N fraction is checked on the original sequence before trimming in update().
    if config['max_n_run'] is not None and any(len(run) > config['max_n_run'] for run in re.findall('N+', seq)):
        return 'serie_de_N_trop_longue'
    ambiguous = sum(base not in 'ACGT' for base in seq) / len(seq)
    if config['max_ambiguous_fraction'] is not None and ambiguous > config['max_ambiguous_fraction']:
        return 'trop_de_bases_ambigues'
    return None


def annotation_rejection(record, config):
    division = (record.get('division') or '').upper()
    if division in config['excluded_divisions']:
        return {'reason': 'division_exclue', 'field': 'division', 'term': division}
    definition = record.get('definition') or ''
    if config['exclude_unverified'] and re.match(r'^\s*UNVERIFIED\s*:', definition, re.I):
        return {'reason': 'ncbi_unverified', 'field': 'definition', 'term': 'UNVERIFIED'}
    qualifiers = record.get('qualifiers', {})
    if config['exclude_lab_host']:
        for value in qualifiers.get('lab_host', []):
            if value.strip() and not value.strip().lower().startswith(MISSING):
                return {'reason': 'hote_de_laboratoire', 'field': 'source.lab_host', 'term': value}
    fields = [('definition', definition)]
    for field in ('isolation_source', 'note', 'strain', 'lab_host', 'clone'):
        fields.extend((f'source.{field}', value) for value in qualifiers.get(field, []))
    for field, value in fields:
        for term in config['excluded_terms']:
            pattern = r'(?<!\w)' + r'[\s_-]+'.join(re.escape(p) for p in re.split(r'[\s_-]+', term.strip())) + r'(?!\w)'
            if re.search(pattern, value, re.I):
                return {'reason': 'annotation_exclue', 'field': field, 'term': term}
    return None


def record_metadata(record, config, known_countries):
    qualifiers = record['qualifiers']
    dates = qualifiers.get('collection_date', [])
    places = qualifiers.get('geo_loc_name') or qualifiers.get('country') or []
    countries = {p.split(':', 1)[0].strip() for p in places if not p.lower().startswith(MISSING)}
    country = next(iter(countries)) if len(countries) == 1 else None
    raw_place = ' | '.join(places) or None
    code = known_countries.get(country.casefold()) if country else None
    alias = config['country_aliases'].get(country, {})
    country, code = alias.get('country', country), alias.get('country_code', code)
    row = {'id': record['accession'], 'collection_date': collection_month(dates[0]) if len(dates) == 1 else None,
            'location': {'raw': raw_place, 'country': country, 'country_code': code},
            'genbank': {'accession': base_accession(record['accession']), 'accession_version': record['accession'],
                        'secondary_accessions': [a for a in record['aliases'] if a not in (record['accession'], base_accession(record['accession']))],
                        'organism': record['organism'], 'definition': record['definition'],
                        'create_date': record['create_date'], 'update_date': record['update_date'],
                        'collection_date_raw': dates, 'source_qualifiers': qualifiers, 'retrieved_at': now()}}
    row['ncbi_unverified'] = unverified_status(row)
    return row


def read_blacklist(path):
    data = read_json(path) if path.exists() else {'schema_version': 1, 'entries': []}
    require(isinstance(data, dict) and data.get('schema_version') == 1 and isinstance(data.get('entries'), list), 'blacklist.json invalide.')
    seen = set()
    for entry in data['entries']:
        acc = entry.get('accession') if isinstance(entry, dict) else None
        require(isinstance(acc, str) and ACCESSION.fullmatch(acc) and acc.upper() not in seen, 'Accession absente, invalide ou dupliquée dans blacklist.json.')
        require(isinstance(entry.get('reason'), str) and entry['reason'].strip(), f'{acc}: motif d’exclusion absent.')
        seen.add(acc.upper())
    return data


def blocked(accession, exclusions):
    return accession.upper() in exclusions or base_accession(accession) in exclusions


def local_state(folder):
    virus = read_json(folder / 'virus.json')
    reference_id = virus.get('reference_id')
    references = read_fasta(folder / 'reference.fasta')
    require(set(references) == {reference_id}, 'reference.fasta doit contenir uniquement reference_id.')
    require('-' not in references[reference_id], 'La référence source doit être sans gaps.')
    alignment = read_fasta(folder / 'alignment.fasta')
    require(reference_id in alignment and len({len(s) for s in alignment.values()}) == 1, 'Alignement local invalide.')
    require(alignment[reference_id].replace('-', '') == references[reference_id], 'La référence de l’alignement diffère de reference.fasta.')
    metadata = read_json(folder / 'metadata.json')
    read_metadata(folder / 'metadata.json', set(alignment) - {reference_id})
    known = {base_accession(key) for key in alignment}
    country_codes = {}
    for row in metadata['records']:
        gb = row.get('genbank') or {}
        for acc in [gb.get('accession'), gb.get('accession_version'), *gb.get('secondary_accessions', [])]:
            if acc:
                known.add(base_accession(acc))
        loc = row.get('location') or {}
        if loc.get('country') and loc.get('country_code'):
            country_codes[loc['country'].casefold()] = loc['country_code']
    blacklist = read_blacklist(folder / 'blacklist.json')
    return virus, alignment, metadata, known, country_codes, blacklist


def fasta_bytes(records, headers=None):
    headers = headers or {}
    return ''.join('>' + headers.get(key, key) + '\n' + '\n'.join(seq[i:i+80] for i in range(0, len(seq), 80)) + '\n'
                   for key, seq in records.items()).encode('utf-8')


def add_to_alignment(old, additions, reference_id, config):
    ordered = [reference_id] + [key for key in old if key != reference_id]
    old_names = {f'OLD{i:09d}': key for i, key in enumerate(ordered)}
    new_names = {f'NEW{i:09d}': key for i, key in enumerate(additions)}
    aliases = {**old_names, **new_names}
    with tempfile.TemporaryDirectory(prefix='nucleye-mafft-') as temporary:
        temp = Path(temporary)
        (temp / 'old.fasta').write_bytes(fasta_bytes({name: old[key] for name, key in old_names.items()}))
        (temp / 'new.fasta').write_bytes(fasta_bytes({name: additions[key] for name, key in new_names.items()}))
        command = [config['mafft'], '--nuc', '--auto', '--thread', str(config['threads']), '--adjustdirection',
                   '--' + config['alignment_mode'], str(temp / 'new.fasta'), str(temp / 'old.fasta')]
        with (temp / 'output.fasta').open('wb') as output, (temp / 'stderr.txt').open('wb') as error:
            result = subprocess.run(command, stdout=output, stderr=error, check=False)
        if result.returncode:
            tail = (temp / 'stderr.txt').read_text(errors='replace')[-3000:]
            raise DataError(f'MAFFT a échoué (code {result.returncode}) ; fichiers conservés.\n{tail}')
        output = read_fasta(temp / 'output.fasta')
    aligned, orientations = {}, {}
    for name, sequence in output.items():
        reversed_ = name.startswith('_R_')
        internal = name[3:] if reversed_ else name
        require(internal in aliases, 'Identifiant inattendu dans la sortie MAFFT.')
        key = aliases[internal]
        if internal in new_names:
            orientations[key] = '-' if reversed_ else '+'
        aligned[key] = sequence
    return aligned, orientations


def configuration(folder, args, virus):
    path = folder / 'genbank.json'
    supplied = read_json(path) if path.exists() else {}
    require(isinstance(supplied, dict), 'genbank.json doit contenir un objet.')
    require(not set(supplied) - set(DEFAULTS), 'Paramètre inconnu dans genbank.json : ' + ', '.join(set(supplied) - set(DEFAULTS)))
    config = {**DEFAULTS, **supplied}
    config['taxid'] = args.taxid or config['taxid'] or virus.get('taxid') or virus.get('orgid')
    config['email'] = args.email or config['email'] or os.environ.get('NCBI_EMAIL')
    for name in ('since', 'lookback_months', 'threads'):
        value = getattr(args, name)
        if value is not None:
            config[name] = value
    if args.lookback_months is not None and args.since is None:
        config['since'] = None
    require(type(config['taxid']) is int and config['taxid'] > 0, 'Renseigner le taxid NCBI avec --taxid ou dans genbank.json.')
    require(isinstance(config['email'], str) and re.fullmatch(r'[^\s@]+@[^\s@]+\.[^\s@]+', config['email']), 'Renseigner un e-mail NCBI avec --email ou NCBI_EMAIL.')
    for name in ('min_length', 'threads', 'lookback_months'):
        require(type(config[name]) is int and config[name] >= 1, f'{name} doit être un entier positif.')
    require(type(config['trim_ends']) is int and config['trim_ends'] >= 0, 'trim_ends doit être un entier positif ou nul.')
    require(config['max_n_run'] is None or type(config['max_n_run']) is int and config['max_n_run'] >= 0, 'max_n_run doit être un entier positif ou nul, ou null.')
    require(config['max_length'] is None or type(config['max_length']) is int and config['max_length'] >= config['min_length'], 'max_length invalide.')
    for name in ('max_ambiguous_fraction', 'max_n_fraction'):
        fraction = config[name]
        require(fraction is None or type(fraction) in (int, float) and 0 <= fraction <= 1, f'{name} doit être entre 0 et 1, ou null.')
    for name in ('exclude_unverified', 'exclude_lab_host'):
        require(type(config[name]) is bool, f'{name} doit être true ou false.')
    for name in ('excluded_divisions', 'excluded_terms'):
        require(isinstance(config[name], list) and all(isinstance(s, str) and s.strip() for s in config[name]), f'{name} doit être une liste de textes non vides.')
    config['excluded_divisions'] = [s.upper() for s in config['excluded_divisions']]
    require(type(config['include_descendants']) is bool, 'include_descendants doit être true ou false.')
    require(config['source'] in ('insdc', 'genbank'), 'source doit être insdc ou genbank.')
    require(config['alignment_mode'] in ('add', 'addfragments'), 'alignment_mode doit être add ou addfragments.')
    require(isinstance(config['mafft'], str) and config['mafft'], 'Chemin MAFFT invalide.')
    require(isinstance(config['country_aliases'], dict), 'country_aliases doit être un objet.')
    for source, alias in config['country_aliases'].items():
        require(isinstance(source, str) and isinstance(alias, dict) and isinstance(alias.get('country'), str), 'Alias de pays invalide.')
        require(alias.get('country_code') is None or re.fullmatch(r'[A-Z]{2}', alias['country_code']), 'Code pays invalide dans country_aliases.')
    return config


def update(folder, config, start, end, *, dry_run=False, client=None, aligner=add_to_alignment):
    expected = {name: (folder / name).read_bytes() if (folder / name).exists() else None
                for name in MUTABLE | {'reference.fasta', 'virus.json', 'genbank.json'}}
    virus, alignment, metadata, known, country_codes, blacklist = local_state(folder)
    require(virus.get('data_kind') != 'synthetic', 'Ce dossier contient une référence artificielle. Configurer un vrai dossier de virus avant une requête GenBank.')
    original_headers = {line[1:].split()[0]: line[1:].strip() for line in expected['alignment.fasta'].decode('utf-8-sig').splitlines() if line.startswith('>')}
    exclusions = {entry['accession'].upper() for entry in blacklist['entries']}
    client = client or EntrezClient(config['email'], os.environ.get('NCBI_API_KEY'))
    accessions, query = client.accessions(config, start, end)
    candidates = [acc for acc in accessions if base_accession(acc) not in known and not blocked(acc, exclusions)]
    print(f'{len(accessions)} résultats ; {len(candidates)} accessions nouvelles à récupérer.', flush=True)
    accepted, originals, new_metadata, rejected = {}, {}, [], []
    processed = set()
    report = {'started_at': now(), 'taxid': config['taxid'], 'date_field': 'PDAT', 'since': start.isoformat(),
              'until': end.isoformat(), 'query': query, 'search_count': len(accessions), 'candidate_count': len(candidates),
              'qc': {key: config[key] for key in ('min_length','max_length','trim_ends','max_n_fraction','max_n_run',
                       'max_ambiguous_fraction','exclude_unverified','exclude_lab_host','excluded_divisions','excluded_terms')},
              'alignment_mode': config['alignment_mode'], 'skipped_aliases': []}
    for record in client.records(candidates):
        acc = record['accession']
        aliases = {base_accession(a) for a in record['aliases']}
        if aliases & (known | processed) or any(blocked(a, exclusions) for a in record['aliases']):
            report['skipped_aliases'].append(acc)
            continue
        processed.update(aliases)
        original = record['sequence']
        rejection = annotation_rejection(record, config)
        if not rejection and (not original or set(original) - DNA):
            rejection = {'reason': 'sequence_absente' if not original else 'alphabet_non_iupac'}
        if not rejection and config['max_n_fraction'] is not None and original.count('N') / len(original) > config['max_n_fraction']:
            rejection = {'reason': 'trop_de_N'}
        trim = config['trim_ends']
        retained = original[trim:-trim] if trim else original
        processed_record = {**record, 'sequence': retained, 'original_length': len(original)}
        if not rejection:
            reason = quality_reason(processed_record, config) if retained else 'sequence_trop_courte_apres_decoupe'
            rejection = {'reason': reason} if reason else None
        if rejection:
            rejected.append({'accession': acc, **rejection, 'recorded_at': now()})
        else:
            accepted[acc] = retained
            originals[acc] = original
            metadata_row = record_metadata(record, config, country_codes)
            metadata_row['sequence_processing'] = {'trim_left': trim, 'trim_right': trim,
                'original_length': len(original), 'retained_length': len(retained)}
            new_metadata.append(metadata_row)
    report.update(accepted=list(accepted), rejected=rejected, finished_at=now(), dry_run=dry_run)
    if dry_run:
        print(json.dumps(report, ensure_ascii=False, indent=2))
        return report
    updates = {}
    if accepted:
        print(f'Ajout de {len(accepted)} séquences à l’alignement existant…', flush=True)
        combined, orientations = aligner(alignment, accepted, virus['reference_id'], config)
        for row in new_metadata:
            row['genbank']['alignment_orientation'] = orientations.get(row['id'])
        updates['alignment.fasta'] = fasta_bytes(combined, original_headers)
        raw_path = folder / 'genbank_raw.fasta'
        raw_sequences = read_fasta(raw_path) if raw_path.exists() else {}
        updates['genbank_raw.fasta'] = fasta_bytes({**raw_sequences, **originals})
    updated_metadata = {**metadata, 'records': [
        {**row, 'ncbi_unverified': unverified_status(row)}
        for row in [*metadata['records'], *new_metadata]]}
    if updated_metadata != metadata:
        updates['metadata.json'] = encoded_json(updated_metadata)
    entries = [{key: entry[key] for key in ('accession', 'reason', 'recorded_at', 'term') if key in entry}
               for entry in [*blacklist['entries'], *rejected]]
    if entries != blacklist['entries']:
        updates['blacklist.json'] = encoded_json({**blacklist, 'entries': entries})
    report['finished_at'] = now()
    updates['virus.json'] = encoded_json({**virus, 'last_updated': report['finished_at']})
    commit(folder, updates, report, expected)
    print(f'Mise à jour terminée : {len(accepted)} ajoutées, {len(rejected)} exclues. JSON des sets inchangés.')
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('virus_dir', type=Path, help='Dossier contenant virus.json, reference.fasta, alignment.fasta et metadata.json')
    parser.add_argument('--taxid', '--orgid', dest='taxid', type=int)
    parser.add_argument('--email')
    period = parser.add_mutually_exclusive_group()
    period.add_argument('--since', help='Premier mois inclus, AAAA-MM (date de mise à disposition publique)')
    parser.add_argument('--until', help='Dernier jour inclus, AAAA-MM-JJ ; aujourd’hui par défaut')
    period.add_argument('--lookback-months', type=int)
    parser.add_argument('--threads', type=int)
    parser.add_argument('--check', action='store_true', help='Valider uniquement les fichiers et la configuration, sans réseau')
    parser.add_argument('--dry-run', action='store_true', help='Rechercher, télécharger et contrôler, sans aligner ni modifier les données')
    args = parser.parse_args()
    try:
        folder = args.virus_dir.resolve()
        require(folder.is_dir(), 'Le dossier de virus n’existe pas.')
        with lock(folder):
            if args.dry_run or args.check:
                require(not (folder / '.genbank-update-pending').exists(), 'Mise à jour interrompue présente ; relancer sans --check/--dry-run pour restaurer.')
            else:
                recover(folder)
            virus = local_state(folder)[0]
            config = configuration(folder, args, virus)
            start, end = date_window(config['since'], config['lookback_months'], args.until)
            if args.check:
                print(f'Fichiers et configuration valides. Fenêtre PDAT : {start} → {end}.')
                return 0
            require(virus.get('data_kind') != 'synthetic', 'Le dossier joint est un exemple artificiel, pas une base GenBank. Renseigner le vrai virus et sa référence.')
            if not args.dry_run:
                require(shutil.which(config['mafft']) is not None, 'MAFFT introuvable ; installer MAFFT 7 ou renseigner son chemin dans genbank.json.')
            update(folder, config, start, end, dry_run=args.dry_run)
    except (DataError, OSError, ValueError, KeyError, TypeError) as exc:
        print(f'Erreur : {exc}', file=sys.stderr)
        return 1
    return 0


if __name__ == '__main__':
    sys.exit(main())
