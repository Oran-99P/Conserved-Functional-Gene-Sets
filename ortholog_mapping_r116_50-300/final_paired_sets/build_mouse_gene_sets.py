#!/usr/bin/env python3
"""Build mouse GO:MF gene sets paired with the approved human 50-300 collection (Ensembl Compara r116).

Mouse set of a GO group = union of ALL Compara r116 mouse orthologs (one2one, one2many, many2many) of its human
genes, de-duplicated by mouse Ensembl gene ID. No filtering, no substitution, no mouse size limit.
Sources (read-only):
  ../../git/results/GO_MF_gene_set_membership_50-300.tsv, GO_MF_gene_sets_50-300.{tsv,gmt}, run_metadata.json
  ../gene_ortholog_status_r116.tsv            per human GO symbol: Ensembl ID + ortholog status
  ../human_mouse_orthologs_all_r116.tsv       all validated human-mouse pairs (source of truth)
  ../GO_set_ortholog_coverage_r116.tsv        previous per-set statistics (cross-check only)
Usage: python3 build_mouse_gene_sets.py
"""
import csv, hashlib, json, re, statistics as st
from collections import Counter, defaultdict
from pathlib import Path

OUT = Path(__file__).resolve().parent; MAP = OUT.parent
RES = OUT.parent.parent / 'GO_MF_gene_sets' / 'results'
TAG = '50-300_r116'
rd = lambda p: list(csv.DictReader(open(p), delimiter='\t'))
memb = rd(RES / 'GO_MF_gene_set_membership_50-300.tsv')
hsets = rd(RES / 'GO_MF_gene_sets_50-300.tsv')
status = {r['human_gene_symbol_GO']: r for r in rd(MAP / 'gene_ortholog_status_r116.tsv')}
pairs = rd(MAP / 'human_mouse_orthologs_all_r116.tsv')
prev = {r['GO_ID']: r for r in rd(MAP / 'GO_set_ortholog_coverage_r116.tsv')}

orth = defaultdict(list); msym, mbio = {}, {}
for p in pairs:
    orth[p['human_ensembl_gene_id']].append(p)
    msym[p['mouse_ensembl_gene_id']] = p['mouse_gene_symbol']; mbio[p['mouse_ensembl_gene_id']] = p['mouse_gene_biotype']
term = {r['GO_ID']: r['GO_term'] for r in hsets}
human = defaultdict(list)
for r in memb: human[r['GO_ID']].append(r['gene_symbol'])

set_rows, mem_rows, pair_rows, gmt_lines = [], [], [], []
for go in [r['GO_ID'] for r in hsets]:                       # keep the order of the human TSV
    src = defaultdict(lambda: {'h_ids': set(), 'h_syms': set(), 'types': set()})
    c = Counter()
    for s in human[go]:
        st_ = status[s]; gid = st_['human_ensembl_gene_id']; cat = st_['ortholog_status']
        c['unresolved' if cat == 'unresolved_ensembl_id' else 'none' if cat == 'no_mouse_ortholog'
          else 'one2one' if cat == 'one2one' else 'non121'] += 1
        for p in orth.get(gid, []) if gid else []:
            d = src[p['mouse_ensembl_gene_id']]; d['h_ids'].add(gid); d['h_syms'].add(s); d['types'].add(p['orthology_type'])
    mice = sorted(src)
    for m in mice:
        d = src[m]
        mem_rows.append({'GO_ID': go, 'GO_term': term[go], 'mouse_ensembl_gene_id': m, 'mouse_gene_symbol': msym[m],
                         'mouse_gene_biotype': mbio[m], 'orthology_types': ';'.join(sorted(d['types'])),
                         'n_source_human_genes': len(d['h_ids']), 'source_human_ensembl_gene_ids': ';'.join(sorted(d['h_ids'])),
                         'source_human_gene_symbols_GO': ';'.join(sorted(d['h_syms']))})
    n = len(human[go]); mapped = c['one2one'] + c['non121']
    set_rows.append({'GO_ID': go, 'GO_term': term[go], 'n_mouse_genes': len(mice), 'n_human_genes_original': n,
                     'mouse_ensembl_gene_ids': ';'.join(mice), 'mouse_gene_symbols': ';'.join(msym[m] for m in mice)})
    pair_rows.append({'GO_ID': go, 'GO_term': term[go], 'n_human_genes': n, 'n_human_with_mouse_ortholog': mapped,
                      'pct_human_with_mouse_ortholog': round(100 * mapped / n, 2), 'n_distinct_mouse_genes': len(mice),
                      'n_human_one2one': c['one2one'], 'n_human_non_one2one_only': c['non121'],
                      'n_human_unmapped': c['none'] + c['unresolved'], 'n_human_no_mouse_ortholog': c['none'],
                      'n_human_unresolved_ensembl_id': c['unresolved']})
    gmt_lines.append('\t'.join([go, f'{term[go]} | mouse Ensembl gene IDs | Compara r116 orthologs of human GO:MF 50-300 set', *mice]))

def write(name, rows):
    with open(OUT / name, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter='\t', lineterminator='\n'); w.writeheader(); w.writerows(rows)
write(f'GO_MF_mouse_gene_sets_{TAG}.tsv', set_rows)
write(f'GO_MF_mouse_gene_set_membership_{TAG}.tsv', mem_rows)
write(f'GO_MF_human_mouse_paired_sets_{TAG}.tsv', pair_rows)
(OUT / f'GO_MF_mouse_gene_sets_{TAG}.gmt').write_text('\n'.join(gmt_lines) + '\n')

# ---------------- validation (re-reads the written files) ----------------
V = []
def check(ok, msg): V.append(('PASS ' if ok else 'FAIL ') + msg)
ms = rd(OUT / f'GO_MF_mouse_gene_sets_{TAG}.tsv'); mm = rd(OUT / f'GO_MF_mouse_gene_set_membership_{TAG}.tsv')
ps = rd(OUT / f'GO_MF_human_mouse_paired_sets_{TAG}.tsv')
gmt = {l.split('\t')[0]: l.rstrip('\n').split('\t') for l in open(OUT / f'GO_MF_mouse_gene_sets_{TAG}.gmt')}
hgmt = {l.split('\t')[0] for l in open(RES / 'GO_MF_gene_sets_50-300.gmt')}
H = set(term)
check(len(H) == 270 and len(ms) == len(ps) == len(gmt) == 270, '270 human groups and 270 mouse groups (TSV, GMT, paired summary)')
check({r['GO_ID'] for r in ms} == {r['GO_ID'] for r in ps} == set(gmt) == H == hgmt == set(human), 'GO IDs identical across human and mouse files')
meta = json.load(open(RES / 'run_metadata.json'))['outputs']
check(all(hashlib.sha256(open(RES / f, 'rb').read()).hexdigest() == h for f, h in meta.items()),
      'human collection unchanged (SHA-256 of all 5 human output files match run_metadata.json)')
pairset = {(p['human_ensembl_gene_id'], p['mouse_ensembl_gene_id']) for p in pairs}
check(all((h, r['mouse_ensembl_gene_id']) in pairset for r in mm for h in r['source_human_ensembl_gene_ids'].split(';')),
      'every mouse gene in every group is backed by a r116 pair with a human gene of that group')
gid_of = {s: status[s]['human_ensembl_gene_id'] for s in status}
expected = {(go, p['mouse_ensembl_gene_id']) for go in H for s in human[go] if gid_of[s] for p in orth.get(gid_of[s], [])}
check(expected == {(r['GO_ID'], r['mouse_ensembl_gene_id']) for r in mm},
      f'all {len(expected):,} expected group x mouse-gene memberships present (all one2many/many2many partners kept, nothing added)')
check(len({(r['GO_ID'], r['mouse_ensembl_gene_id']) for r in mm}) == len(mm)
      and all(len(v[2:]) == len(set(v[2:])) for v in gmt.values()), 'no duplicate mouse Ensembl IDs within any group')
cnt = Counter(r['GO_ID'] for r in mm)
check(all(int(r['n_mouse_genes']) == cnt[r['GO_ID']] == len(gmt[r['GO_ID']][2:]) == len(r['mouse_ensembl_gene_ids'].split(';') if r['mouse_ensembl_gene_ids'] else [])
          for r in ms) and all(int(r['n_distinct_mouse_genes']) == cnt[r['GO_ID']] for r in ps),
      'per-group mouse counts agree between set TSV, membership TSV, GMT and paired summary')
check(all(set(gmt[r['GO_ID']][2:]) == {x['mouse_ensembl_gene_id'] for x in mm if x['GO_ID'] == r['GO_ID']} for r in ms),
      'GMT gene lists equal membership TSV gene lists')
agree = all(int(r['n_human_genes']) == int(prev[r['GO_ID']]['n_human_genes'])
            and int(r['n_human_with_mouse_ortholog']) == int(prev[r['GO_ID']]['n_with_mouse_ortholog'])
            and int(r['n_human_one2one']) == int(prev[r['GO_ID']]['n_one2one'])
            and int(r['n_human_non_one2one_only']) == int(prev[r['GO_ID']]['n_non_one2one_only'])
            and int(r['n_human_no_mouse_ortholog']) == int(prev[r['GO_ID']]['n_no_mouse_ortholog'])
            and int(r['n_human_unresolved_ensembl_id']) == int(prev[r['GO_ID']]['n_unresolved_ensembl_id'])
            and int(r['n_distinct_mouse_genes']) == int(prev[r['GO_ID']]['n_distinct_mouse_orthologs_all']) for r in ps)
check(agree, 'paired summary agrees with GO_set_ortholog_coverage_r116.tsv for every group')
allm = {r['mouse_ensembl_gene_id'] for r in mm}
check(allm <= {p['mouse_ensembl_gene_id'] for p in pairs} and all(re.fullmatch(r'ENSMUSG\d{11}', m) for m in allm),
      'only Ensembl mouse IDs from the r116 mapping are used (no substitutions, no symbols as IDs)')
check(all(int(r['n_human_genes']) == int(r['n_human_with_mouse_ortholog']) + int(r['n_human_unmapped']) for r in ps),
      'per-group human counts are exhaustive (mapped + unmapped = original)')

# symbol issues among mouse genes used
nosym = sorted(m for m in allm if not msym[m])
bysym = defaultdict(set)
for m in allm:
    if msym[m]: bysym[msym[m]].add(m)
dupsym = {s: sorted(v) for s, v in bysym.items() if len(v) > 1}
sizes = sorted(int(r['n_mouse_genes']) for r in ms)
ratio = sorted(int(r['n_mouse_genes']) / int(r['n_human_genes_original']) for r in ms)
summ = {'groups': len(ms), 'unique_mouse_genes': len(allm), 'memberships': len(mm),
        'mouse_size': {'min': sizes[0], 'median': st.median(sizes), 'max': sizes[-1], 'mean': round(st.mean(sizes), 2)},
        'mouse_to_human_size_ratio': {'min': round(ratio[0], 3), 'median': round(st.median(ratio), 3), 'max': round(ratio[-1], 3)},
        'groups_zero_mouse_genes': [r['GO_ID'] for r in ms if r['n_mouse_genes'] == '0'],
        'groups_mouse_larger_than_human': sum(int(r['n_mouse_genes']) > int(r['n_human_genes_original']) for r in ms),
        'groups_mouse_over_300': [(r['GO_ID'], r['GO_term'], r['n_mouse_genes']) for r in ms if int(r['n_mouse_genes']) > 300],
        'groups_mouse_under_50': [(r['GO_ID'], r['GO_term'], r['n_mouse_genes']) for r in ms if int(r['n_mouse_genes']) < 50],
        'mouse_genes_without_symbol': nosym, 'mouse_symbols_shared_by_several_ids': dupsym, 'validation': V}
json.dump(summ, open(OUT / 'validation_summary.json', 'w'), indent=1)
print(json.dumps({k: v for k, v in summ.items() if k != 'mouse_genes_without_symbol'}, indent=1))
print('mouse genes without symbol:', len(nosym))
