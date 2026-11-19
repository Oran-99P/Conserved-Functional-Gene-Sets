#!/usr/bin/env python3
"""Human-mouse ortholog mapping of the approved GO:MF 50-300 collection with Ensembl Compara release 116.

Inputs (read-only)
  GO:   ../git/results/GO_MF_gene_set_membership_50-300.tsv  (270 sets, 11,086 GO gene symbols)
  Ens.: ../ortholog_pilot_compara_20261118/downloads/release-116__*  (official Ensembl FTP files, downloaded 2026-11-18)
        - Compara.116.protein_default.homologies.tsv.gz from BOTH homologies/homo_sapiens/ and homologies/mus_musculus/
          (since r110 each human-mouse pair is stored in only one of the two per-genome files; in r116 all are in the
          mouse file), merged and de-duplicated by gene pair
        - Homo_sapiens.GRCh38.116.uniprot.tsv.gz, Homo_sapiens.GRCh38.116.gtf.gz, Mus_musculus.GRCm39.116.gtf.gz
Human ID resolution (same rules as ../ortholog_coverage_r116_20261118/)
  1. UniProt accession(s) -> Ensembl gene IDs (r116 xref), restricted to genes in the r116 GTF;
     several IDs -> keep the one whose GTF gene_name equals the GO symbol, if exactly one; otherwise ambiguous.
  2. No UniProt hit -> exact GO symbol == GTF gene_name, only if that name is unique in the GTF ("symbol_fallback").
  Otherwise the gene is unresolved (not_found / ambiguous). Symbols are annotations only; IDs are the keys.
No orthology filters are applied: every Compara r116 human-mouse ortholog pair is kept with its type and confidence.
Usage: python3 map_orthologs_r116.py   (writes next to this script)
"""
import csv, gzip, hashlib, json, re, statistics as st
from collections import Counter, defaultdict
from pathlib import Path

OUT = Path(__file__).resolve().parent
GO = Path('/home/user/oranpo/PHD/GO_gene_sets')
MEMB = GO / 'git/results/GO_MF_gene_set_membership_50-300.tsv'
SETS = GO / 'git/results/GO_MF_gene_sets_50-300.tsv'
DL = GO / 'ortholog_pilot_compara_20261118/downloads'
P = lambda s: DL / f'release-116__{s}'
F = {'compara_human': P('tsv__ensembl-compara__homologies__homo_sapiens__Compara.116.protein_default.homologies.tsv.gz'),
     'compara_mouse': P('tsv__ensembl-compara__homologies__mus_musculus__Compara.116.protein_default.homologies.tsv.gz'),
     'uniprot_xref': P('tsv__homo_sapiens__Homo_sapiens.GRCh38.116.uniprot.tsv.gz'),
     'human_gtf': P('gtf__homo_sapiens__Homo_sapiens.GRCh38.116.gtf.gz'),
     'mouse_gtf': P('gtf__mus_musculus__Mus_musculus.GRCm39.116.gtf.gz')}
TYPE = {'ortholog_one2one': 'one2one', 'ortholog_one2many': 'one2many', 'ortholog_many2many': 'many2many'}
LOW_121 = 50.0   # review flag only: sets with < 50% one-to-one coverage (never used to exclude anything)

# ---------------- inputs ----------------
memb = list(csv.DictReader(open(MEMB), delimiter='\t'))
set_name = {r['GO_ID']: r['GO_term'] for r in csv.DictReader(open(SETS), delimiter='\t')}
acc, sets_of = defaultdict(set), defaultdict(set)
for r in memb:
    acc[r['gene_symbol']] |= {a for a in r['uniprot_accessions'].split(';') if a}
    sets_of[r['gene_symbol']].add(r['GO_ID'])

def gtf(path):
    g = {}
    with gzip.open(path, 'rt') as fh:
        for l in fh:
            if '\tgene\t' not in l: continue
            a = dict(re.findall(r'(\S+) "([^"]*)"', l.split('\t')[8]))
            g[a['gene_id']] = (a.get('gene_name', ''), a.get('gene_biotype', ''))
    return g
hg, mg = gtf(F['human_gtf']), gtf(F['mouse_gtf'])
byname = defaultdict(list)
for gid, (n, _) in hg.items():
    if n: byname[n].append(gid)
uni = defaultdict(set)
with gzip.open(F['uniprot_xref'], 'rt') as fh:
    next(fh)
    for l in fh:
        c = l.split('\t'); uni[c[3].split('.')[0].split('-')[0]].add(c[0])

pairs, rows_seen, file_rows = {}, 0, Counter()
for key in ('compara_human', 'compara_mouse'):
    with gzip.open(F[key], 'rt') as fh:
        hdr = fh.readline().rstrip('\n').split('\t'); ix = {h: i for i, h in enumerate(hdr)}
        for l in fh:
            c = l.rstrip('\n').split('\t')
            s1, s2 = c[ix['species']], c[ix['homology_species']]
            if {s1, s2} != {'homo_sapiens', 'mus_musculus'}: continue
            file_rows[key] += 1
            hs = s1 == 'homo_sapiens'
            h, m = (c[ix['gene_stable_id']], c[ix['homology_gene_stable_id']]) if hs else (c[ix['homology_gene_stable_id']], c[ix['gene_stable_id']])
            hid, mid = (c[ix['identity']], c[ix['homology_identity']]) if hs else (c[ix['homology_identity']], c[ix['identity']])
            rec = {'type': TYPE[c[ix['homology_type']]], 'human_pct_identity': hid, 'mouse_pct_identity': mid,
                   'goc_score': c[ix['goc_score']], 'wga_coverage': c[ix['wga_coverage']],
                   'is_high_confidence': c[ix['is_high_confidence']], 'compara_homology_id': c[ix['homology_id']]}
            if (h, m) in pairs and pairs[(h, m)]['type'] != rec['type']:
                raise SystemExit(f'conflicting Compara types for {(h, m)}')
            pairs[(h, m)] = rec
orth = defaultdict(dict); mouse_partners = defaultdict(set)
for (h, m), rec in pairs.items():
    orth[h][m] = rec; mouse_partners[m].add(h)

# ---------------- per GO symbol ----------------
gene = {}
for sym in sorted(acc):
    accs = acc[sym]
    cand = set().union(*(uni.get(a, set()) for a in accs)) & hg.keys() if accs else set()
    if len(cand) > 1:
        t = {g for g in cand if hg[g][0] == sym}
        cand_all, cand = cand, (t if len(t) == 1 else cand)
    else:
        cand_all = cand
    if len(cand) == 1:
        gid = next(iter(cand)); how = 'uniprot' if len(cand_all) == 1 else 'uniprot_multi_symbol_tiebreak'
    elif not cand and len(byname.get(sym, [])) == 1:
        gid, how = byname[sym][0], 'symbol_fallback'
    else:
        gid = None
        how = (f'ambiguous_uniprot({len(cand)} genes)' if cand else
               f'ambiguous_symbol({len(byname[sym])} genes)' if byname.get(sym) else 'not_found')
    if gid is None: status = 'unresolved_ensembl_id'
    elif not orth.get(gid): status = 'no_mouse_ortholog'
    else:
        types = {r['type'] for r in orth[gid].values()}
        if len(types) != 1: raise SystemExit(f'mixed ortholog types for {gid}: {types}')
        status = types.pop()
    cands = sorted(cand_all or (set(byname.get(sym, [])) if gid is None else set()))
    gene[sym] = {'gid': gid, 'how': how, 'status': status, 'cands': cands}
syms_of_gid = defaultdict(set)
for s, g in gene.items():
    if g['gid']: syms_of_gid[g['gid']].add(s)

# ---------------- outputs ----------------
def write(name, rows):
    with open(OUT / name, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter='\t', lineterminator='\n'); w.writeheader(); w.writerows(rows)

status_rows = []
for s, g in gene.items():
    gid = g['gid']; o = orth.get(gid, {}) if gid else {}
    status_rows.append({
        'human_gene_symbol_GO': s, 'uniprot_accessions': ';'.join(sorted(acc[s])),
        'human_ensembl_gene_id': gid or '', 'human_ensembl_gene_name': hg[gid][0] if gid else '',
        'human_gene_biotype': hg[gid][1] if gid else '', 'id_resolution': g['how'],
        'candidate_ensembl_ids_if_unresolved': ';'.join(g['cands']) if not gid else '',
        'ensembl_id_shared_with_GO_symbols': ';'.join(sorted(syms_of_gid[gid] - {s})) if gid else '',
        'ortholog_status': g['status'], 'n_mouse_orthologs': len(o),
        'mouse_ensembl_gene_ids': ';'.join(sorted(o)), 'mouse_gene_symbols': ';'.join(mg.get(m, ('',))[0] for m in sorted(o)),
        'n_GO_sets': len(sets_of[s]), 'GO_IDs': ';'.join(sorted(sets_of[s]))})
write('gene_ortholog_status_r116.tsv', status_rows)

all_rows = []
for gid in sorted(syms_of_gid):
    for m, rec in sorted(orth.get(gid, {}).items()):
        all_rows.append({'human_ensembl_gene_id': gid, 'human_ensembl_gene_name': hg[gid][0],
                         'human_gene_symbols_GO': ';'.join(sorted(syms_of_gid[gid])),
                         'mouse_ensembl_gene_id': m, 'mouse_gene_symbol': mg.get(m, ('',))[0],
                         'mouse_gene_biotype': mg.get(m, ('', ''))[1], 'orthology_type': rec['type'],
                         'n_mouse_orthologs_of_human_gene': len(orth[gid]), 'n_human_orthologs_of_mouse_gene': len(mouse_partners[m]),
                         **{k: v for k, v in rec.items() if k != 'type'}, 'ensembl_release': 116})
write('human_mouse_orthologs_all_r116.tsv', all_rows)
one_rows = [r for r in all_rows if r['orthology_type'] == 'one2one']
write('human_mouse_orthologs_one2one_r116.tsv', one_rows)

# ---------------- per GO set ----------------
members = defaultdict(list)
for r in memb: members[r['GO_ID']].append(r['gene_symbol'])
set_rows = []
pct = lambda a, b: round(100 * a / b, 2)
for go in sorted(members):
    ss = members[go]; n = len(ss); c = Counter(gene[s]['status'] for s in ss)
    anyo = c['one2one'] + c['one2many'] + c['many2many']
    mice = {m for s in ss if gene[s]['gid'] for m in orth.get(gene[s]['gid'], {})}
    mice121 = {m for s in ss if gene[s]['status'] == 'one2one' for m in orth[gene[s]['gid']]}
    set_rows.append({'GO_ID': go, 'GO_term': set_name[go], 'n_human_genes': n,
                     'n_unique_human_ensembl_ids': len({gene[s]['gid'] for s in ss if gene[s]['gid']}),
                     'n_with_mouse_ortholog': anyo, 'pct_with_mouse_ortholog': pct(anyo, n),
                     'n_one2one': c['one2one'], 'pct_one2one': pct(c['one2one'], n),
                     'n_one2many': c['one2many'], 'n_many2many': c['many2many'],
                     'n_non_one2one_only': c['one2many'] + c['many2many'],
                     'n_no_mouse_ortholog': c['no_mouse_ortholog'], 'n_unresolved_ensembl_id': c['unresolved_ensembl_id'],
                     'n_distinct_mouse_orthologs_all': len(mice), 'n_distinct_mouse_orthologs_one2one': len(mice121),
                     'review_flag': f'low_one2one(<{LOW_121:g}%)' if pct(c['one2one'], n) < LOW_121 else ''})
write('GO_set_ortholog_coverage_r116.tsv', set_rows)

# ---------------- validation ----------------
V = []
def check(cond, msg):
    V.append(('PASS' if cond else 'FAIL') + ' ' + msg)
syms = {r['gene_symbol'] for r in memb}
check(len(gene) == len(syms) == 11086 and len(status_rows) == 11086, 'all 11,086 input symbols have one status row')
check(len(set_rows) == len(set_name) == 270 and set(members) == set(set_name), 'all 270 GO sets represented')
check(sum(r['n_human_genes'] for r in set_rows) == len(memb), f'per-set gene counts sum to {len(memb)} membership rows')
check(all(r['n_human_genes'] == r['n_with_mouse_ortholog'] + r['n_no_mouse_ortholog'] + r['n_unresolved_ensembl_id'] for r in set_rows)
      and all(r['n_with_mouse_ortholog'] == r['n_one2one'] + r['n_non_one2one_only'] for r in set_rows), 'per-set categories are exhaustive and exclusive')
check(all((r['human_ensembl_gene_id'], r['mouse_ensembl_gene_id']) in pairs for r in all_rows), 'every mapped pair exists in Compara r116')
check(len({(r['human_ensembl_gene_id'], r['mouse_ensembl_gene_id']) for r in all_rows}) == len(all_rows), 'no duplicate pairs in complete table')
c121h = Counter(r['human_ensembl_gene_id'] for r in one_rows); c121m = Counter(r['mouse_ensembl_gene_id'] for r in one_rows)
check(max(c121h.values()) == 1 and max(c121m.values()) == 1
      and all(r['n_mouse_orthologs_of_human_gene'] == 1 and r['n_human_orthologs_of_mouse_gene'] == 1 for r in one_rows),
      'one-to-one table: each human and each mouse gene appears once; neither has any other mouse/human ortholog in r116')
check(all(r['n_mouse_orthologs_of_human_gene'] > 1 or r['n_human_orthologs_of_mouse_gene'] > 1 for r in all_rows if r['orthology_type'] != 'one2one'),
      'every non-one-to-one pair really involves >1 partner on at least one side (no collapsing)')
expect_pairs = sum(len(orth.get(g, {})) for g in syms_of_gid)
check(len(all_rows) == expect_pairs, f'complete table keeps all {expect_pairs} Compara pairs of the resolved human genes')
HID, MID = re.compile(r'ENSG\d{11}$'), re.compile(r'ENSMUSG\d{11}$')
check(all(HID.match(r['human_ensembl_gene_id']) for r in all_rows) and all(MID.match(r['mouse_ensembl_gene_id']) for r in all_rows)
      and all(HID.match(g['gid']) for g in gene.values() if g['gid']) and not any(HID.match(s) or MID.match(s) for s in gene),
      'all IDs have Ensembl stable-ID format; no GO symbol is an Ensembl ID (the real symbol ENSA is a gene name)')
check(all(gene[s]['status'] != 'unresolved_ensembl_id' or not gene[s]['gid'] for s in gene)
      and all(gene[s]['status'] != 'no_mouse_ortholog' or (gene[s]['gid'] and not orth.get(gene[s]['gid'])) for s in gene),
      'unresolved genes have no ID; no-ortholog genes have an ID and no Compara mouse ortholog')
fb = [s for s in gene if gene[s]['how'] == 'symbol_fallback']
check(all(byname[s] == [gene[s]['gid']] for s in fb), f'all {len(fb)} symbol-fallback matches are unique in the r116 GTF')
check(file_rows['compara_human'] + file_rows['compara_mouse'] == len(pairs), 'Compara rows = unique pairs (no double counting between files)')

# sanity check vs previous 50-200 run
prev = {r['gene_symbol']: r for r in csv.DictReader(open(GO / 'ortholog_coverage_r116_20261118/gene_ortholog_status_r116.tsv'), delimiter='\t')}
pmap = {'one2one': 'one2one', 'one2many/many2many only': 'non121', 'no mouse ortholog': 'no_mouse_ortholog', 'unresolved': 'unresolved_ensembl_id'}
cur = lambda s: 'non121' if gene[s]['status'] in ('one2many', 'many2many') else gene[s]['status']
diff = [s for s in prev if s not in gene or cur(s) != pmap[prev[s]['category']] or (gene[s]['gid'] or '') != prev[s]['ensembl_gene_r116']]
check(set(prev) <= set(gene) and not diff, f'all {len(prev)} genes of the 50-200 run reproduce the same Ensembl ID and category ({len(diff)} differences)')

# ---------------- summary numbers ----------------
cs = Counter(g['status'] for g in gene.values())
by_gid = {}
for gid, ss in syms_of_gid.items(): by_gid[gid] = gene[next(iter(ss))]['status']
cg = Counter(by_gid.values())
res = Counter(g['how'].split('(')[0] for g in gene.values())
shared = {gid: sorted(ss) for gid, ss in syms_of_gid.items() if len(ss) > 1}
p121 = sorted(r['pct_one2one'] for r in set_rows); pany = sorted(r['pct_with_mouse_ortholog'] for r in set_rows)
bins = Counter('>=90' if x >= 90 else '75-90' if x >= 75 else '50-75' if x >= 50 else '<50' for x in p121)
flagged = sorted((r for r in set_rows if r['review_flag']), key=lambda r: r['pct_one2one'])
new_syms = sorted(set(gene) - set(prev)); cn = Counter(gene[s]['status'] for s in new_syms)
md5 = {k: hashlib.md5(v.read_bytes()).hexdigest() for k, v in F.items()}
manifest = {m['file']: m for m in json.load(open(DL / 'download_manifest.json')) if 'file' in m}
summary = {'gene_symbols': cs, 'ensembl_ids': cg, 'resolution': res, 'shared': shared, 'pairs_all': len(all_rows),
           'pairs_one2one': len(one_rows), 'distinct_mouse_all': len({r['mouse_ensembl_gene_id'] for r in all_rows}),
           'compara_rows': file_rows, 'validation': V}
json.dump({**summary, 'gene_symbols': dict(cs), 'ensembl_ids': dict(cg), 'resolution': dict(res), 'compara_rows': dict(file_rows),
           'input_md5': md5}, open(OUT / 'run_summary.json', 'w'), indent=1)

N = len(gene); NG = len(by_gid) + cs['unresolved_ensembl_id']
L = lambda k, d=N: f"{cs[k]:>6,}  {100 * cs[k] / d:6.2f}%"
anyo = cs['one2one'] + cs['one2many'] + cs['many2many']; anyg = cg['one2one'] + cg['one2many'] + cg['many2many']
rep = f"""HUMAN-MOUSE ORTHOLOG MAPPING OF THE GO:MF 50-300 COLLECTION (Ensembl Compara release 116)
Generated by map_orthologs_r116.py - 2026-11-18

INPUT
  {MEMB}
  270 GO:MF gene sets, {len(memb):,} memberships, {N:,} unique human GO gene symbols (unchanged).

SOURCE (Ensembl release 116, https://ftp.ensembl.org/pub/release-116/, downloaded {manifest[F['compara_mouse'].name]['downloaded_utc'][:10]})
  Compara.116.protein_default.homologies.tsv.gz  homologies/homo_sapiens/ ({file_rows['compara_human']} human-mouse rows)
                                                 homologies/mus_musculus/ ({file_rows['compara_mouse']:,} human-mouse rows)
  Both per-genome files are always read: since r110 each pair is written to only one of them (r116: all in the mouse file).
  Homo_sapiens.GRCh38.116.uniprot.tsv.gz, Homo_sapiens.GRCh38.116.gtf.gz, Mus_musculus.GRCm39.116.gtf.gz
  MD5 of the files used are in run_summary.json.

METHOD
  1. Human ID: UniProt accession -> Ensembl gene (r116 xref, genes in the r116 GTF); several hits -> the one whose
     GTF name equals the GO symbol; no UniProt hit -> exact, unique GTF gene_name ('symbol_fallback').
     Anything else is unresolved (not found, or ambiguous). Ensembl stable IDs are the keys; symbols are annotations.
  2. Orthology: all Compara r116 human-mouse ortholog pairs (one2one, one2many, many2many), no extra filters.
     Confidence columns (identity, GOC, WGA, is_high_confidence) are kept for information only.
  3. Status per gene: one2one / one2many / many2many / no_mouse_ortholog / unresolved_ensembl_id (exclusive).
     No mouse gene is chosen arbitrarily for one2many or many2many genes.

OVERALL RESULT                           by GO symbol (n={N:,})   by unique Ensembl gene (n={len(by_gid):,} + {cs['unresolved_ensembl_id']} unresolved)
  With >=1 mouse ortholog                {anyo:>6,}  {100*anyo/N:6.2f}%        {anyg:>6,}
    one-to-one                           {L('one2one')}        {cg['one2one']:>6,}
    one-to-many                          {L('one2many')}        {cg['one2many']:>6,}
    many-to-many                         {L('many2many')}        {cg['many2many']:>6,}
  No mouse ortholog in Compara r116      {L('no_mouse_ortholog')}        {cg['no_mouse_ortholog']:>6,}
  No reliable human Ensembl ID           {L('unresolved_ensembl_id')}        (not applicable)
  Human ID resolution: {res['uniprot']:,} UniProt, {res['uniprot_multi_symbol_tiebreak']} UniProt with symbol tie-break, {res['symbol_fallback']} symbol fallback,
                       {res['not_found']} not found, {res['ambiguous_uniprot'] + res['ambiguous_symbol']} ambiguous.
  Pairs: {len(all_rows):,} human-mouse ortholog pairs in total; {len(one_rows):,} one-to-one pairs;
         {summary['distinct_mouse_all']:,} distinct mouse genes.
  {len(shared)} Ensembl genes are shared by two GO symbols (counted once per Ensembl gene):
    {'; '.join(f"{g}: {'/'.join(s)}" for g, s in sorted(shared.items()))}

GO-SET COVERAGE (270 sets)
  % one-to-one per set:   median {st.median(p121):.1f}%, range {p121[0]:.1f}-{p121[-1]:.1f}%
                          >=90%: {bins['>=90']} sets | 75-90%: {bins['75-90']} | 50-75%: {bins['50-75']} | <50%: {bins['<50']}
  % any ortholog per set: median {st.median(pany):.1f}%, range {pany[0]:.1f}-{pany[-1]:.1f}%
  Flagged for review (one-to-one < {LOW_121:g}%; flag only, nothing excluded):
""" + ''.join(f"    {r['GO_ID']}  {r['GO_term'][:60]:<60} n={r['n_human_genes']:<4} 1:1 {r['pct_one2one']:5.1f}%  any {r['pct_with_mouse_ortholog']:5.1f}%\n" for r in flagged) + f"""
SANITY CHECK VS. THE 50-200 RUN (ortholog_coverage_r116_20261118)
  All 9,742 previous genes keep the same Ensembl ID and category. The {len(new_syms):,} genes new in 50-300:
  {cn['one2one']} one-to-one, {cn['one2many'] + cn['many2many']} non-one-to-one only, {cn['no_mouse_ortholog']} without ortholog, {cn['unresolved_ensembl_id']} unresolved.

VALIDATION
""" + ''.join(f'  {v}\n' for v in V) + f"""
LIMITATIONS
  - Compara assigns some conserved genes no ortholog or a questionable partner (the release pilot found e.g. CKS1B
    and SEM1 without a call in r116); these are reported as they are, not corrected.
  - Human genes not in Ensembl r116 or with ambiguous IDs (pseudogenes, X/Y PAR genes, ERVK loci, UniProt-only
    entries) are listed as unresolved, with candidate IDs where available.
  - Lineage-expanded families (immunoglobulin/TCR, HLA, olfactory receptors, KIR, PRAME, ...) have low one-to-one
    coverage; the affected sets are flagged above, not removed.
  - One GO symbol can share an Ensembl gene with another symbol (listed above); per-gene counts are given both ways.
  - Mapping is by Ensembl r116 IDs; datasets built on other Ensembl versions need ID matching (the release pilot
    found >99.4% of r116 one-to-one IDs present in Ensembl 93/98/110 annotations).

FILES
  human_mouse_orthologs_all_r116.tsv      all human-mouse ortholog pairs for the resolved genes (one row per pair)
  human_mouse_orthologs_one2one_r116.tsv  one-to-one pairs only
  gene_ortholog_status_r116.tsv           one row per GO gene symbol, incl. no-ortholog and unresolved genes
  GO_set_ortholog_coverage_r116.tsv       one row per GO set
  run_summary.json                        counts, validation results and input MD5s
  map_orthologs_r116.py                   this pipeline (reruns in under a minute)
"""
(OUT / 'summary_report.txt').write_text(rep)
print(rep)
