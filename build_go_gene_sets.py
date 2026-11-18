#!/usr/bin/env python3
"""Conserved-Functional-Gene-Sets: human GO Molecular Function (GO:MF) gene sets.

 1. Ontology  : DB/go-basic.obo -- is_a edges, alt_id mapping, obsolete terms removed.
 2. Annotation: DB/goa_human.gaf.gz -- UniProtKB human MF annotations; NOT qualifier and
                ND evidence removed; all other evidence codes (including IEA) kept.
 3. Propagate : every gene annotated to a term is added to all of its is_a ancestors.
 4. Size      : keep terms with MIN..MAX genes (default 50-200), at any GO depth.
 5. Redundancy: terms are visited most-specific first (descending maximum is_a depth, so every
                descendant is visited before its ancestors). A term is removed when an already
                retained term has an identical gene set (exact duplicate), or when a retained
                descendant has Jaccard >= threshold (hierarchical near-duplicate; the more
                specific descendant is kept). All other overlaps are kept and documented.
                There is no coverage objective and no target number of groups.

Annotation-quality and gene-family columns are informational only; they never change membership.
"""
import argparse
import csv
import gzip
import hashlib
import json
import platform
import re
import statistics
import sys
from collections import Counter, defaultdict, deque
from datetime import datetime, timezone
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = 'GO:0003674'                     # molecular_function
NAMESPACE, ASPECT = 'molecular_function', 'F'
BINDING = 'GO:0005488'
UNINFORMATIVE = {ROOT, BINDING, 'GO:0005515'}   # root, binding, protein binding
EXPERIMENTAL = frozenset({'EXP', 'IDA', 'IPI', 'IMP', 'IGI', 'IEP', 'HTP', 'HDA', 'HMP', 'HGI', 'HEP'})
# Heuristic symbol patterns for lineage-expanded / poorly conserved families (informational flag only).
FAMILY_RE = re.compile(r'^(IG[HKL][VDJ]|TR[ABGD][VDJ]|HLA-|KIR\d|LILR|OR\d+[A-Z]+\d|TAS2R|KRTAP|ZNF\d|ZSCAN|'
                       r'PRAMEF|SPANX|USP17L|DEFB\d|NBPF|GAGE|MAGE[ABC]|CT45|SSX\d|TBC1D3|LCE\d|SIGLEC)')
REPORT_JACCARD, REPORT_OVERLAP = 0.5, 0.95   # overlap pairs written to the overlap table


# ----------------------------------------------------------------------------- input parsing
def read_obo(path):
    """Return (header, terms, alt_id -> primary ID). Only [Term] stanzas and is_a edges are used."""
    header, terms, current, section = {}, {}, None, None
    with open(path, encoding='utf-8') as fh:
        for raw in fh:
            line = raw.strip()
            if line.startswith('['):
                section, current = line, None
                continue
            key, _, value = line.partition(': ')
            if section is None and value:
                header.setdefault(key, value)
            if section != '[Term]' or not line or line.startswith('!'):
                continue
            if key == 'id':
                current = terms[value] = {'name': '', 'namespace': '', 'parents': [], 'alt_ids': [], 'obsolete': False}
            elif current is not None:
                if key == 'name':
                    current['name'] = value
                elif key == 'namespace':
                    current['namespace'] = value
                elif key == 'is_a':
                    current['parents'].append(value.split()[0])
                elif key == 'alt_id':
                    current['alt_ids'].append(value)
                elif key == 'is_obsolete':
                    current['obsolete'] = value == 'true'
    aliases = {alt: term for term, data in terms.items() for alt in data['alt_ids']}
    return header, terms, aliases


def read_gaf(path, aliases, terms):
    """Return ({(symbol, GO ID): evidence codes}, {symbol: UniProt accessions}, QC counts, header)."""
    pairs, accessions, counts, header = defaultdict(set), defaultdict(set), Counter(), {}
    with gzip.open(path, 'rt', encoding='utf-8') as fh:
        for line in fh:
            if line.startswith('!'):
                key, _, value = line[1:].strip().partition(': ')
                if value:
                    header.setdefault(key, value)   # first occurrence = this file's own header
                continue
            f = line.rstrip('\n').split('\t')
            if len(f) < 15:
                counts['malformed'] += 1
                continue
            db, acc, sym, qualifier, go, evidence, aspect, obj_type, taxon = (
                f[0], f[1], f[2], f[3], f[4], f[6], f[8], f[11], f[12])
            if aspect != ASPECT:
                continue
            counts['aspect_lines'] += 1
            if db != 'UniProtKB':
                counts['non_uniprotkb'] += 1
                continue
            if taxon.split('|')[0] != 'taxon:9606':          # first entry = taxon of the gene product
                counts['non_human'] += 1
                continue
            if 'NOT' in qualifier.split('|') or evidence == 'ND':
                counts['NOT_or_ND'] += 1
                continue
            if obj_type == 'protein_complex' or not sym or not go.startswith('GO:'):
                counts['invalid_object'] += 1
                continue
            go = aliases.get(go, go)
            if go not in terms or terms[go]['obsolete']:
                counts['missing_or_obsolete_term'] += 1
                continue
            pairs[sym, go].add(evidence)
            accessions[sym].add(acc)
            counts['kept_lines'] += 1
    return pairs, accessions, counts, header


# ----------------------------------------------------------------------------- ontology + propagation
def build_graph(terms):
    """is_a DAG of the MF namespace with minimum (BFS) and maximum (longest path) depth."""
    nodes = sorted(go for go, d in terms.items() if d['namespace'] == NAMESPACE and not d['obsolete'])
    parents = {go: terms[go]['parents'] for go in nodes}
    children = {go: [] for go in nodes}
    for go in nodes:
        for parent in parents[go]:
            if parent not in children:
                raise ValueError(f'Parent {parent} of {go} missing from the MF is_a graph')
            children[parent].append(go)
    min_depth, queue = {ROOT: 0}, deque([ROOT])
    while queue:
        node = queue.popleft()
        for child in children[node]:
            if child not in min_depth:
                min_depth[child] = min_depth[node] + 1
                queue.append(child)
    if len(min_depth) != len(nodes):
        raise ValueError(f'{len(nodes) - len(min_depth)} MF terms unreachable from {ROOT}')
    indegree = {go: len(parents[go]) for go in nodes}
    order, queue = [], deque([ROOT])
    while queue:
        node = queue.popleft()
        order.append(node)
        for child in children[node]:
            indegree[child] -= 1
            if indegree[child] == 0:
                queue.append(child)
    if len(order) != len(nodes):
        raise ValueError('MF is_a graph contains a cycle')
    max_depth = {ROOT: 0}
    for node in order:
        for child in children[node]:
            max_depth[child] = max(max_depth.get(child, 0), max_depth[node] + 1)
    return {'nodes': nodes, 'parents': parents, 'children': children, 'order': order,
            'min_depth': min_depth, 'max_depth': max_depth}


def propagate(graph, pairs, gene_index):
    """Bitmask per term = genes annotated to the term or to any is_a descendant."""
    direct = defaultdict(int)
    for symbol, go in pairs:
        direct[go] |= 1 << gene_index[symbol]
    masks = {}
    for node in reversed(graph['order']):
        mask = direct[node]
        for child in graph['children'][node]:
            mask |= masks[child]
        masks[node] = mask
    return masks, direct


def closure(start, edges):
    """All nodes reachable from start via edges (start excluded)."""
    seen, stack = set(), list(edges[start])
    while stack:
        node = stack.pop()
        if node not in seen:
            seen.add(node)
            stack.extend(edges[node])
    return seen


def members(mask, genes):
    out = []
    while mask:
        bit = mask & -mask
        out.append(genes[bit.bit_length() - 1])
        mask ^= bit
    return out


# ----------------------------------------------------------------------------- selection
def reduce_redundancy(candidates, masks, graph, ancestors, jaccard_threshold):
    """Most-specific-first visit; remove exact duplicates and ancestors of near-duplicate retained descendants."""
    order = sorted(candidates, key=lambda t: (-graph['max_depth'][t], masks[t].bit_count(), t))
    kept, decisions = [], {}
    for term in order:
        decision = None
        for other in kept:
            if masks[other] == masks[term]:
                decision = ('removed_exact_duplicate', other)
                break
        if decision is None:
            best = None
            for other in kept:
                if term in ancestors[other]:                       # other is a descendant of term
                    j = jaccard(masks[term], masks[other])
                    if j >= jaccard_threshold and (best is None or j > best[0]):
                        best = (j, other)
            if best:
                decision = ('removed_hierarchical_near_duplicate', best[1])
        if decision is None:
            kept.append(term)
            decision = ('retained', '')
        decisions[term] = decision
    return sorted(kept), decisions


def jaccard(a, b):
    shared = (a & b).bit_count()
    return shared / (a.bit_count() + b.bit_count() - shared) if shared else 0.0


def relation(a, b, ancestors):
    if a in ancestors[b]:
        return 'A_is_ancestor_of_B'
    if b in ancestors[a]:
        return 'A_is_descendant_of_B'
    return 'not_hierarchical'


def size_summary(sizes):
    q1, median, q3 = statistics.quantiles(sizes, n=4, method='inclusive')
    return {'min': min(sizes), 'q1': q1, 'median': median, 'mean': round(statistics.mean(sizes), 2),
            'q3': q3, 'max': max(sizes)}


# ----------------------------------------------------------------------------- output helpers
def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as fh:
        for block in iter(lambda: fh.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def write_tsv(path, rows, fields):
    with open(path, 'w', newline='', encoding='utf-8') as fh:
        writer = csv.DictWriter(fh, fieldnames=fields, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


def pct(n, d):
    return round(100 * n / d, 2) if d else 0.0


# ----------------------------------------------------------------------------- main
def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--db', type=Path, default=HERE / 'DB', help='folder with go-basic.obo and goa_human.gaf.gz')
    p.add_argument('--outdir', type=Path, default=HERE / 'results')
    p.add_argument('--min-genes', type=int, default=50)
    p.add_argument('--max-genes', type=int, default=200)
    p.add_argument('--jaccard', type=float, default=0.8, help='hierarchical near-duplicate threshold')
    a = p.parse_args()
    if not 1 <= a.min_genes <= a.max_genes:
        p.error('invalid gene-count range')
    if not 0 < a.jaccard <= 1:
        p.error('--jaccard must be in (0, 1]')
    obo, gaf = a.db / 'go-basic.obo', a.db / 'goa_human.gaf.gz'
    for f in (obo, gaf):
        if not f.is_file():
            p.error(f'missing input: {f}')
    a.outdir.mkdir(parents=True, exist_ok=True)
    tag = f'{a.min_genes}-{a.max_genes}'

    # 1-3: parse and propagate
    obo_header, terms, aliases = read_obo(obo)
    pair_evidence, accessions, gaf_qc, gaf_header = read_gaf(gaf, aliases, terms)
    graph = build_graph(terms)
    nodes = set(graph['nodes'])
    pairs = sorted(pr for pr in pair_evidence if pr[1] in nodes)
    genes = sorted({sym for sym, _ in pairs})
    gene_index = {sym: i for i, sym in enumerate(genes)}
    masks, direct = propagate(graph, pairs, gene_index)
    if masks[ROOT].bit_count() != len(genes):
        raise AssertionError('root does not contain every annotated gene')
    names = {go: terms[go]['name'] for go in nodes}
    terms_by_gene = defaultdict(set)
    for sym, go in pairs:
        terms_by_gene[sym].add(go)
    informative_universe = sum(1 for g in genes if not terms_by_gene[g] <= UNINFORMATIVE)

    # 4-5: size filter and redundancy reduction
    candidates = sorted(t for t in graph['nodes'] if a.min_genes <= masks[t].bit_count() <= a.max_genes)
    ancestors = {t: closure(t, graph['parents']) for t in candidates}
    final, decisions = reduce_redundancy(candidates, masks, graph, ancestors, a.jaccard)
    final_set = set(final)

    # candidate audit
    merged_into = defaultdict(list)
    cand_rows = []
    for t in candidates:
        status, kept_by = decisions[t]
        row = {'GO_ID': t, 'GO_term': names[t], 'min_depth': graph['min_depth'][t], 'max_depth': graph['max_depth'][t],
               'n_genes': masks[t].bit_count(), 'status': status, 'retained_GO_ID': kept_by,
               'retained_GO_term': names.get(kept_by, ''), 'jaccard': '', 'overlap_coefficient': '', 'relation': ''}
        if kept_by:
            merged_into[kept_by].append(t)
            shared = (masks[t] & masks[kept_by]).bit_count()
            row.update(jaccard=round(jaccard(masks[t], masks[kept_by]), 4),
                       overlap_coefficient=round(shared / min(masks[t].bit_count(), masks[kept_by].bit_count()), 4),
                       relation={'A_is_ancestor_of_B': 'removed_is_ancestor_of_retained',
                                 'A_is_descendant_of_B': 'removed_is_descendant_of_retained',
                                 'not_hierarchical': 'not_hierarchical'}[relation(t, kept_by, ancestors)])
        cand_rows.append(row)

    # pairwise overlap among final sets
    overlap_rows, stats = [], Counter()
    best_partner = {t: (0.0, '') for t in final}
    contained_in, contains = Counter(), Counter()
    for i, x in enumerate(final):
        for y in final[i + 1:]:
            shared = (masks[x] & masks[y]).bit_count()
            if not shared:
                continue
            nx, ny = masks[x].bit_count(), masks[y].bit_count()
            j, oc = shared / (nx + ny - shared), shared / min(nx, ny)
            rel = relation(x, y, ancestors)
            stats['pairs_sharing_genes'] += 1
            stats['ancestor_descendant_pairs'] += rel != 'not_hierarchical'
            stats['pairs_jaccard_ge_0.5'] += j >= 0.5
            stats['pairs_jaccard_ge_threshold'] += j >= a.jaccard
            stats['hierarchical_pairs_jaccard_ge_threshold'] += j >= a.jaccard and rel != 'not_hierarchical'
            stats['exact_duplicate_pairs'] += masks[x] == masks[y]
            stats['pairs_overlap_coefficient_ge_0.95'] += oc >= 0.95
            stats['fully_nested_pairs'] += shared == min(nx, ny)
            if shared == min(nx, ny):
                small, big = (x, y) if nx <= ny else (y, x)
                contained_in[small] += 1
                contains[big] += 1
            for s, o in ((x, y), (y, x)):           # highest Jaccard; ties -> lowest GO ID
                best_j, best_o = best_partner[s]
                if j > best_j or (j == best_j and o < best_o):
                    best_partner[s] = (j, o)
            if j >= REPORT_JACCARD or oc >= REPORT_OVERLAP:
                overlap_rows.append({'GO_ID_A': x, 'GO_term_A': names[x], 'n_genes_A': nx, 'GO_ID_B': y,
                                     'GO_term_B': names[y], 'n_genes_B': ny, 'shared_genes': shared,
                                     'jaccard': round(j, 4), 'overlap_coefficient': round(oc, 4), 'relation': rel})

    # per-set table, membership table, GMT
    set_rows, member_rows = [], []
    top_level = {t: sorted(names[x] for x in ancestors[t] | {t} if graph['min_depth'][x] == 1) for t in final}
    for t in final:
        desc = closure(t, graph['children']) | {t}
        gene_list = members(masks[t], genes)
        support = {g: set().union(*(pair_evidence[g, go] for go in terms_by_gene[g] & desc)) for g in gene_list}
        n = len(gene_list)
        accession_only = [g for g in gene_list if g in accessions[g]]
        for g in gene_list:
            member_rows.append({'GO_ID': t, 'gene_symbol': g, 'uniprot_accessions': ';'.join(sorted(accessions[g])),
                                'annotation': 'direct' if direct[t] >> gene_index[g] & 1 else 'propagated',
                                'evidence_codes': ';'.join(sorted(support[g]))})
        set_rows.append({
            'GO_ID': t, 'GO_term': names[t], 'min_depth': graph['min_depth'][t], 'max_depth': graph['max_depth'][t],
            'n_genes': n, 'top_level_MF': ';'.join(top_level[t]), 'binding_descendant': BINDING in ancestors[t],
            'merged_redundant_terms': ';'.join(sorted(merged_into[t])),
            'max_jaccard_other_set': round(best_partner[t][0], 4), 'most_similar_set': best_partner[t][1],
            'n_sets_containing_this': contained_in[t], 'n_sets_contained_in_this': contains[t],
            'pct_IEA_only': pct(sum(s == {'IEA'} for s in support.values()), n),
            'pct_no_experimental': pct(sum(not s & EXPERIMENTAL for s in support.values()), n),
            'pct_IPI_only': pct(sum(s == {'IPI'} for s in support.values()), n),
            'pct_family_heuristic': pct(sum(bool(FAMILY_RE.match(g)) for g in gene_list), n),
            'n_accession_only_ids': len(accession_only),
            'genes': ';'.join(gene_list)})

    out = {
        'gmt': a.outdir / f'GO_MF_gene_sets_{tag}.gmt',
        'sets': a.outdir / f'GO_MF_gene_sets_{tag}.tsv',
        'membership': a.outdir / f'GO_MF_gene_set_membership_{tag}.tsv',
        'candidates': a.outdir / f'GO_MF_candidates_{tag}_redundancy_audit.tsv',
        'overlaps': a.outdir / f'GO_MF_gene_set_overlaps_{tag}.tsv',
    }
    with open(out['gmt'], 'w', encoding='utf-8', newline='\n') as fh:
        for r in set_rows:
            fh.write('\t'.join([r['GO_ID'], r['GO_term'], *r['genes'].split(';')]) + '\n')
    write_tsv(out['sets'], set_rows, list(set_rows[0]))
    write_tsv(out['membership'], member_rows, list(member_rows[0]))
    write_tsv(out['candidates'], cand_rows, list(cand_rows[0]))
    write_tsv(out['overlaps'], overlap_rows, ['GO_ID_A', 'GO_term_A', 'n_genes_A', 'GO_ID_B', 'GO_term_B', 'n_genes_B',
                                              'shared_genes', 'jaccard', 'overlap_coefficient', 'relation'])

    covered = 0
    for t in final:
        covered |= masks[t]
    covered_candidates = 0
    for t in candidates:
        covered_candidates |= masks[t]
    status_counts = Counter(s for s, _ in decisions.values())
    summary = {
        'GO_MF_terms': len(nodes),
        'annotated_gene_universe': len(genes),
        'universe_excluding_protein_binding_only_genes': informative_universe,
        'candidate_terms_in_size_range': len(candidates),
        'candidate_unique_genes': covered_candidates.bit_count(),
        'removed_exact_duplicates': status_counts['removed_exact_duplicate'],
        'removed_hierarchical_near_duplicates': status_counts['removed_hierarchical_near_duplicate'],
        'final_gene_sets': len(final),
        'final_unique_genes': covered.bit_count(),
        'coverage_pct_of_universe': pct(covered.bit_count(), len(genes)),
        'coverage_pct_excluding_protein_binding_only': pct(covered.bit_count(), informative_universe),
        'size_distribution': size_summary([masks[t].bit_count() for t in final]),
        'min_depth_distribution': dict(sorted(Counter(graph['min_depth'][t] for t in final).items())),
        'overlap': {'pairs_total': len(final) * (len(final) - 1) // 2, **dict(sorted(stats.items()))},
    }
    metadata = {
        'generated_utc': datetime.now(timezone.utc).isoformat(timespec='seconds'),
        'python': platform.python_version(),
        'script_sha256': sha256(Path(__file__)),
        'command': ' '.join([Path(sys.argv[0]).name, *sys.argv[1:]]),
        'sources': {
            'go-basic.obo': {'sha256': sha256(obo), 'bytes': obo.stat().st_size,
                             'data_version': obo_header.get('data-version'),
                             'format_version': obo_header.get('format-version')},
            'goa_human.gaf.gz': {'sha256': sha256(gaf), 'bytes': gaf.stat().st_size,
                                 'gaf_version': gaf_header.get('gaf-version'),
                                 'date_generated': gaf_header.get('date-generated')},
        },
        'parameters': {'min_genes': a.min_genes, 'max_genes': a.max_genes, 'jaccard_threshold': a.jaccard,
                       'overlap_report_thresholds': {'jaccard': REPORT_JACCARD, 'overlap_coefficient': REPORT_OVERLAP},
                       'target_group_count': None},
        'method': {
            'annotations': 'UniProtKB, taxon 9606, aspect F; NOT qualifier and ND evidence removed; all other '
                           'evidence codes kept (incl. IEA); alt_id mapped to primary ID; obsolete terms dropped',
            'propagation': 'is_a edges only; genes propagated to all is_a ancestors',
            'depth': 'min_depth = shortest is_a path from root; max_depth = longest is_a path from root',
            'selection': f'terms with {a.min_genes}-{a.max_genes} genes at any depth; visited by descending '
                         'max_depth, then size, then GO ID; removed if a retained term has an identical gene set, '
                         f'or a retained descendant has Jaccard >= {a.jaccard}; no coverage objective, no group cap',
            'flags': 'pct_* and family columns are informational only and do not affect membership',
        },
        'gaf_qc': dict(sorted(gaf_qc.items())),
        'summary': summary,
        'outputs': {p.name: sha256(p) for p in out.values()},
    }
    (a.outdir / 'run_metadata.json').write_text(json.dumps(metadata, indent=2, ensure_ascii=False) + '\n', encoding='utf-8')

    print(f"GO:MF terms {len(nodes):,} | annotated genes {len(genes):,}")
    print(f"Terms with {tag} genes: {len(candidates)} | removed exact duplicates: {summary['removed_exact_duplicates']} "
          f"| removed hierarchical near-duplicates: {summary['removed_hierarchical_near_duplicates']}")
    print(f"Final gene sets: {len(final)} | unique genes {covered.bit_count():,} ({summary['coverage_pct_of_universe']}%)")
    print(f"Outputs: {a.outdir}")


if __name__ == '__main__':
    main()
