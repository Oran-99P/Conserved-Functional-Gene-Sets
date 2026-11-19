# Paired human–mouse GO:MF gene sets (Ensembl Compara release 116)

This directory holds 270 functional gene groups, each defined in human and mouse under the same GO ID. They are
ready for cross-species scRNA-seq analysis.

## Sources

**Human sets (unchanged, not duplicated here)**

The approved GO:MF collection is available in `../../GO_MF_gene_sets/results/`. The size range of 50–300 human genes per set was applied only when selecting the human GO terms.

| File | Content | Gene identifier |
|---|---|---|
| `GO_MF_gene_sets_50-300.gmt` | Human sets | HGNC gene symbol from GOA |
| `GO_MF_gene_sets_50-300.tsv` | One row per set | HGNC gene symbol from GOA |
| `GO_MF_gene_set_membership_50-300.tsv` | One row per set × gene; adds UniProt accessions | HGNC gene symbol from GOA |

- GO release 2026-07-26; GOA human annotations generated 2026-05-21.
- 270 sets, 11,086 genes.

**Orthology**

All orthology mapping files are located in the parent directory (`../`):

- `human_mouse_orthologs_all_r116.tsv`: all 14,067 human–mouse ortholog pairs from Ensembl Compara release 116. **This is the source of truth for every relationship.**
- `gene_ortholog_status_r116.tsv`: Ensembl gene ID and orthology status for each human gene symbol.
- `GO_set_ortholog_coverage_r116.tsv`: ortholog coverage statistics for each of the 270 GO groups.
- `run_summary.json`: mapping statistics, source-file checksums, and validation results.
- `map_orthologs_r116.py`: reproducible ortholog mapping pipeline.

## How the mouse sets were built

For each GO group:
1. Take every human gene in the original set.
2. Collect **all** of its Compara r116 mouse orthologs: one-to-one, one-to-many and many-to-many.
3. List each mouse gene once, even when several human genes map to it.

The build is deliberately minimal:
- Human genes without a mouse ortholog, or without a reliable Ensembl ID, contribute nothing. No mouse gene is
  substituted for them.
- No filter is applied on orthology confidence, sequence identity, gene family or mouse set size.

**Mouse set sizes differ from human set sizes and are not restricted to 50–300:**

| Result | Value |
|---|---|
| Mouse set size, min / median / max | 41 / 92 / 298 |
| Mouse/human size ratio, min / median / max | 0.76 / 0.99 / 2.31 |
| Groups with 0 mouse genes | none |
| Groups with fewer than 50 mouse genes | 8 |
| Groups with a larger mouse set than human set | 88 |

- **Smaller mouse sets:** these occur when human genes have no mouse ortholog. Examples: chemokine activity, and
  immunoglobulin V segments in antigen binding.
- **Larger mouse sets:** these occur when one human gene has several mouse orthologs. Odorant binding has 74 human
  genes but 171 mouse genes, because the olfactory receptor family is larger in mouse. Other examples are the PRAME,
  serpin and P450 families.

## Output files (this directory)

| File | Content | Gene identifier |
|---|---|---|
| `GO_MF_mouse_gene_sets_50-300_r116.gmt` | Mouse sets in GMT format, one line per GO group (details below) | **Mouse Ensembl gene ID** (`ENSMUSG…`, unversioned, GRCm39 / Ensembl 116) |
| `GO_MF_mouse_gene_sets_50-300_r116.tsv` | One row per GO group: GO ID, term, number of mouse genes, original number of human genes, mouse Ensembl IDs, mouse symbols (same order) | Mouse Ensembl ID; symbols as annotation |
| `GO_MF_mouse_gene_set_membership_50-300_r116.tsv` | One row per group × mouse gene: symbol, biotype, orthology type(s), source human Ensembl IDs and GO symbols | Mouse Ensembl ID; human Ensembl ID and GO symbol for traceability |
| `GO_MF_human_mouse_paired_sets_50-300_r116.tsv` | One row per GO group: human set size, human genes with an ortholog (n, %), distinct mouse genes, one-to-one / non-one-to-one-only / unmapped human genes (split into no ortholog vs unresolved ID) | — |
| `validation_summary.json` | Validation results, size distributions, symbol issues | — |
| `build_mouse_gene_sets.py` | Rebuilds every file above in a few seconds | — |

The `50-300` in the filenames refers to the size range of the **human** selection. It is not a limit on mouse set size.

GMT line layout: `GO_ID <tab> "GO term | mouse Ensembl gene IDs | Compara r116 orthologs of human GO:MF 50-300 set" <tab> ENSMUSG… <tab> ENSMUSG…`

## Which files to use downstream

- **Paired analysis:** human `../../GO_MF_gene_sets/results/GO_MF_gene_sets_50-300.gmt` together with mouse
  `GO_MF_mouse_gene_sets_50-300_r116.gmt`, matched by GO ID.
  - Match mouse datasets by Ensembl gene ID.
  - Human genes are listed by symbol; their Ensembl IDs are in `../gene_ortholog_status_r116.tsv`.
- **Strict gene-to-gene correspondence:** use `../human_mouse_orthologs_one2one_r116.tsv` (9,615 pairs).
- **Coverage per group:** use `GO_MF_human_mouse_paired_sets_50-300_r116.tsv`.

## Summary and validation

**Coverage**
- 30,172 group × mouse-gene memberships; 10,942 unique mouse Ensembl genes.
- 10,233 of 11,086 human genes (92.31%) have at least one mouse ortholog.

**All 11 checks pass** (details in `validation_summary.json`):
- 270 human and 270 mouse groups, with identical GO IDs.
- The human files are unchanged; their SHA-256 values match `../../GO_MF_gene_sets/results/run_metadata.json`.
- Every mouse gene is backed by a release-116 pair.
- All non-one-to-one partners are kept.
- No mouse gene is duplicated within a group.
- The set TSV, membership TSV, GMT and paired summary agree with each other.
- The figures match the earlier ortholog coverage table.

**Mouse gene symbols**
- 28 mouse genes have no gene name in the Ensembl 116 GTF. Most are PRAME-family orthologs. They are kept by Ensembl
  ID with an empty symbol.
- No symbol is shared by different mouse Ensembl IDs.
