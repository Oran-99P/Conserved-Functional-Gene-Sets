# Conserved-Functional-Gene-Sets

A reproducible pipeline for constructing human **Gene Ontology Molecular Function (GO:MF)** gene sets using the official GO ontology and GOA human annotations.

The resulting gene sets are intended for downstream functional analyses, including human–mouse ortholog mapping and cross-species comparisons.

## Quick start

Requires Python 3.10 or later. No external Python packages are needed.

```bash
python3 build_go_gene_sets.py
```

The script reads the source files in `DB/` and generates the gene sets and validation outputs in `results/`.

## Method

Gene sets are constructed using the following procedure:

1. **Ontology:** Parse `go-basic.obo`, using `is_a` relationships within the Molecular Function namespace. Alternative GO IDs are mapped to their primary IDs, and obsolete terms are excluded.
2. **Annotations:** Parse `goa_human.gaf.gz`, retaining human UniProtKB annotations for Molecular Function. Annotations with `NOT` qualifiers or `ND` evidence are excluded. Electronic annotations (`IEA`) are retained.
3. **Propagation:** Propagate each gene annotation to all `is_a` ancestors according to the GO true-path rule.
4. **Size filtering:** Retain GO terms containing **50–200 genes**, regardless of their depth in the ontology.
5. **Redundancy reduction:** Remove exact duplicate gene sets and highly overlapping hierarchical terms (Jaccard ≥ 0.8), retaining the more specific descendant term.

Non-hierarchical overlaps are preserved. The selection is deterministic and does not use greedy coverage optimization or impose a predefined number of gene sets.

## Results

**GO release:** 2026-07-26  
**GOA human annotations:** 2026-05-21

| Metric | Result |
|---|---:|
| Human genes with GO:MF annotations | 18,293 |
| Candidate gene sets (50–200 genes) | 286 |
| Exact duplicates removed | 2 |
| Hierarchical near-duplicates removed | 38 |
| **Final gene sets** | **246** |
| **Unique genes covered** | **9,742** |
| **Gene coverage** | **53.26%** |
| Median gene-set size | 84 |
| Gene-set size range | 50–200 |

### Main output files

Located in `results/`:

| File | Description |
|---|---|
| `GO_MF_gene_sets_50-200.tsv` | Final gene sets with GO IDs, functional descriptions, sizes, annotation metadata, and gene membership |
| `GO_MF_gene_sets_50-200.gmt` | Final gene sets in standard GMT format |
| `GO_MF_gene_set_membership_50-200.tsv` | Gene-level membership table, including UniProt identifiers and annotation evidence |
| `run_metadata.json` | Input versions, parameters, checksums, and execution metadata |

The gene-level membership table is suitable for downstream identifier normalization and human–mouse ortholog mapping. The GMT file can be used with standard gene-set analysis tools.

### Validation files

Located in `results/validation/`:

| File | Description |
|---|---|
| `GO_MF_gene_set_overlaps_50-200.tsv` | Pairwise overlap statistics and GO hierarchical relationships between retained gene sets |
| `GO_MF_candidates_50-200_redundancy_audit.tsv` | Audit of all 286 candidate terms, documenting which were retained or removed and the reasons for removal |

These files provide transparency regarding redundancy reduction and gene-set overlap. They are not required as inputs for downstream ortholog mapping.

## Limitations

- Shared GO:MF annotations indicate a common annotated molecular function, but do not necessarily imply co-expression or conserved function across species.
- GO annotation coverage and experimental evidence vary between genes and functional categories.
- Some gene sets contain lineage-expanded gene families that may have limited one-to-one human–mouse ortholog representation.
- Gene-set overlap is retained where biologically justified and should be considered when interpreting downstream analyses.
- Gene symbols and UniProt identifiers should be normalized to stable gene identifiers before cross-species mapping.

No gene sets are excluded based on predicted ortholog availability or expression detectability. These properties will be evaluated in downstream analyses.

## Source data and attribution

The `DB/` directory contains the original source files used to generate the results.

| File | Version | Source |
|---|---|---|
| `go-basic.obo` | GO release 2026-07-26 | [Gene Ontology](https://purl.obolibrary.org/obo/go/go-basic.obo) |
| `goa_human.gaf.gz` | Generated 2026-05-21 | [GOA human annotations](https://current.geneontology.org/annotations/goa_human.gaf.gz) |

The source URLs may point to newer releases. The included files and their SHA-256 checksums, recorded in `run_metadata.json`, identify the exact inputs used for this analysis.

Gene Ontology data are provided by the Gene Ontology Consortium under the [Creative Commons Attribution 4.0 International License](https://creativecommons.org/licenses/by/4.0/).

**References**

- Ashburner M. et al. (2000). Gene Ontology: tool for the unification of biology. *Nature Genetics*, 25, 25–29.
- The Gene Ontology Consortium (2023). The Gene Ontology knowledgebase in 2023. *Genetics*, 224, iyad031.