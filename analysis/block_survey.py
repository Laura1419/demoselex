#!/usr/bin/env python3
"""[inference] Exploratory fixed-tile survey; protocol: docs/block-survey-prereg.md.

Gene proximity is not causal attribution. Akbari labels share genotype inputs.
No linked block null or modern functional evidence is available locally.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import pickle
import re
import sys
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

TILE = 250_000
NEAREST_LIMIT = 100_000
ANNOTATION = Path('data/raw/annotation/hg19_refGene.txt.gz')
ANNOTATION_META = Path('data/raw/annotation/hg19_refGene.provenance.json')
URL = 'https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/refGene.txt.gz'
PDF = Path('data/raw/akbari/41586_2026_10358_MOESM1_ESM.pdf')
WORKBOOK = Path('data/raw/akbari/41586_2026_10358_MOESM3_ESM.xlsx')
CACHE = Path('data/derived/tier2_pool_cache.pkl')
OUT = Path('data/derived')
CONTROLS = {
    'LCT/MCM6': {'LCT', 'MCM6'}, 'SLC24A5': {'SLC24A5'},
    'SLC45A2': {'SLC45A2'}, 'TYR': {'TYR'}, 'OCA2/HERC2': {'OCA2', 'HERC2'},
    'FADS1/FADS2': {'FADS1', 'FADS2'}, 'TLR1/TLR6': {'TLR1', 'TLR6'},
    'AMY1/AMY2 region': {'AMY1A', 'AMY1B', 'AMY1C', 'AMY2A', 'AMY2B'},
}
PRIOR_GENES = set.union(*CONTROLS.values())
SHARDS = [Path(f'data/derived/prio_panels/shard{i}_out/priority_loci.tsv') for i in range(4)]


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for part in iter(lambda: f.read(1024 * 1024), b''):
            h.update(part)
    return h.hexdigest()


def ensure_annotation(fetch=False):
    if not ANNOTATION.exists():
        if not fetch:
            raise FileNotFoundError('Annotation unavailable; rerun with --fetch-annotation. No coordinate guesses.')
        ANNOTATION.parent.mkdir(parents=True, exist_ok=True)
        with urllib.request.urlopen(URL, timeout=90) as response:
            content = response.read()
        raw = gzip.decompress(content)  # Verify CRC/truncation before writing raw data.
        ANNOTATION.write_bytes(content)
        ANNOTATION_META.write_text(json.dumps(dict(
            url=URL, retrieved_at_utc=datetime.now(timezone.utc).isoformat(),
            sha256=sha256(ANNOTATION), gzip_verified=True,
            uncompressed_bytes=len(raw), compressed_bytes=len(content),
            coordinate_system='hg19/GRCh37, zero-based half-open txStart/txEnd'), indent=2) + '\n')
    metadata = json.loads(ANNOTATION_META.read_text())
    if metadata['sha256'] != sha256(ANNOTATION) or metadata['url'] != URL:
        raise ValueError('Annotation provenance mismatch; stop rather than guess.')
    with gzip.open(ANNOTATION, 'rb') as f:
        f.read()  # Verify gzip on each run, without editing immutable raw files.
    return metadata


def read_refgene_transcripts(path):
    """[observation] UCSC refGene: chrom=2, strand=3, txStart=4, txEnd=5, name2=12."""
    raw = pd.read_csv(path, sep='\t', header=None, usecols=[2, 3, 4, 5, 12],
                      names=['chrom', 'strand', 'start', 'end', 'gene_symbol'])
    raw = raw[raw.chrom.str.fullmatch(r'chr(?:[1-9]|1[0-9]|2[0-2])')].copy()
    raw['chrom'] = raw.chrom.str[3:].astype(int)
    if raw.gene_symbol.isna().any() or (raw.end <= raw.start).any():
        raise ValueError('Invalid refGene span/symbol')
    return raw


def merge_transcript_intervals(raw):
    """[inference] Approved correction: zero-gap interval union, never bounding gaps."""
    records = []
    for (symbol, chrom, strand), group in raw.groupby(['gene_symbol', 'chrom', 'strand']):
        start = end = None
        count = 0
        for row in group.sort_values(['start', 'end']).itertuples():
            if start is not None and row.start > end:
                records.append(dict(gene_symbol=symbol, chrom=chrom, strand=strand,
                                    start=start, end=end, n_transcripts=count))
                start = end = None
                count = 0
            if start is None:
                start, end = int(row.start), int(row.end)
            else:
                end = max(end, int(row.end))
            count += 1
        if start is not None:
            records.append(dict(gene_symbol=symbol, chrom=chrom, strand=strand,
                                start=start, end=end, n_transcripts=count))
    merged = pd.DataFrame(records).sort_values(['chrom', 'start', 'end', 'gene_symbol', 'strand']).reset_index(drop=True)
    merged['excluded_pathological_span'] = (merged.end - merged.start).gt(4_000_000)
    return merged


def parse_refgene(path):
    return merge_transcript_intervals(read_refgene_transcripts(path))


def build_blocks(variants):
    variants = variants.copy()
    if variants.locus.duplicated().any():
        raise ValueError('Duplicate scored locus: refusing silent double-counting')
    if not variants.chrom.between(1, 22).all() or not variants.pos.ge(1).all():
        raise ValueError('Expected one-based autosomal coordinates')
    if not np.isfinite(variants.total_score).all():
        raise ValueError('Nonfinite scores: cannot silently change scored denominator')
    variants['block_start'] = ((variants.pos - 1) // TILE) * TILE
    variants['abs_X'] = variants.akbari_X.abs()
    ordered = variants.sort_values(['total_score', 'locus'], ascending=[False, True])
    reps = ordered.drop_duplicates(['chrom', 'block_start'])
    blocks = reps[['chrom', 'block_start', 'locus', 'pos', 'total_score']].rename(
        columns={'locus': 'representative_locus', 'pos': 'representative_pos'})
    summary = variants.groupby(['chrom', 'block_start'], as_index=False).agg(
        best_abs_X=('abs_X', 'max'), n_variants_scored_in_block=('locus', 'size'),
        n_variants_missing_X=('abs_X', lambda a: int(a.isna().sum())))
    blocks = blocks.merge(summary, on=['chrom', 'block_start'], validate='one_to_one')
    blocks = blocks.sort_values(['total_score', 'representative_locus'], ascending=[False, True]).reset_index(drop=True)
    blocks['block_end'] = blocks.block_start + TILE
    blocks['rank'] = np.arange(1, len(blocks) + 1)
    for name, threshold, strict in [('strong', 5.45, True), ('moderate', 3.16, False)]:
        flag = (blocks.best_abs_X.gt(threshold) if strict else blocks.best_abs_X.ge(threshold)).astype('boolean')
        # If any missing X could change a false result, it is unknown, not negative.
        flag.loc[~flag & blocks.n_variants_missing_X.gt(0)] = pd.NA
        blocks[f'is_akbari_{name}'] = flag
    blocks['akbari_label'] = np.select(
        [blocks.is_akbari_strong.fillna(False).to_numpy(bool),
         blocks.is_akbari_moderate.fillna(False).to_numpy(bool)], ['strong', 'moderate_only'], default='neither')
    blocks.loc[blocks.is_akbari_moderate.isna(), 'akbari_label'] = 'unknown'
    blocks['is_mhc'] = blocks.chrom.eq(6) & blocks.block_start.lt(33_500_000) & blocks.block_end.gt(28_500_000)
    return blocks


def annotate_blocks(blocks, genes):
    annotations = []
    if 'excluded_pathological_span' in genes:
        genes = genes[~genes.excluded_pathological_span]
    by_chrom = {c: g for c, g in genes.groupby('chrom')}
    for row in blocks.itertuples():
        g = by_chrom.get(row.chrom, genes.iloc[:0])
        overlap = g[g.start.lt(row.block_end) & g.end.gt(row.block_start)]
        symbols = sorted(overlap.gene_symbol.unique())
        nearest = None
        distance = None
        if not symbols and len(g):
            gaps = np.maximum(np.maximum(g.start - row.block_end, row.block_start - g.end), 0)
            choices = g.assign(distance=gaps).sort_values(['distance', 'gene_symbol', 'start', 'end', 'strand'])
            best = choices.iloc[0]
            if best.distance <= NEAREST_LIMIT:
                nearest, distance = best.gene_symbol, int(best.distance)
        prior = sorted((set(symbols) | ({nearest} if nearest else set())) & PRIOR_GENES)
        annotations.append(dict(genes_overlapping=';'.join(symbols) or None,
            nearest_gene=nearest, nearest_gene_distance=distance,
            annotation_status='overlap' if symbols else ('nearest_within_100kb' if nearest else 'no_gene_within_100kb'),
            prior_selection_genes=';'.join(prior) or None,
            has_prior_selection_gene=bool(prior),
            prior_selection_relation=('overlap' if symbols else 'nearest_context_only') if prior else 'no_match',
            prior_selection_evidence=('task_fixed_canonical_panel_overlap' if symbols else 'nearest_context_not_block_support') if prior else 'no_match_in_fixed_panel'))
    return pd.concat([blocks.reset_index(drop=True), pd.DataFrame(annotations)], axis=1)


def extract_independent_loci():
    """[observation] PDF pages 52–91: S5.1–S5.40, the 479 trajectory tag coordinates.

    No LD interval is inferred from these tag positions; the 347 pruned set is
    not reconstructed from the count-only Table S5.1.
    """
    from pypdf import PdfReader
    reader = PdfReader(PDF)
    texts = [page.extract_text() for page in reader.pages]
    (OUT / 'block_survey_supplement_text.txt').write_text('\n'.join(
        f'\n--- PAGE {i+1} ---\n{text}' for i, text in enumerate(texts)))
    pattern = re.compile(r'\n(\d+)\s*\n\s*\|X\|\s*=\s*([\d.]+)\s*\n([^\n]+)\nchr(\d+):(\d+)\s*\(([^>]+)>([^)]+)\)')
    records = []
    for i in range(51, 91):
        if f'Supplementary Figure S5.{i - 50}:' not in texts[i]:
            raise ValueError('Unexpected supplement page/figure mapping')
        for index, abs_x, label, chrom, pos, ref, alt in pattern.findall(texts[i]):
            records.append(dict(trajectory_index=int(index), chrom=int(chrom), pos=int(pos),
                pdf_abs_X=float(abs_x), pdf_label=label.strip(), pdf_ref=ref, pdf_alt=alt, pdf_page=i+1))
    loci = pd.DataFrame(records)
    if len(loci) != 479 or set(loci.trajectory_index) != set(range(1, 480)) or loci.duplicated(['chrom', 'pos']).any():
        raise ValueError('Incomplete/ambiguous 479-tag extraction')
    # Local cache is derived from GRCh37 summary statistics. Rounded X consistency
    # checks all coordinates, including indels with ambiguous printed alleles.
    with open(CACHE, 'rb') as f:
        pool = pickle.load(f)
    pool = pool.rename(columns={'CHROM': 'chrom', 'POS': 'pos', 'ID': 'cached_locus', 'X': 'cached_X'})
    # Restrict to extracted coordinates before checking uniqueness: the cache
    # contains multiallelic positions elsewhere, but each extracted tag is unique.
    pool = pool.merge(loci[['chrom', 'pos']], on=['chrom', 'pos'], how='inner', validate='many_to_one')
    loci = loci.merge(pool[['chrom', 'pos', 'cached_locus', 'cached_X']], on=['chrom', 'pos'], validate='one_to_one', how='left')
    if loci.cached_locus.isna().any() or (loci.cached_X.abs() - loci.pdf_abs_X).abs().gt(.015).any():
        raise ValueError('PDF tag coordinate/X check against local GRCh37 cache failed')
    loci['block_start'] = ((loci.pos - 1) // TILE) * TILE
    return loci


def add_evidence(blocks, loci):
    blocks = blocks.copy()
    if loci is None:
        blocks['akbari_479_tag_in_block'] = pd.NA
        blocks['representative_is_akbari_479_tag'] = pd.NA
        blocks['akbari_479_tag_loci'] = None
        blocks['akbari_independent_evidence_status'] = 'unavailable_extraction_failed'
    else:
        tags = {(c, s): ';'.join(g.cached_locus.sort_values()) for (c, s), g in loci.groupby(['chrom', 'block_start'])}
        ids = set(loci.cached_locus)
        blocks['akbari_479_tag_loci'] = [tags.get((r.chrom, r.block_start)) for r in blocks.itertuples()]
        blocks['akbari_479_tag_in_block'] = blocks.akbari_479_tag_loci.notna()
        blocks['representative_is_akbari_479_tag'] = blocks.representative_locus.isin(ids)
        blocks['akbari_independent_evidence_status'] = '479_tags_only_no_LD_intervals_or_347_subset'
    blocks['akbari_347_locus_overlap'] = pd.NA
    blocks['table_s5_locus_phenotype'] = None
    blocks['table_s5_evidence_status'] = 'polygenic_aggregate_only_no_locus_mapping'
    for name in ['gwas_catalog', 'eqtl', 'functional_annotation']:
        blocks[f'{name}_evidence_status'] = 'unavailable_locally_not_checked'
    return blocks


def cache_evidence(blocks, variants):
    """[inference] Post-registration evidence only: no reranking or moderate enrichment."""
    with open(CACHE, 'rb') as f:
        cache = pickle.load(f)
    cache = cache.rename(columns={'CHROM': 'chrom', 'POS': 'pos', 'X': 'akbari_X'})
    blocks = blocks.copy()
    for source, frame in [('scored', variants), ('cache', cache)]:
        frame = frame.copy()
        frame['block_start'] = ((frame.pos - 1) // TILE) * TILE
        frame['abs_X'] = frame.akbari_X.abs()
        frame['strong'] = frame.abs_X.gt(5.45)
        frame['moderate_band'] = frame.abs_X.ge(3.16) & frame.abs_X.lt(5.45)
        grouped = frame.groupby(['chrom', 'block_start'], as_index=False).agg(
            maximum=('abs_X', 'max'), strong=('strong', 'any'),
            moderate=('moderate_band', 'any'), count=('abs_X', 'size'))
        grouped = grouped.rename(columns={'maximum': f'block_max_absX_{source}',
            'strong': f'block_has_strong_{source}', 'moderate': f'block_has_moderate_{source}',
            'count': f'n_{source}_variants_in_block'})
        blocks = blocks.merge(grouped, on=['chrom', 'block_start'], how='left', validate='one_to_one')
        for flag in ['strong', 'moderate']:
            blocks[f'block_has_{flag}_{source}'] = blocks[f'block_has_{flag}_{source}'].fillna(False).astype(bool)
        blocks[f'block_has_any_akbari_{source}'] = blocks[f'block_has_strong_{source}'] | blocks[f'block_has_moderate_{source}']
        blocks[f'n_{source}_variants_in_block'] = blocks[f'n_{source}_variants_in_block'].fillna(0).astype(int)
    blocks['cache_evidence_status'] = 'post_registration_QC_cache_annotation_only_no_rescoring'
    return blocks, dict(cache_variants=len(cache), cache_strong_variants=int(cache.akbari_X.abs().gt(5.45).sum()),
        cache_moderate_band_variants=int((cache.akbari_X.abs().ge(3.16) & cache.akbari_X.abs().lt(5.45)).sum()),
        cache_exact_5_45_variants=int(cache.akbari_X.abs().eq(5.45).sum()),
        no_evidence_meaning='No >=3.16 evidence in available QC cache, not proof of neutrality or absence from all data')


def annotation_correction_diagnostics(blocks, genes):
    baseline = pd.read_csv(OUT / 'block_survey_initial_bounding.tsv', sep='\t', float_precision='round_trip')
    old = baseline.set_index(['chrom', 'block_start'])
    new = blocks.set_index(['chrom', 'block_start']).loc[old.index]
    changed = pd.Series(False, index=old.index)
    for col in ['genes_overlapping', 'nearest_gene', 'nearest_gene_distance']:
        changed |= old[col].fillna('NA').ne(new[col].fillna('NA'))
    previous_novel = old.akbari_label.eq('neither')
    previous_reported = old[previous_novel].head(10).index
    losses = []
    for index in previous_reported:
        before_value, after_value = old.loc[index, 'genes_overlapping'], new.loc[index, 'genes_overlapping']
        before = set() if pd.isna(before_value) else set(before_value.split(';'))
        after = set() if pd.isna(after_value) else set(after_value.split(';'))
        losses.append(dict(chrom=int(index[0]), block_start=int(index[1]), rank=int(old.loc[index, 'rank']),
            lost_genes=sorted(before - after), gained_genes=sorted(after - before),
            before_genes=sorted(before), after_genes=sorted(after)))
    groups = genes.groupby(['gene_symbol', 'chrom', 'strand']).size()
    pathological = genes[genes.excluded_pathological_span]
    return dict(merge_gap_tolerance_bp=0,
        disjoint_multi_interval_symbol_chrom_strand_groups=int(groups.gt(1).sum()),
        symbols_with_disjoint_intervals=int(len({index[0] for index in groups[groups.gt(1)].index})),
        pathological_length_limit_bp=4_000_000,
        excluded_pathological_symbols=sorted(pathological.gene_symbol.unique()),
        excluded_pathological_interval_count=len(pathological),
        changed_blocks=int(changed.sum()), changed_registered_neither_blocks=int((changed & previous_novel).sum()),
        changed_previously_reported_novel_blocks=int(changed.loc[previous_reported].sum()),
        previously_reported_novel_annotation_delta=losses,
        mhc_flags_unchanged=bool(old.is_mhc.equals(new.is_mhc)),
        ranks_scores_unchanged=bool(old['rank'].equals(new['rank']) and np.allclose(old.total_score, new.total_score, rtol=0, atol=0)))


def phenotype_context():
    frame = pd.read_excel(WORKBOOK, sheet_name='STable5(696.polygenic.tests)', header=1)
    column = 'gamma_z (beta without HLA)'
    z = pd.to_numeric(frame[column], errors='coerce')
    columns = ['phenotype ID', 'phenotype name', 'coding_description (Pan-UKBB)', column, 'gamma_z (beta with HLA)']
    strongest = frame.loc[z.abs().sort_values(ascending=False).head(5).index, columns]
    return dict(rows=len(frame), locus_mapping_available=False,
        summary_statistic=column, strongest_absolute_z=strongest.replace({np.nan: None}).to_dict('records'),
        interpretation='Aggregate polygenic tests only; no block/gene-to-phenotype attribution')


def simulation_diagnostics(blocks):
    """[inference] Existing independent-locus draws cannot form a linked-block null."""
    result = dict(block_level_null_rate=None,
        gap='Simulated columns have no genomic coordinates or linked haplotypes; block maxima/multiplicity cannot be calibrated.',
        per_variant_tail_at_block_threshold=[], legacy={})
    metrics = json.loads(Path('data/derived/prioritization/metrics.json').read_text())
    result['legacy'] = dict(scored_variants=metrics['scored_variants'],
        random_walk_jump_top100=next(x for x in metrics['drift_simulation']['random_walk_jump']['rates'] if x['N'] == 100),
        bootstrap_top100=next(x for x in metrics['bootstrap']['metrics'] if x['N'] == 100),
        bootstrap_replicates=metrics['bootstrap']['replicates'])
    for mechanism in ['random_walk', 'random_walk_jump', 'pure_replacement']:
        path = Path(f'data/derived/prioritization/simulation_{mechanism}.tsv')
        draws = pd.read_csv(path, sep='\t')
        draws = draws[draws.valid]
        for n in [100, 500, 1000]:
            threshold = float(blocks.iloc[n-1].total_score)
            exceedances = int(draws.total_score.ge(threshold).sum())
            result['per_variant_tail_at_block_threshold'].append(dict(mechanism=mechanism, N=n,
                threshold=threshold, draws=len(draws), exceedances=exceedances,
                rate=exceedances / len(draws) if len(draws) else None))
    return result


def summarize(blocks, variants, annotation_metadata, independent_status):
    counts = []
    for n in [100, 500, 1000]:
        top = blocks.head(n)
        strong, moderate = int(top.is_akbari_strong.sum()), int(top.is_akbari_moderate.sum())
        counts.append(dict(N=n, strong=strong, moderate_including_strong=moderate,
            moderate_only=moderate-strong, non_strong=int(top.is_akbari_strong.eq(False).sum()),
            neither=int(top.akbari_label.eq('neither').sum()), unknown=int(top.akbari_label.eq('unknown').sum()),
            strong_enrichment=strong / len(top) / blocks.is_akbari_strong.mean(),
            moderate_enrichment=moderate / len(top) / blocks.is_akbari_moderate.mean(),
            mhc=int(top.is_mhc.sum()), prior_gene_blocks=int(top.has_prior_selection_gene.sum()),
            tag_blocks=int(top.akbari_479_tag_in_block.sum()),
            cache_strong=int(top.block_has_strong_cache.sum()),
            cache_moderate_only=int((~top.block_has_strong_cache & top.block_has_moderate_cache).sum()),
            cache_no_evidence=int((~top.block_has_any_akbari_cache).sum()),
            registered_novel_with_cache_moderate=int((top.is_akbari_strong.eq(False) & top.block_has_moderate_cache).sum()),
            registered_novel_with_no_cache_evidence=int((top.is_akbari_strong.eq(False) & ~top.block_has_any_akbari_cache).sum())))
    return dict(annotation_provenance=annotation_metadata, scored_variants=len(variants),
        strong_variants=int(variants.akbari_X.abs().gt(5.45).sum()), scored_blocks=len(blocks),
        strong_blocks=int(blocks.is_akbari_strong.sum()), moderate_blocks=int(blocks.is_akbari_moderate.sum()),
        strong_block_prevalence=float(blocks.is_akbari_strong.mean()),
        moderate_block_prevalence=float(blocks.is_akbari_moderate.mean()),
        annotation_status_counts=blocks.annotation_status.value_counts().to_dict(),
        independent_locus_status=independent_status, top_N=counts, simulation=simulation_diagnostics(blocks))


def gene_display(row):
    if pd.notna(row.genes_overlapping):
        return row.genes_overlapping
    if pd.notna(row.nearest_gene):
        return f'{row.nearest_gene} (nearest, {int(row.nearest_gene_distance)} bp)'
    return 'null (no gene within 100 kb)'


def candidate_table(frame):
    lines = ['|Rank|GRCh37 tile [start,end)|Representative|Score|Best abs X|Scored label|MHC|Genes (overlap unless stated)|Prior panel|479 tag in tile / rep|Cache evidence (POST-REG)|',
             '|---:|---|---|---:|---:|---|---|---|---|---|---|']
    for r in frame.itertuples():
        tag = 'unknown' if pd.isna(r.akbari_479_tag_in_block) else ('yes' if r.akbari_479_tag_in_block else 'no')
        rep = 'unknown' if pd.isna(r.representative_is_akbari_479_tag) else ('yes' if r.representative_is_akbari_479_tag else 'no')
        prior = r.prior_selection_genes if pd.notna(r.prior_selection_genes) else 'no fixed-panel match'
        cache = 'strong' if r.block_has_strong_cache else ('moderate_only' if r.block_has_moderate_cache else 'none_in_QC_cache')
        lines.append(f'|{r.rank}|chr{r.chrom}:[{r.block_start},{r.block_end})|{r.representative_locus}|{r.total_score:.3f}|{r.best_abs_X:.3f}|{r.akbari_label}|{"yes" if r.is_mhc else "no"}|{gene_display(r)}|{prior}|{tag} / {rep}|{cache}|')
    return '\n'.join(lines)


def write_results(blocks, summary, phenotype):
    target = Path('docs/block-survey-results.md')
    lines = ['# Exploratory block-level gene survey — 2026-10-09', '',
        '[inference] **An exploratory hypothesis list for downstream modern-phenotype follow-up, NOT a selection discovery.** False positives are expected. Gene proximity and Akbari agreement do not establish function, causality, or independent validation.', '',
        '## Registration and reproducibility', '',
        '[observation] Protocol persisted first in `docs/block-survey-prereg.md` at 2026-10-09T02:18:36Z; annotation retrieval followed at '+summary['annotation_provenance']['retrieved_at_utc']+'. Inputs are the four scored shard TSVs, not the earlier small panel. Reproduce with `.venv/bin/python analysis/block_survey.py` (add `--fetch-annotation` only if the raw annotation is absent). Explicit TSV null token: `NA`.', '',
        '[observation] UCSC hg19 refGene URL: '+URL+'; gzip CRC verified; SHA256 `'+summary['annotation_provenance']['sha256']+'`. Provenance: `data/raw/annotation/hg19_refGene.provenance.json`. Parsed symbol/chromosome/interval-union/strand: `data/derived/block_survey_gene_spans.tsv`. All gene names in the candidate/control tables below are generated from annotated rows of `data/derived/block_survey.tsv`, never from PDF gene labels.', '',
        f"[observation] {summary['scored_variants']:,} finite scored variants, {summary['strong_variants']:,} strong variants; {summary['scored_blocks']:,} occupied 250 kb tiles, {summary['strong_blocks']:,} strong tiles and {summary['moderate_blocks']:,} moderate-including-strong tiles. No duplicate loci or missing input X. Annotation statuses: `{summary['annotation_status_counts']}`.", '',
        f"[observation] The task's 181,933 and 31,939 variant counts agree. Approximately 723 strong blocks compares with the measured {summary['strong_blocks']:,}; tiling is explicitly one-based-pos to zero-based-start. Tiles are LD-reduced, not independent: adjacent tiles and MHC can share one signal.", '',
        '## Top 20 blocks and local evidence', '',
        '[observation] Strong means any scored variant |X|>5.45; moderate_only means max |X|>=3.16 but not strong. Moderate flags in TSV include strong. Prior panel = task-supplied canonical ancient-DNA support panel, also represented by the existing prioritizer controls; sparse illustrative controls, not an exhaustive literature review. A non-match is not evidence of no prior support. Coordinates/gene spans are zero-based half-open; representative locus strings contain one-based positions.', '',
        candidate_table(blocks.head(20)), '',
        '[inference] Overlapping genes are positional candidates only (union of actual transcript spans, not coding-exon intersection). Nearest fallback is allowed only when no gene overlaps, at <=100 kb. Nearest context is NOT positive evidence for a block. MHC flag = overlap chr6:[28500000,33500000). Multiple high-ranked MHC tiles remain one extended LD complex, not multiple independent biological confirmations.', '',
        '## Prior ancient-DNA support panel: what surfaced?', '',
        '|Control group|Best annotated block rank|Block representative|Annotated panel genes in that row|Label|',
        '|---|---:|---|---|---|']
    all_found = set()
    for group, names in CONTROLS.items():
        hits = blocks[blocks.prior_selection_genes.fillna('').map(lambda x: bool(set(x.split(';')) & names))]
        if hits.empty:
            # Do not name untraceable symbols in the results; group name only describes protocol.
            lines.append(f'|Protocol control group (see preregistration)|NA|NA|not annotated in scored tiles|unknown|')
        else:
            row = hits.iloc[0]
            found = ';'.join(sorted(set(row.prior_selection_genes.split(';')) & names))
            all_found.update(found.split(';'))
            lines.append(f'|{group}|{row["rank"]}|{row.representative_locus}|{found}|{row.akbari_label}|')
    top100 = summary['top_N'][0]
    lines += ['',
        f"[observation] Pre-registered predictions: top-100 strong enrichment is {top100['strong_enrichment']:.2f}× (prediction >1: {'met' if top100['strong_enrichment'] > 1 else 'falsified'}); {top100['prior_gene_blocks']} top-100 tiles have a fixed-panel gene (prediction >=1: {'met' if top100['prior_gene_blocks'] else 'falsified'}). Neither success supplies independent validation.", '',
        '[observation] Inherited variant-level limits remain: SLC24A5 ranked approximately 9,161 and the AMY1 SNP proxy approximately 698 in the earlier 10,641-variant panel. The table above reports NEW block ranks in the much denser merged survey and does not replace those prior variant ranks. The within-lineage residual design is structurally blind to carriage-driven change predominantly BETWEEN demes; it should not be presented as capturing those mechanisms. The amylase region SNP proxy is not a copy-number selection test.', '',
        '## Akbari reported independent loci and phenotype evidence', '',
        '[observation] '+summary['independent_locus_status'], '',
        '[inference] “479 tag in tile” checks whether an extracted tag coordinate lies in the tile, and “rep” checks exact representative locus ID against its verified cached tag ID. It is NOT overlap with an LD-defined independent-locus interval. Such intervals and membership of the more strictly pruned 347 set are not available from this extraction; `akbari_347_locus_overlap=NA` throughout. Absence of a tag from a tile does not imply independence or absence of the reported signal. Section 5 describes 479 tags (410 outside HLA +69 inside); Table S5.1 supplies D-prime pruning COUNTS, not the identities of a 347-tag subset. This is a difference in list/definition, not a claim that 347 is incorrect.', '',
        f"[observation] Table S5 sheet `STable5(696.polygenic.tests)` contains {phenotype['rows']} aggregate phenotype tests, with beta/sign-beta analyses with and without HLA. It contains no locus or gene mapping. The strongest five absolute z values in the beta-without-HLA column are shown solely to describe available context, not to connect these phenotypes to the ranked tiles:", '',
        '|Phenotype ID|Phenotype / coding description|z beta without HLA|z beta with HLA|',
        '|---|---|---:|---:|']
    for r in phenotype['strongest_absolute_z']:
        coding = '' if r['coding_description (Pan-UKBB)'] is None else ' / '+str(r['coding_description (Pan-UKBB)'])
        lines.append(f"|{r['phenotype ID']}|{r['phenotype name']}{coding}|{r['gamma_z (beta without HLA)']:.3f}|{r['gamma_z (beta with HLA)']:.3f}|")
    lines += ['',
        '[inference] Aggregate Table S5 results cannot adjudicate any candidate gene; all locus-phenotype cells are explicitly null. **GWAS Catalog, eQTL, and functional annotation are NOT available locally and were not checked.** Their evidence statuses are explicit in every TSV row. RefGene transcript spans are positional annotations, not functional evidence.', '',
        '## False-positive / follow-up cost at block level', '',
        f"[observation] Base prevalence in the {summary['scored_blocks']:,} scored tiles: strong {summary['strong_blocks']:,}/{summary['scored_blocks']:,} = {100*summary['strong_block_prevalence']:.2f}%; moderate including strong {summary['moderate_blocks']:,}/{summary['scored_blocks']:,} = {100*summary['moderate_block_prevalence']:.2f}%. Denominator is the ascertained scored panel, not all genomic tiles.", '',
        '|Top N|Strong (count / N)|Strong enrichment|Moderate incl. strong (count / N)|Moderate enrichment|Moderate only|Non-strong|Neither (strict novel)|MHC|479 tag tiles|',
        '|---:|---|---:|---|---:|---:|---:|---:|---:|---:|']
    for r in summary['top_N']:
        lines.append(f"|{r['N']}|{r['strong']}/{r['N']} ({100*r['strong']/r['N']:.1f}%)|{r['strong_enrichment']:.2f}×|{r['moderate_including_strong']}/{r['N']} ({100*r['moderate_including_strong']/r['N']:.1f}%)|{r['moderate_enrichment']:.2f}×|{r['moderate_only']}|{r['non_strong']}|{r['neither']}|{r['mhc']}|{r['tag_blocks']}|")
    sim = summary['simulation']; legacy = sim['legacy']
    cost = legacy['random_walk_jump_top100']
    retention = legacy['bootstrap_top100']['mean_original_top_membership']
    lines += ['',
        '[inference] Counts not strong are novel relative to the strong threshold; counts neither are the stricter novel subset, used below. These are follow-up costs, NOT proven false positives: non-significance is not neutrality. Strong and moderate agreement reuses Akbari data (same imputed genotypes), so enrichment is NOT independent validation; original variant ascertainment and tile maxima further favor strong/dense regions.', '',
        f"[observation] A calibrated LD-honest **block drift-null FP rate could not be computed**. Existing `analysis/prioritization_score.py::drift_simulation` / `--simulations` generates independent simulated single-locus genotype columns with no genomic coordinates or linked haplotypes. Arbitrary grouping would invent an LD model. The saved 300-draw drift+jump tail above the EARLIER top-100 variant threshold is {cost['exceedances']}/300 = {100*cost['rate']:.2f}% (agrees with quoted 10.67%). That legacy run scored {legacy['scored_variants']:,} variants, NOT 181,933; all four merged shard runs used zero simulation replicates. Therefore the quoted value is not a merged-survey or block FDR.", '',
        '[observation] Reusing the existing simulation draws, the following SINGLE-LOCUS tails exceed each observed BLOCK-rank score threshold. This is only a threshold diagnostic; it omits within-tile multiplicity and LD, so it is NOT the requested calibrated block false-positive probability or list precision:', '',
        '|Top-block N|Block score threshold|Brownian drift only exceedances / 300|Drift + jump exceedances / 300|Replacement-only exceedances / 300|',
        '|---:|---:|---|---|---|']
    for n in [100, 500, 1000]:
        rows = [r for r in sim['per_variant_tail_at_block_threshold'] if r['N'] == n]
        costs = [f"{r['exceedances']}/{r['draws']} ({100*r['rate']:.2f}%)" for r in rows]
        lines.append(f"|{n}|{rows[0]['threshold']:.4f}|"+'|'.join(costs)+'|')
    lines += ['',
        '[inference] Finite simulation size matters: even 0/300 has a one-sided binomial 95% upper tail probability about 0.994%, not zero. No observed-block FP count is inferable from these single-locus rates. A linked genotype/path simulation with observed block sizes and missingness would be needed to estimate block maxima under drift.', '',
        f"[observation] Earlier within-node bootstrap top-100 variant membership retention = {100*retention:.2f}% across {legacy['bootstrap_replicates']} resamples (agrees with approximately 43.7%). This has not been recomputed at block level or on the merged panel; it warns that ranks near cuts are unstable, not that block retention equals that number.", '',
        '## Highest-ranked novel blocks: bounded follow-up hypotheses', '',
        '[observation] The first ten tiles with neither strong nor moderate scored-variant X are listed mechanically (not selected by appealing gene stories). “Interesting” here means rank/proximity under the pre-registered exploratory protocol only. The broader non-strong candidate set additionally includes moderate-only tiles.', '',
        candidate_table(blocks[blocks.akbari_label.eq('neither')].head(10)), '',
        '[speculation] These tiles could motivate a small modern-phenotype follow-up pilot. Their ranks may instead reflect drift/jitter, carriage residuals, LD with a neighboring signal, or dense-panel maxima; **the local evidence base cannot adjudicate them**. Novel means unlabeled in this scored panel, not globally unreported selection. No causal gene or trait association is claimed.', '',
        '## Validation and bounded recommendation', '',
        '[observation] `analysis/test_block_survey.py` covers tile-boundary/tie rules, any-variant labels and missing-X unknowns, MHC boundaries, overlap/nearest distance limits, and output provenance/gene traceability. Run `.venv/bin/python -m unittest analysis.test_block_survey -v` after producing the outputs. Run manifest hashes the protocol, script, all shard tables, raw annotation/provenance, PDF, workbook, cache and saved simulations; git metadata intentionally disabled because the operator prohibits all git commands.', '',
        '[inference] Use this table as an exploratory, LD-reduced queue, not a selection result or validated production discovery ranking. Budget explicitly for the non-Akbari counts above, cross-tile LD and the uncalibrated block drift cost. First downstream step is independent modern-phenotype evidence acquisition and a linked block-null/stability audit before broad follow-up. Earlier authors\' conclusions are untouched.', '']
    correction = summary['annotation_correction']
    cache = summary['post_registration_cache_evidence']
    lines += ['## POST-REGISTRATION CORRECTIONS AND EVIDENCE SENSITIVITY', '',
        '[observation] Approved after the initial registered run, at approximately 2026-10-09T02:28Z. Earlier text above this appended survey is historical, NOT the current annotation. Initial TSV/spans/summary are retained as `data/derived/block_survey_initial_bounding*`. The scored labels, ranks, scores, tiling and MHC flags are unchanged; this appended annotated survey and current TSV supersede the first annotation attempt.', '',
        '[observation] **Annotation bug:** a single min-start/max-end span per symbol/chrom/strand bridged disjoint refGene mappings. For example, RNU1-1, RNU1-2, RNU1-3, RNU1-4 and RNVU1-18 acquired an artificial chr1 interval [16840616,149224221), spanning ~132 Mb. The corrected rule is a SET of transcript intervals per symbol/chrom/strand, merging only overlapping or touching intervals (gap tolerance 0 bp). Every overlap must meet an actual transcript interval, not an invented inter-transcript gap.', '',
        f"[observation] Disjoint-interval diagnostics: {correction['symbols_with_disjoint_intervals']} symbols / {correction['disjoint_multi_interval_symbol_chrom_strand_groups']} symbol-chrom-strand groups have >1 merged interval. Spans >4 Mb are marked and excluded from overlap/nearest annotation: {correction['excluded_pathological_interval_count']} intervals; omitted pathological symbols: {', '.join(correction['excluded_pathological_symbols']) or 'none'}. The interval file retains these flags for audit. Pathology/gap tests in `analysis/test_block_survey.py` check a disjoint synthetic interval and actual refGene union coverage; no merged interval may bridge an uncovered source gap.", '',
        f"[observation] Annotation damage: {correction['changed_blocks']}/{len(blocks)} blocks change overlap/nearest annotation; {correction['changed_registered_neither_blocks']} registered neither-labelled blocks change; {correction['changed_previously_reported_novel_blocks']}/10 previously listed novel tiles change. MHC flags unchanged: {correction['mhc_flags_unchanged']}; rank/score unchanged: {correction['ranks_scores_unchanged']}. Detailed lost/gained symbols are in `block_survey_summary.json` annotation_correction; corrected top lists appear in this appended survey.", '',
        '[observation] **Empty moderate band is a panel CONSTRUCTION ARTIFACT:** cases were selected at |X|>5.45 and controls at |X|<3.16. The 3.16–5.45 band was never scored in these shards. Consequently registered scored-only moderate==strong; registered “neither” must NOT be read as no evidence anywhere in Akbari.', '',
        f"[observation] Supplemental POST-REGISTRATION evidence uses the local QC cache ({cache['cache_variants']:,} variants, {cache['cache_strong_variants']:,} strong, {cache['cache_moderate_band_variants']:,} in 3.16<=|X|<5.45; exact |X|=5.45 count {cache['cache_exact_5_45_variants']}). Cache moderate denotes the BAND, not moderate-including-strong. TSV fields: `block_max_absX_cache`, `block_has_strong_cache` (>5.45), `block_has_moderate_cache` (3.16<=|X|<5.45), `block_has_any_akbari_cache` (strong OR band), `n_cache_variants_in_block`; matching `_scored` fields distinguish the original panel. Multiple evidence classes can occur in one block. Null maximum if no cache rows; no-cache flags false means absence of evidence in this cache only.", '',
        '|Frozen top N|Cache strong evidence|Cache moderate-only evidence|No Akbari cache evidence|Registered non-strong novel|Of those: cache moderate evidence|Of those: no cache evidence|',
        '|---:|---:|---:|---:|---:|---:|---:|']
    for r in summary['top_N']:
        lines.append(f"|{r['N']}|{r['cache_strong']}|{r['cache_moderate_only']}|{r['cache_no_evidence']}|{r['non_strong']}|{r['registered_novel_with_cache_moderate']}|{r['registered_novel_with_no_cache_evidence']}|")
    lines += ['',
        '[inference] The genuinely no-Akbari-cache-evidence count is the honest novelty estimate in this available resource, smaller than the registered scored-only novel count wherever moderate evidence appears. This is evidence adjudication, NOT an enrichment analysis of the unscored moderate band: its variants have no measured ranking here. No variant was rescored or block reranked. Novel still does not mean neutral or globally unreported.', '',
        '### Corrected highest-ranked no-cache-evidence hypotheses', '',
        '[observation] Additional post-registration display: first ten frozen-ranked blocks with no strong or moderate evidence in the QC cache. This does not change the registered first-ten scored-neither display above; it distinguishes which of those survive cache adjudication. Gene/trait support remains unadjudicated.', '',
        candidate_table(blocks[~blocks.block_has_any_akbari_cache].head(10)), '']
    if target.exists():
        with target.open('a') as f:
            f.write('\n\n---\n\n## Appended regenerated survey — '+datetime.now(timezone.utc).isoformat()+'\n\n'+'\n'.join(lines))
    else:
        target.write_text('\n'.join(lines))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--fetch-annotation', action='store_true')
    args = parser.parse_args()
    if not Path('docs/block-survey-prereg.md').exists():
        raise FileNotFoundError('Pre-registration must be persisted before the survey')
    metadata = ensure_annotation(args.fetch_annotation)
    genes = parse_refgene(ANNOTATION)
    genes.to_csv(OUT / 'block_survey_gene_spans.tsv', sep='\t', index=False)
    variants = pd.concat([pd.read_csv(path, sep='\t', usecols=['locus', 'chrom', 'pos', 'akbari_X', 'total_score'])
                          for path in SHARDS], ignore_index=True)
    blocks = annotate_blocks(build_blocks(variants), genes)
    loci = None
    try:
        loci = extract_independent_loci()
        loci.to_csv(OUT / 'block_survey_akbari_479_tags.tsv', sep='\t', index=False)
        independent_status = ('Extracted all 479 unique indexed tag coordinates from PDF pages 52–91 (Figures S5.1–S5.40); '
            'all 479 match local GRCh37 cache coordinates and rounded |X| within 0.015. Output: '
            '`data/derived/block_survey_akbari_479_tags.tsv`. The 347 pruned-list identities / LD intervals were not extracted.')
    except (ValueError, OSError, KeyError) as error:
        independent_status = f'Independent-locus list unavailable: {error}; no substitute coordinates guessed.'
    blocks = add_evidence(blocks, loci)
    blocks, cache_summary = cache_evidence(blocks, variants)
    correction = annotation_correction_diagnostics(blocks, genes)
    blocks.to_csv(OUT / 'block_survey.tsv', sep='\t', index=False, na_rep='NA')
    phenotype = phenotype_context()
    (OUT / 'block_survey_table_s5_context.json').write_text(json.dumps(phenotype, indent=2) + '\n')
    summary = summarize(blocks, variants, metadata, independent_status)
    summary['post_registration_cache_evidence'] = cache_summary
    summary['annotation_correction'] = correction
    (OUT / 'block_survey_summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    print('[observation] Incremental table and summary saved:', len(blocks), 'blocks', flush=True)
    write_results(blocks, summary, phenotype)
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    import scripts.provenance as provenance
    provenance.get_git_info = lambda: {'disabled': 'operator prohibits all git commands'}
    recorder = provenance.RunRecorder('block-survey', command=['.venv/bin/python', 'analysis/block_survey.py', *sys.argv[1:]])
    for path in [Path('docs/block-survey-prereg.md'), Path('analysis/block_survey.py'), *SHARDS,
                 ANNOTATION, ANNOTATION_META, PDF, WORKBOOK, CACHE,
                 OUT / 'block_survey_initial_bounding.tsv',
                 Path('analysis/prioritization_score.py'), Path('data/derived/prioritization/metrics.json'),
                 *[Path(f'data/derived/prioritization/simulation_{m}.tsv') for m in ['random_walk', 'random_walk_jump', 'pure_replacement']]]:
        recorder.record_input(path)
    recorder.record_parameter('tile_size_bp', TILE)
    recorder.record_parameter('nearest_gene_limit_bp', NEAREST_LIMIT)
    recorder.record_parameter('block_null_available', False)
    recorder.record_parameter('independent_locus_status', independent_status)
    for path in [*OUT.glob('block_survey*'), Path('docs/block-survey-results.md')]:
        if path.is_file():
            recorder.record_output(path)
    recorder.save()
    print('[observation] Survey complete; no discovery claim. Top-N counts:', json.dumps(summary['top_N']))


if __name__ == '__main__':
    main()
