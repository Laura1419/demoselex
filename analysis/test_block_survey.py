"""[observation] Source-named tests for the exploratory block survey and its artifacts."""
import gzip
import json
import re
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

from analysis import block_survey as survey


def variants(rows):
    return pd.DataFrame(rows, columns=['locus', 'chrom', 'pos', 'akbari_X', 'total_score'])


def gene_frame(rows):
    return pd.DataFrame(rows, columns=['gene_symbol', 'chrom', 'strand', 'start', 'end'])


class BlockSurveyRules(unittest.TestCase):
    def test_tile_boundaries_representative_ties_and_any_variant_X(self):
        frame = variants([('z', 1, 1, 0, 5), ('a', 1, 250000, 3.16, 5),
                          ('b', 1, 250001, 5.45, 4), ('c', 1, 250002, 5.46, 1)])
        result = survey.build_blocks(frame)
        self.assertEqual(result.representative_locus.tolist(), ['a', 'b'])
        self.assertEqual(result.block_start.tolist(), [0, 250000])
        self.assertEqual(result.n_variants_scored_in_block.tolist(), [2, 2])
        self.assertEqual(result.akbari_label.tolist(), ['moderate_only', 'strong'])
        self.assertEqual(result.best_abs_X.tolist(), [3.16, 5.46])

    def test_missing_X_is_unknown_unless_observed_threshold_crossed(self):
        result = survey.build_blocks(variants([('a', 1, 1, np.nan, 5),
            ('b', 1, 2, 0, 4), ('c', 1, 250001, 6, 3), ('d', 1, 250002, np.nan, 2)]))
        self.assertEqual(result.akbari_label.tolist(), ['unknown', 'strong'])
        self.assertTrue(pd.isna(result.iloc[0].is_akbari_strong))
        self.assertTrue(result.iloc[1].is_akbari_strong)

    def test_duplicate_and_nonfinite_scores_fail_loudly(self):
        with self.assertRaisesRegex(ValueError, 'Duplicate'):
            survey.build_blocks(variants([('a', 1, 1, 0, 1), ('a', 1, 2, 0, 0)]))
        with self.assertRaisesRegex(ValueError, 'Nonfinite'):
            survey.build_blocks(variants([('a', 1, 1, 0, -np.inf)]))

    def test_mhc_half_open_boundary(self):
        result = survey.build_blocks(variants([('a', 6, 28500000, 0, 4),
            ('b', 6, 28500001, 0, 3), ('c', 6, 33500000, 0, 2), ('d', 6, 33500001, 0, 1)]))
        self.assertEqual(result.is_mhc.tolist(), [False, True, True, False])

    def test_overlap_nearest_ties_limits_and_explicit_null(self):
        block = survey.build_blocks(variants([('a', 1, 250001, 0, 1)]))
        genes = gene_frame([('left', 1, '+', 249990, 250000),
                            ('right', 1, '+', 500000, 500010),
                            ('inside', 1, '+', 300000, 300010)])
        result = survey.annotate_blocks(block, genes).iloc[0]
        self.assertEqual(result.genes_overlapping, 'inside')
        self.assertIsNone(result.nearest_gene)
        result = survey.annotate_blocks(block, genes.iloc[:2]).iloc[0]
        self.assertEqual(result.nearest_gene, 'left')
        self.assertEqual(result.nearest_gene_distance, 0)
        result = survey.annotate_blocks(block, gene_frame([('edge', 1, '+', 600000, 600010)])).iloc[0]
        self.assertEqual(result.nearest_gene_distance, 100000)
        result = survey.annotate_blocks(block, gene_frame([('outside', 1, '+', 600001, 600010)])).iloc[0]
        self.assertIsNone(result.nearest_gene)
        self.assertEqual(result.annotation_status, 'no_gene_within_100kb')

    def test_disjoint_RNU_pathology_gap_not_fabricated(self):
        raw = gene_frame([('RNU_test', 1, '-', 16840616, 16840780),
                          ('RNU_test', 1, '-', 149224000, 149224221),
                          ('touch', 1, '+', 10, 20), ('touch', 1, '+', 20, 30)])
        merged = survey.merge_transcript_intervals(raw)
        rnu = merged[merged.gene_symbol.eq('RNU_test')]
        self.assertEqual(len(rnu), 2)
        self.assertEqual(merged[merged.gene_symbol.eq('touch')][['start', 'end']].values.tolist(), [[10, 30]])
        result = survey.annotate_blocks(survey.build_blocks(variants([('a', 1, 52000001, 0, 1)])), merged).iloc[0]
        self.assertIsNone(result.genes_overlapping)
        self.assertIsNone(result.nearest_gene)

    def test_pathological_span_marked_and_excluded(self):
        merged = survey.merge_transcript_intervals(gene_frame([('huge', 1, '+', 0, 4000001)]))
        self.assertTrue(merged.iloc[0].excluded_pathological_span)
        result = survey.annotate_blocks(survey.build_blocks(variants([('a', 1, 1, 0, 1)])), merged).iloc[0]
        self.assertIsNone(result.genes_overlapping)
        self.assertIsNone(result.nearest_gene)


class BlockSurveyArtifacts(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.blocks = pd.read_csv('data/derived/block_survey.tsv', sep='\t', float_precision='round_trip')
        cls.genes = pd.read_csv('data/derived/block_survey_gene_spans.tsv', sep='\t')
        cls.summary = json.loads(Path('data/derived/block_survey_summary.json').read_text())

    def test_registered_denominators_labels_and_frozen_ranking(self):
        b = self.blocks
        source = pd.concat([pd.read_csv(path, sep='\t', usecols=['locus', 'chrom', 'pos', 'akbari_X', 'total_score'])
                            for path in survey.SHARDS], ignore_index=True)
        expected = survey.build_blocks(source)
        self.assertEqual(len(source), 181933)
        self.assertEqual(int(source.akbari_X.abs().gt(5.45).sum()), 31939)
        self.assertEqual(len(b), 10615)
        self.assertEqual(int(b.is_akbari_strong.sum()), 723)
        self.assertFalse(b.duplicated(['chrom', 'block_start']).any())
        self.assertEqual(b.representative_locus.tolist(), expected.representative_locus.tolist())
        self.assertEqual(b.akbari_label.tolist(), expected.akbari_label.tolist())
        self.assertTrue(np.allclose(b.total_score, expected.total_score, rtol=0, atol=0))
        self.assertEqual(int(b.n_variants_scored_in_block.sum()), len(source))
        self.assertTrue(self.summary['annotation_correction']['ranks_scores_unchanged'])
        self.assertTrue(self.summary['annotation_correction']['mhc_flags_unchanged'])

    def test_all_gene_overlap_calls_match_actual_raw_transcripts(self):
        raw = survey.read_refgene_transcripts(survey.ANNOTATION)
        pathological = self.genes[self.genes.excluded_pathological_span]
        self.assertEqual(len(pathological), 0, 'If real pathological spans appear, explicitly adapt raw-union audit')
        by_chrom = {c: g for c, g in raw.groupby('chrom')}
        for r in self.blocks.itertuples():
            g = by_chrom[r.chrom]
            actual = sorted(g[g.start.lt(r.block_end) & g.end.gt(r.block_start)].gene_symbol.unique())
            annotated = [] if pd.isna(r.genes_overlapping) else r.genes_overlapping.split(';')
            self.assertEqual(annotated, actual, f'chr{r.chrom}:{r.block_start}')
            if pd.notna(r.nearest_gene):
                self.assertFalse(actual)
                spans = self.genes[(self.genes.chrom == r.chrom) & (self.genes.gene_symbol == r.nearest_gene)]
                gap = np.maximum(np.maximum(spans.start - r.block_end, r.block_start - spans.end), 0).min()
                self.assertEqual(r.nearest_gene_distance, gap)
                self.assertLessEqual(gap, 100000)

    def test_real_merged_intervals_have_no_uncovered_source_gaps(self):
        raw = survey.read_refgene_transcripts(survey.ANNOTATION)
        groups = {key: g.sort_values(['start', 'end']) for key, g in raw.groupby(['gene_symbol', 'chrom', 'strand'])}
        for r in self.genes.itertuples():
            source = groups[(r.gene_symbol, r.chrom, r.strand)]
            inside = source[source.start.lt(r.end) & source.end.gt(r.start)]
            self.assertEqual(int(inside.start.min()), r.start)
            reach = r.start
            for transcript in inside.itertuples():
                self.assertLessEqual(transcript.start, reach, f'gap fabricated for {r.gene_symbol}')
                reach = max(reach, transcript.end)
            self.assertEqual(reach, r.end)

    def test_provenance_gzip_registration_order_and_evidence_nulls(self):
        meta = json.loads(survey.ANNOTATION_META.read_text())
        self.assertEqual(meta['sha256'], survey.sha256(survey.ANNOTATION))
        with gzip.open(survey.ANNOTATION, 'rb') as f:
            self.assertGreater(len(f.read()), 0)
        self.assertGreater(meta['retrieved_at_utc'], '2026-10-09T02:18:36Z')
        prereg = Path('docs/block-survey-prereg.md').read_text()
        self.assertIn('before computing block rankings', prereg)
        for flag in ['gwas_catalog', 'eqtl', 'functional_annotation']:
            self.assertTrue(self.blocks[f'{flag}_evidence_status'].eq('unavailable_locally_not_checked').all())
        self.assertTrue(self.blocks.table_s5_locus_phenotype.isna().all())
        self.assertTrue(self.blocks.akbari_347_locus_overlap.isna().all())
        tags = pd.read_csv('data/derived/block_survey_akbari_479_tags.tsv', sep='\t')
        self.assertEqual(len(tags), 479)
        self.assertEqual(set(tags.trajectory_index), set(range(1, 480)))
        self.assertTrue((tags.cached_X.abs() - tags.pdf_abs_X).abs().lt(.015).all())

    def test_cache_adjudication_counts_and_no_moderate_enrichment(self):
        b = self.blocks
        self.assertEqual(int(b.block_has_moderate_scored.sum()), 0)
        self.assertTrue(b.block_has_any_akbari_cache.equals(b.block_has_strong_cache | b.block_has_moderate_cache))
        self.assertTrue(b.block_has_strong_cache.equals(b.block_has_strong_scored))
        for row in self.summary['top_N']:
            top = b.head(row['N'])
            self.assertEqual(row['cache_strong'] + row['cache_moderate_only'] + row['cache_no_evidence'], row['N'])
            self.assertEqual(row['registered_novel_with_cache_moderate'] + row['registered_novel_with_no_cache_evidence'], row['non_strong'])
            self.assertNotIn('cache_moderate_enrichment', row)
            self.assertEqual(row['cache_no_evidence'], int((~top.block_has_any_akbari_cache).sum()))
        self.assertEqual(self.summary['post_registration_cache_evidence']['cache_exact_5_45_variants'], 0)

    def test_current_results_candidate_rows_gene_traceability_and_FP_numbers(self):
        text = Path('docs/block-survey-results.md').read_text()
        # The most recently appended survey is the corrected authoritative version.
        current = text.split('## Appended regenerated survey —')[-1]
        b = self.blocks.set_index('rank')
        raw_symbols = set(survey.read_refgene_transcripts(survey.ANNOTATION).gene_symbol)
        tables = 0
        for line in current.splitlines():
            if re.match(r'^\|\d+\|chr\d+:', line):
                fields = line.split('|')
                rank, locus = int(fields[1]), fields[3]
                row = b.loc[rank]
                self.assertEqual(locus, row.representative_locus)
                self.assertEqual(fields[8], survey.gene_display(row))
                symbols = [] if pd.isna(row.genes_overlapping) else row.genes_overlapping.split(';')
                if pd.notna(row.nearest_gene):
                    symbols.append(row.nearest_gene)
                self.assertTrue(set(symbols) <= raw_symbols)
                tables += 1
        self.assertEqual(tables, 40)
        self.assertIn('32/300 = 10.67%', current)
        self.assertIn('43.73%', current)
        self.assertIn('POST-REGISTRATION', current)
        self.assertIsNone(self.summary['simulation']['block_level_null_rate'])


if __name__ == '__main__':
    unittest.main()
