# Exploratory block-level gene survey — 2026-10-09

[inference] **An exploratory hypothesis list for downstream modern-phenotype follow-up, NOT a selection discovery.** False positives are expected. Gene proximity and Akbari agreement do not establish function, causality, or independent validation.

## Registration and reproducibility

[observation] Protocol persisted first in `docs/block-survey-prereg.md` at 2026-10-09T02:18:36Z; annotation retrieval followed at 2026-10-09T02:19:41.016138+00:00. Inputs are the four scored shard TSVs, not the earlier small panel. Reproduce with `.venv/bin/python analysis/block_survey.py` (add `--fetch-annotation` only if the raw annotation is absent). Explicit TSV null token: `NA`.

[observation] UCSC hg19 refGene URL: https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/refGene.txt.gz; gzip CRC verified; SHA256 `28e7fe22cf1923a6887c7b60e6004c75987f810b8374d4d396021a34ac0410d1`. Provenance: `data/raw/annotation/hg19_refGene.provenance.json`. Parsed symbol/chromosome/transcript-union span/strand: `data/derived/block_survey_gene_spans.tsv`. All gene names in the candidate/control tables below are generated from annotated rows of `data/derived/block_survey.tsv`, never from PDF gene labels.

[observation] 181,933 finite scored variants, 31,939 strong variants; 10,615 occupied 250 kb tiles, 723 strong tiles and 723 moderate-including-strong tiles. No duplicate loci or missing input X. Annotation statuses: `{'overlap': 9165, 'no_gene_within_100kb': 862, 'nearest_within_100kb': 588}`.

[observation] The task's 181,933 and 31,939 variant counts agree. Approximately 723 strong blocks compares with the measured 723; tiling is explicitly one-based-pos to zero-based-start. Tiles are LD-reduced, not independent: adjacent tiles and MHC can share one signal.

## Top 20 blocks and local evidence

[observation] Strong means any scored variant |X|>5.45; moderate_only means max |X|>=3.16 but not strong. Moderate flags in TSV include strong. Prior panel = task-supplied canonical ancient-DNA support panel, also represented by the existing prioritizer controls; sparse illustrative controls, not an exhaustive literature review. A non-match is not evidence of no prior support. Coordinates/gene spans are zero-based half-open; representative locus strings contain one-based positions.

|Rank|GRCh37 tile [start,end)|Representative|Score|Best abs X|Label|MHC|Genes (overlap unless stated)|Prior panel|479 tag in tile / rep|
|---:|---|---|---:|---:|---|---|---|---|---|
|1|chr6:[31250000,31500000)|6_31359394_C_T|13.419|13.733|strong|yes|ATP6V1G2-DDX39B;DDX39B;HCG26;HCP5;HLA-B;LINC01149;LINC02571;MCCD1;MICA;MICA-AS1;MICB;MICB-DT;MIR6891|no fixed-panel match|unknown / unknown|
|2|chr6:[32500000,32750000)|6_32627830_C_T|12.900|13.704|strong|yes|HLA-DQA1;HLA-DQA2;HLA-DQB1;HLA-DQB1-AS1;HLA-DQB2;HLA-DRB1;HLA-DRB6;MIR3135B|no fixed-panel match|unknown / unknown|
|3|chr6:[28500000,28750000)|6_28689672_C_A|11.405|10.190|strong|yes|GPX5;LINC00533;ZBED9|no fixed-panel match|unknown / unknown|
|4|chr6:[30750000,31000000)|6_30789647_G_C|11.125|10.596|strong|yes|DDR1;GTF2H4;HCG20;HCG21;LINC00243;LINC02570;MIR4640;MUC21;MUC22;MUCL3;SFTA2;VARS2|no fixed-panel match|unknown / unknown|
|5|chr6:[31500000,31750000)|6_31577099_A_T|11.010|11.446|strong|yes|ABHD16A;AIF1;APOM;ATP6V1G2;ATP6V1G2-DDX39B;BAG6;C6orf47;CLIC1;CSNK2B;DDAH2;DDX39B;DDX39B-AS1;GPANK1;LOC100287329;LST1;LTA;LTB;LY6G5B;LY6G5C;LY6G6C;LY6G6D;LY6G6E;LY6G6F;LY6G6F-LY6G6D;MIR4646;MIR6832;MPIG6B;MSH5;MSH5-SAPCD1;NCR3;NFKBIL1;PRRC2A;SAPCD1;SAPCD1-AS1;SNORA38;SNORD117;SNORD84;TNF;VARS1;VWA7|no fixed-panel match|unknown / unknown|
|6|chr6:[32000000,32250000)|6_32158319_A_G|10.304|12.460|strong|yes|AGER;AGPAT1;ATF6B;C4B;C4B_2;CYP21A2;EGFL8;FKBPL;GPSM3;LOC100507547;MIR6721;MIR6833;NOTCH4;PBX2;PPT2;PPT2-EGFL8;PRRT1;RNF5;TNXB;TSBP1-AS1|no fixed-panel match|unknown / unknown|
|7|chr6:[32750000,33000000)|6_32778571_A_ATT|10.078|9.699|strong|yes|BRD2;HLA-DMA;HLA-DMB;HLA-DOA;HLA-DOB;LOC100294145;PSMB8;PSMB8-AS1;PSMB9;TAP1;TAP2|no fixed-panel match|unknown / unknown|
|8|chr6:[31750000,32000000)|6_31757881_TA_T|9.797|13.835|strong|yes|C2;C2-AS1;C4A;C4B;C4B_2;CFB;CYP21A1P;DXO;EHMT2;HSPA1A;HSPA1B;HSPA1L;LOC110384692;LSM2;MIR1236;NELFE;NEU1;SKIV2L;SLC44A4;SNHG32;SNORD48;SNORD52;STK19;TNXA;TNXB;VARS1;ZBTB12|no fixed-panel match|unknown / unknown|
|9|chr6:[29000000,29250000)|6_29144532_C_T|9.733|8.248|strong|yes|LOC100129636;LOC101929006;OR2B3;OR2J1;OR2J2;OR2J3;OR2W1|no fixed-panel match|unknown / unknown|
|10|chr2:[136500000,136750000)|2_136608646_G_A|9.682|16.523|strong|no|DARS-AS1;DARS1;LCT;LCT-AS1;LOC107985946;MCM6;UBXN4|LCT;MCM6|unknown / unknown|
|11|chr2:[136000000,136250000)|2_136161092_C_G|9.450|15.257|strong|no|ZRANB3|no fixed-panel match|unknown / unknown|
|12|chr6:[32250000,32500000)|6_32343369_T_C|9.317|11.222|strong|yes|BTNL2;HCG23;HLA-DRA;HLA-DRB5;TSBP1;TSBP1-AS1|no fixed-panel match|unknown / unknown|
|13|chr2:[135250000,135500000)|2_135498924_G_C|9.269|14.895|strong|no|TMEM163|no fixed-panel match|unknown / unknown|
|14|chr6:[26250000,26500000)|6_26338128_T_C|8.832|8.353|strong|no|BTN2A1;BTN2A2;BTN2A3P;BTN3A1;BTN3A2;BTN3A3;H2BC10;H2BC9;H3C7;H3C8;H4C8;LOC285819|no fixed-panel match|unknown / unknown|
|15|chr6:[30250000,30500000)|6_30345118_C_T|8.734|11.967|strong|yes|HCG17;HCG18;HLA-E;LINC02569;RPP21;TRIM39;TRIM39-RPP21|no fixed-panel match|unknown / unknown|
|16|chr2:[136250000,136500000)|2_136381348_G_A|8.682|15.815|strong|no|LOC107985946;MIR128-1;R3HDM1;UBXN4;ZRANB3|no fixed-panel match|unknown / unknown|
|17|chr6:[31000000,31250000)|6_31196862_G_A|8.649|12.082|strong|yes|C6orf15;CCHCR1;CDSN;HCG22;HCG27;HLA-C;MUC22;POU5F1;PSORS1C1;PSORS1C2;PSORS1C3;TCF19|no fixed-panel match|unknown / unknown|
|18|chr6:[33000000,33250000)|6_33029081_CCT_C|8.504|11.391|strong|yes|B3GALT4;COL11A2;HCG24;HCG25;HLA-DPA1;HLA-DPB1;HLA-DPB2;HSD17B8;MIR219A1;RING1;RPS18;RXRB;SLC39A7;VPS52;WDR46|no fixed-panel match|unknown / unknown|
|19|chr2:[135500000,135750000)|2_135689188_G_T|8.474|14.365|strong|no|ACMSD;CCNT2;CCNT2-AS1;MAP3K19;MIR5590|no fixed-panel match|unknown / unknown|
|20|chr6:[29500000,29750000)|6_29542947_T_C|8.332|10.709|strong|yes|GABBR1;HLA-F;HLA-F-AS1;IFITM4P;LINC01015;MOG;OR2H2;SNORD32B;UBD;ZFP57|no fixed-panel match|unknown / unknown|

[inference] Overlapping genes are positional candidates only (full transcript bounding spans, not coding-exon intersection). Nearest fallback is allowed only when no gene overlaps, at <=100 kb. MHC flag = overlap chr6:[28500000,33500000). Multiple high-ranked MHC tiles remain one extended LD complex, not multiple independent biological confirmations.

## Prior ancient-DNA support panel: what surfaced?

|Control group|Best annotated block rank|Block representative|Annotated panel genes in that row|Label|
|---|---:|---|---|---|
|LCT/MCM6|10|2_136608646_G_A|LCT;MCM6|strong|
|SLC24A5|10510|15_48269895_T_G|SLC24A5|strong|
|SLC45A2|49|5_33951693_C_G|SLC45A2|strong|
|TYR|99|11_88911696_C_A|TYR|strong|
|OCA2/HERC2|2995|15_28225348_G_A|OCA2|neither|
|FADS1/FADS2|89|11_61618012_C_T|FADS1;FADS2|strong|
|TLR1/TLR6|34|4_38798648_C_A|TLR1;TLR6|strong|
|AMY1/AMY2 region|6814|1_104499856_C_T|AMY1A;AMY1B;AMY1C|neither|

[observation] Pre-registered predictions: top-100 strong enrichment is 12.63× (prediction >1: met); 5 top-100 tiles have a fixed-panel gene (prediction >=1: met). Neither success supplies independent validation.

[observation] Inherited variant-level limits remain: SLC24A5 ranked approximately 9,161 and the AMY1 SNP proxy approximately 698 in the earlier 10,641-variant panel. The table above reports NEW block ranks in the much denser merged survey and does not replace those prior variant ranks. The within-lineage residual design is structurally blind to carriage-driven change predominantly BETWEEN demes; it should not be presented as capturing those mechanisms. The amylase region SNP proxy is not a copy-number selection test.

## Akbari reported independent loci and phenotype evidence

[observation] Independent-locus list unavailable: Merge keys are not unique in right dataset; not a one-to-one merge; no substitute coordinates guessed.

[inference] “479 tag in tile” checks whether an extracted tag coordinate lies in the tile, and “rep” checks exact representative locus ID against its verified cached tag ID. It is NOT overlap with an LD-defined independent-locus interval. Such intervals and membership of the more strictly pruned 347 set are not available from this extraction; `akbari_347_locus_overlap=NA` throughout. Absence of a tag from a tile does not imply independence or absence of the reported signal. Section 5 describes 479 tags (410 outside HLA +69 inside); Table S5.1 supplies D-prime pruning COUNTS, not the identities of a 347-tag subset. This is a difference in list/definition, not a claim that 347 is incorrect.

[observation] Table S5 sheet `STable5(696.polygenic.tests)` contains 696 aggregate phenotype tests, with beta/sign-beta analyses with and without HLA. It contains no locus or gene mapping. The strongest five absolute z values in the beta-without-HLA column are shown solely to describe available context, not to connect these phenotypes to the ranked tiles:

|Phenotype ID|Phenotype / coding description|z beta without HLA|z beta with HLA|
|---|---|---:|---:|
|Pan-UKBB.4004|Skin colour|-18.195|-18.249|
|Pan-UKBB.137|Hair colour (natural, before greying) / Black|-15.784|-15.783|
|UKB_460K.pigment_SUNBURN|Sunburn|14.070|14.023|
|Pan-UKBB.136|Hair colour (natural, before greying) / Dark brown|-13.005|-12.968|
|Pan-UKBB.4453|Average total household income before tax|9.462|9.514|

[inference] Aggregate Table S5 results cannot adjudicate any candidate gene; all locus-phenotype cells are explicitly null. **GWAS Catalog, eQTL, and functional annotation are NOT available locally and were not checked.** Their evidence statuses are explicit in every TSV row. RefGene transcript spans are positional annotations, not functional evidence.

## False-positive / follow-up cost at block level

[observation] Base prevalence in the 10,615 scored tiles: strong 723/10,615 = 6.81%; moderate including strong 723/10,615 = 6.81%. Denominator is the ascertained scored panel, not all genomic tiles.

|Top N|Strong (count / N)|Strong enrichment|Moderate incl. strong (count / N)|Moderate enrichment|Moderate only|Non-strong|Neither (strict novel)|MHC|479 tag tiles|
|---:|---|---:|---|---:|---:|---:|---:|---:|---:|
|100|86/100 (86.0%)|12.63×|86/100 (86.0%)|12.63×|0|14|14|19|0|
|500|188/500 (37.6%)|5.52×|188/500 (37.6%)|5.52×|0|312|312|20|0|
|1000|244/1000 (24.4%)|3.58×|244/1000 (24.4%)|3.58×|0|756|756|20|0|

[inference] Counts not strong are novel relative to the strong threshold; counts neither are the stricter novel subset, used below. These are follow-up costs, NOT proven false positives: non-significance is not neutrality. Strong and moderate agreement reuses Akbari data (same imputed genotypes), so enrichment is NOT independent validation; original variant ascertainment and tile maxima further favor strong/dense regions.

[observation] A calibrated LD-honest **block drift-null FP rate could not be computed**. Existing `analysis/prioritization_score.py::drift_simulation` / `--simulations` generates independent simulated single-locus genotype columns with no genomic coordinates or linked haplotypes. Arbitrary grouping would invent an LD model. The saved 300-draw drift+jump tail above the EARLIER top-100 variant threshold is 32/300 = 10.67% (agrees with quoted 10.67%). That legacy run scored 10,641 variants, NOT 181,933; all four merged shard runs used zero simulation replicates. Therefore the quoted value is not a merged-survey or block FDR.

[observation] Reusing the existing simulation draws, the following SINGLE-LOCUS tails exceed each observed BLOCK-rank score threshold. This is only a threshold diagnostic; it omits within-tile multiplicity and LD, so it is NOT the requested calibrated block false-positive probability or list precision:

|Top-block N|Block score threshold|Brownian drift only exceedances / 300|Drift + jump exceedances / 300|Replacement-only exceedances / 300|
|---:|---:|---|---|---|
|100|1.9056|0/300 (0.00%)|6/300 (2.00%)|0/300 (0.00%)|
|500|0.0863|0/300 (0.00%)|35/300 (11.67%)|19/300 (6.33%)|
|1000|-0.4146|1/300 (0.33%)|59/300 (19.67%)|43/300 (14.33%)|

[inference] Finite simulation size matters: even 0/300 has a one-sided binomial 95% upper tail probability about 0.994%, not zero. No observed-block FP count is inferable from these single-locus rates. A linked genotype/path simulation with observed block sizes and missingness would be needed to estimate block maxima under drift.

[observation] Earlier within-node bootstrap top-100 variant membership retention = 43.73% across 30 resamples (agrees with approximately 43.7%). This has not been recomputed at block level or on the merged panel; it warns that ranks near cuts are unstable, not that block retention equals that number.

## Highest-ranked novel blocks: bounded follow-up hypotheses

[observation] The first ten tiles with neither strong nor moderate scored-variant X are listed mechanically (not selected by appealing gene stories). “Interesting” here means rank/proximity under the pre-registered exploratory protocol only. The broader non-strong candidate set additionally includes moderate-only tiles.

|Rank|GRCh37 tile [start,end)|Representative|Score|Best abs X|Label|MHC|Genes (overlap unless stated)|Prior panel|479 tag in tile / rep|
|---:|---|---|---:|---:|---|---|---|---|---|
|58|chr1:[52000000,52250000)|1_52025112_T_TACAC|3.746|3.153|neither|no|CALR4P;OSBPL9;RNU1-1;RNU1-2;RNU1-3;RNU1-4;RNVU1-18|no fixed-panel match|unknown / unknown|
|60|chr1:[31000000,31250000)|1_31033851_A_G|3.580|2.930|neither|no|LAPTM5;MATN1;MATN1-AS1;MIR4420;RNU1-1;RNU1-2;RNU1-3;RNU1-4;RNVU1-18|no fixed-panel match|unknown / unknown|
|68|chr16:[4750000,5000000)|16_4958961_A_C|2.950|2.538|neither|no|ANKS3;C16orf71;GLYR1;MIR6511B1;MIR6511B2;PPL;ROGDI;SEPTIN12;SMIM22;UBN1;ZNF500|no fixed-panel match|unknown / unknown|
|69|chr12:[114750000,115000000)|12_114917493_C_T|2.947|1.784|neither|no|TBX5;TBX5-AS1|no fixed-panel match|unknown / unknown|
|76|chr7:[155750000,156000000)|7_155996962_C_T|2.712|2.574|neither|no|LOC389602|no fixed-panel match|unknown / unknown|
|79|chr3:[138500000,138750000)|3_138651944_G_A|2.468|2.596|neither|no|FOXL2;FOXL2NB;LINC01391;PIK3CB;PRR23A;PRR23B|no fixed-panel match|unknown / unknown|
|81|chr4:[1000000,1250000)|4_1213988_C_CAA|2.444|2.217|neither|no|CTBP1;CTBP1-AS;CTBP1-DT;FGFRL1;LOC100130872;LOC105374344;RNF212;SPON2;TMED11P|no fixed-panel match|unknown / unknown|
|82|chr5:[161000000,161250000)|5_161229995_G_A|2.426|2.959|neither|no|GABRA6|no fixed-panel match|unknown / unknown|
|85|chr9:[105000000,105250000)|9_105230080_C_A|2.175|2.503|neither|no|LINC00587 (nearest, 31918 bp)|no fixed-panel match|unknown / unknown|
|87|chr10:[55500000,55750000)|10_55669459_G_T|2.153|2.908|neither|no|PCDH15|no fixed-panel match|unknown / unknown|

[speculation] These tiles could motivate a small modern-phenotype follow-up pilot. Their ranks may instead reflect drift/jitter, carriage residuals, LD with a neighboring signal, or dense-panel maxima; **the local evidence base cannot adjudicate them**. Novel means unlabeled in this scored panel, not globally unreported selection. No causal gene or trait association is claimed.

## Validation and bounded recommendation

[observation] `analysis/test_block_survey.py` covers tile-boundary/tie rules, any-variant labels and missing-X unknowns, MHC boundaries, overlap/nearest distance limits, and output provenance/gene traceability. Run `.venv/bin/python -m unittest analysis.test_block_survey -v` after producing the outputs. Run manifest hashes the protocol, script, all shard tables, raw annotation/provenance, PDF, workbook, cache and saved simulations; git metadata intentionally disabled because the operator prohibits all git commands.

[inference] Use this table as an exploratory, LD-reduced queue, not a selection result or validated production discovery ranking. Budget explicitly for the non-Akbari counts above, cross-tile LD and the uncalibrated block drift cost. First downstream step is independent modern-phenotype evidence acquisition and a linked block-null/stability audit before broad follow-up. Earlier authors' conclusions are untouched.


---

## Appended regenerated survey — 2026-10-09T02:33:18.757156+00:00

# Exploratory block-level gene survey — 2026-10-09

[inference] **An exploratory hypothesis list for downstream modern-phenotype follow-up, NOT a selection discovery.** False positives are expected. Gene proximity and Akbari agreement do not establish function, causality, or independent validation.

## Registration and reproducibility

[observation] Protocol persisted first in `docs/block-survey-prereg.md` at 2026-10-09T02:18:36Z; annotation retrieval followed at 2026-10-09T02:19:41.016138+00:00. Inputs are the four scored shard TSVs, not the earlier small panel. Reproduce with `.venv/bin/python analysis/block_survey.py` (add `--fetch-annotation` only if the raw annotation is absent). Explicit TSV null token: `NA`.

[observation] UCSC hg19 refGene URL: https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/refGene.txt.gz; gzip CRC verified; SHA256 `28e7fe22cf1923a6887c7b60e6004c75987f810b8374d4d396021a34ac0410d1`. Provenance: `data/raw/annotation/hg19_refGene.provenance.json`. Parsed symbol/chromosome/interval-union/strand: `data/derived/block_survey_gene_spans.tsv`. All gene names in the candidate/control tables below are generated from annotated rows of `data/derived/block_survey.tsv`, never from PDF gene labels.

[observation] 181,933 finite scored variants, 31,939 strong variants; 10,615 occupied 250 kb tiles, 723 strong tiles and 723 moderate-including-strong tiles. No duplicate loci or missing input X. Annotation statuses: `{'overlap': 8966, 'no_gene_within_100kb': 979, 'nearest_within_100kb': 670}`.

[observation] The task's 181,933 and 31,939 variant counts agree. Approximately 723 strong blocks compares with the measured 723; tiling is explicitly one-based-pos to zero-based-start. Tiles are LD-reduced, not independent: adjacent tiles and MHC can share one signal.

## Top 20 blocks and local evidence

[observation] Strong means any scored variant |X|>5.45; moderate_only means max |X|>=3.16 but not strong. Moderate flags in TSV include strong. Prior panel = task-supplied canonical ancient-DNA support panel, also represented by the existing prioritizer controls; sparse illustrative controls, not an exhaustive literature review. A non-match is not evidence of no prior support. Coordinates/gene spans are zero-based half-open; representative locus strings contain one-based positions.

|Rank|GRCh37 tile [start,end)|Representative|Score|Best abs X|Scored label|MHC|Genes (overlap unless stated)|Prior panel|479 tag in tile / rep|Cache evidence (POST-REG)|
|---:|---|---|---:|---:|---|---|---|---|---|---|
|1|chr6:[31250000,31500000)|6_31359394_C_T|13.419|13.733|strong|yes|ATP6V1G2-DDX39B;DDX39B;HCG26;HCP5;HLA-B;LINC01149;LINC02571;MCCD1;MICA;MICA-AS1;MICB;MICB-DT;MIR6891|no fixed-panel match|yes / no|strong|
|2|chr6:[32500000,32750000)|6_32627830_C_T|12.900|13.704|strong|yes|HLA-DQA1;HLA-DQA2;HLA-DQB1;HLA-DQB1-AS1;HLA-DQB2;HLA-DRB1;HLA-DRB6;MIR3135B|no fixed-panel match|yes / no|strong|
|3|chr6:[28500000,28750000)|6_28689672_C_A|11.405|10.190|strong|yes|GPX5;LINC00533;ZBED9|no fixed-panel match|no / no|strong|
|4|chr6:[30750000,31000000)|6_30789647_G_C|11.125|10.596|strong|yes|DDR1;GTF2H4;HCG20;HCG21;LINC00243;LINC02570;MIR4640;MUC21;MUC22;MUCL3;SFTA2;VARS2|no fixed-panel match|yes / no|strong|
|5|chr6:[31500000,31750000)|6_31577099_A_T|11.010|11.446|strong|yes|ABHD16A;AIF1;APOM;ATP6V1G2;ATP6V1G2-DDX39B;BAG6;C6orf47;CLIC1;CSNK2B;DDAH2;DDX39B;DDX39B-AS1;GPANK1;LOC100287329;LST1;LTA;LTB;LY6G5B;LY6G5C;LY6G6C;LY6G6D;LY6G6E;LY6G6F;LY6G6F-LY6G6D;MIR4646;MIR6832;MPIG6B;MSH5;MSH5-SAPCD1;NCR3;NFKBIL1;PRRC2A;SAPCD1;SAPCD1-AS1;SNORA38;SNORD117;SNORD84;TNF;VARS1;VWA7|no fixed-panel match|yes / no|strong|
|6|chr6:[32000000,32250000)|6_32158319_A_G|10.304|12.460|strong|yes|AGER;AGPAT1;ATF6B;C4B;C4B_2;CYP21A2;EGFL8;FKBPL;GPSM3;LOC100507547;MIR6721;MIR6833;NOTCH4;PBX2;PPT2;PPT2-EGFL8;PRRT1;RNF5;TNXB;TSBP1-AS1|no fixed-panel match|yes / no|strong|
|7|chr6:[32750000,33000000)|6_32778571_A_ATT|10.078|9.699|strong|yes|BRD2;HLA-DMA;HLA-DMB;HLA-DOA;HLA-DOB;LOC100294145;PSMB8;PSMB8-AS1;PSMB9;TAP1;TAP2|no fixed-panel match|yes / no|strong|
|8|chr6:[31750000,32000000)|6_31757881_TA_T|9.797|13.835|strong|yes|C2;C2-AS1;C4A;C4B;C4B_2;CFB;CYP21A1P;DXO;EHMT2;HSPA1A;HSPA1B;HSPA1L;LOC110384692;LSM2;MIR1236;NELFE;NEU1;SKIV2L;SLC44A4;SNHG32;SNORD48;SNORD52;STK19;TNXA;TNXB;VARS1;ZBTB12|no fixed-panel match|yes / no|strong|
|9|chr6:[29000000,29250000)|6_29144532_C_T|9.733|8.248|strong|yes|LOC100129636;LOC101929006;OR2B3;OR2J1;OR2J2;OR2J3;OR2W1|no fixed-panel match|yes / no|strong|
|10|chr2:[136500000,136750000)|2_136608646_G_A|9.682|16.523|strong|no|DARS-AS1;DARS1;LCT;LCT-AS1;LOC107985946;MCM6;UBXN4|LCT;MCM6|yes / yes|strong|
|11|chr2:[136000000,136250000)|2_136161092_C_G|9.450|15.257|strong|no|ZRANB3|no fixed-panel match|no / no|strong|
|12|chr6:[32250000,32500000)|6_32343369_T_C|9.317|11.222|strong|yes|BTNL2;HCG23;HLA-DRA;HLA-DRB5;TSBP1;TSBP1-AS1|no fixed-panel match|yes / no|strong|
|13|chr2:[135250000,135500000)|2_135498924_G_C|9.269|14.895|strong|no|TMEM163|no fixed-panel match|no / no|strong|
|14|chr6:[26250000,26500000)|6_26338128_T_C|8.832|8.353|strong|no|BTN2A1;BTN2A2;BTN2A3P;BTN3A1;BTN3A2;BTN3A3;H2BC10;H2BC9;H3C7;H3C8;H4C8;LOC285819|no fixed-panel match|yes / no|strong|
|15|chr6:[30250000,30500000)|6_30345118_C_T|8.734|11.967|strong|yes|HCG17;HCG18;HLA-E;LINC02569;RPP21;TRIM39;TRIM39-RPP21|no fixed-panel match|yes / no|strong|
|16|chr2:[136250000,136500000)|2_136381348_G_A|8.682|15.815|strong|no|LOC107985946;MIR128-1;R3HDM1;UBXN4;ZRANB3|no fixed-panel match|yes / no|strong|
|17|chr6:[31000000,31250000)|6_31196862_G_A|8.649|12.082|strong|yes|C6orf15;CCHCR1;CDSN;HCG22;HCG27;HLA-C;MUC22;POU5F1;PSORS1C1;PSORS1C2;PSORS1C3;TCF19|no fixed-panel match|yes / no|strong|
|18|chr6:[33000000,33250000)|6_33029081_CCT_C|8.504|11.391|strong|yes|B3GALT4;COL11A2;HCG24;HCG25;HLA-DPA1;HLA-DPB1;HLA-DPB2;HSD17B8;MIR219A1;RING1;RPS18;RXRB;SLC39A7;VPS52;WDR46|no fixed-panel match|yes / no|strong|
|19|chr2:[135500000,135750000)|2_135689188_G_T|8.474|14.365|strong|no|ACMSD;CCNT2;CCNT2-AS1;MAP3K19;MIR5590|no fixed-panel match|no / no|strong|
|20|chr6:[29500000,29750000)|6_29542947_T_C|8.332|10.709|strong|yes|GABBR1;HLA-F;HLA-F-AS1;IFITM4P;LINC01015;MOG;OR2H2;SNORD32B;UBD;ZFP57|no fixed-panel match|yes / no|strong|

[inference] Overlapping genes are positional candidates only (union of actual transcript spans, not coding-exon intersection). Nearest fallback is allowed only when no gene overlaps, at <=100 kb. Nearest context is NOT positive evidence for a block. MHC flag = overlap chr6:[28500000,33500000). Multiple high-ranked MHC tiles remain one extended LD complex, not multiple independent biological confirmations.

## Prior ancient-DNA support panel: what surfaced?

|Control group|Best annotated block rank|Block representative|Annotated panel genes in that row|Label|
|---|---:|---|---|---|
|LCT/MCM6|10|2_136608646_G_A|LCT;MCM6|strong|
|SLC24A5|10510|15_48269895_T_G|SLC24A5|strong|
|SLC45A2|49|5_33951693_C_G|SLC45A2|strong|
|TYR|99|11_88911696_C_A|TYR|strong|
|OCA2/HERC2|2995|15_28225348_G_A|OCA2|neither|
|FADS1/FADS2|89|11_61618012_C_T|FADS1;FADS2|strong|
|TLR1/TLR6|34|4_38798648_C_A|TLR1;TLR6|strong|
|AMY1/AMY2 region|6814|1_104499856_C_T|AMY1A;AMY1B;AMY1C|neither|

[observation] Pre-registered predictions: top-100 strong enrichment is 12.63× (prediction >1: met); 5 top-100 tiles have a fixed-panel gene (prediction >=1: met). Neither success supplies independent validation.

[observation] Inherited variant-level limits remain: SLC24A5 ranked approximately 9,161 and the AMY1 SNP proxy approximately 698 in the earlier 10,641-variant panel. The table above reports NEW block ranks in the much denser merged survey and does not replace those prior variant ranks. The within-lineage residual design is structurally blind to carriage-driven change predominantly BETWEEN demes; it should not be presented as capturing those mechanisms. The amylase region SNP proxy is not a copy-number selection test.

## Akbari reported independent loci and phenotype evidence

[observation] Extracted all 479 unique indexed tag coordinates from PDF pages 52–91 (Figures S5.1–S5.40); all 479 match local GRCh37 cache coordinates and rounded |X| within 0.015. Output: `data/derived/block_survey_akbari_479_tags.tsv`. The 347 pruned-list identities / LD intervals were not extracted.

[inference] “479 tag in tile” checks whether an extracted tag coordinate lies in the tile, and “rep” checks exact representative locus ID against its verified cached tag ID. It is NOT overlap with an LD-defined independent-locus interval. Such intervals and membership of the more strictly pruned 347 set are not available from this extraction; `akbari_347_locus_overlap=NA` throughout. Absence of a tag from a tile does not imply independence or absence of the reported signal. Section 5 describes 479 tags (410 outside HLA +69 inside); Table S5.1 supplies D-prime pruning COUNTS, not the identities of a 347-tag subset. This is a difference in list/definition, not a claim that 347 is incorrect.

[observation] Table S5 sheet `STable5(696.polygenic.tests)` contains 696 aggregate phenotype tests, with beta/sign-beta analyses with and without HLA. It contains no locus or gene mapping. The strongest five absolute z values in the beta-without-HLA column are shown solely to describe available context, not to connect these phenotypes to the ranked tiles:

|Phenotype ID|Phenotype / coding description|z beta without HLA|z beta with HLA|
|---|---|---:|---:|
|Pan-UKBB.4004|Skin colour|-18.195|-18.249|
|Pan-UKBB.137|Hair colour (natural, before greying) / Black|-15.784|-15.783|
|UKB_460K.pigment_SUNBURN|Sunburn|14.070|14.023|
|Pan-UKBB.136|Hair colour (natural, before greying) / Dark brown|-13.005|-12.968|
|Pan-UKBB.4453|Average total household income before tax|9.462|9.514|

[inference] Aggregate Table S5 results cannot adjudicate any candidate gene; all locus-phenotype cells are explicitly null. **GWAS Catalog, eQTL, and functional annotation are NOT available locally and were not checked.** Their evidence statuses are explicit in every TSV row. RefGene transcript spans are positional annotations, not functional evidence.

## False-positive / follow-up cost at block level

[observation] Base prevalence in the 10,615 scored tiles: strong 723/10,615 = 6.81%; moderate including strong 723/10,615 = 6.81%. Denominator is the ascertained scored panel, not all genomic tiles.

|Top N|Strong (count / N)|Strong enrichment|Moderate incl. strong (count / N)|Moderate enrichment|Moderate only|Non-strong|Neither (strict novel)|MHC|479 tag tiles|
|---:|---|---:|---|---:|---:|---:|---:|---:|---:|
|100|86/100 (86.0%)|12.63×|86/100 (86.0%)|12.63×|0|14|14|19|54|
|500|188/500 (37.6%)|5.52×|188/500 (37.6%)|5.52×|0|312|312|20|118|
|1000|244/1000 (24.4%)|3.58×|244/1000 (24.4%)|3.58×|0|756|756|20|148|

[inference] Counts not strong are novel relative to the strong threshold; counts neither are the stricter novel subset, used below. These are follow-up costs, NOT proven false positives: non-significance is not neutrality. Strong and moderate agreement reuses Akbari data (same imputed genotypes), so enrichment is NOT independent validation; original variant ascertainment and tile maxima further favor strong/dense regions.

[observation] A calibrated LD-honest **block drift-null FP rate could not be computed**. Existing `analysis/prioritization_score.py::drift_simulation` / `--simulations` generates independent simulated single-locus genotype columns with no genomic coordinates or linked haplotypes. Arbitrary grouping would invent an LD model. The saved 300-draw drift+jump tail above the EARLIER top-100 variant threshold is 32/300 = 10.67% (agrees with quoted 10.67%). That legacy run scored 10,641 variants, NOT 181,933; all four merged shard runs used zero simulation replicates. Therefore the quoted value is not a merged-survey or block FDR.

[observation] Reusing the existing simulation draws, the following SINGLE-LOCUS tails exceed each observed BLOCK-rank score threshold. This is only a threshold diagnostic; it omits within-tile multiplicity and LD, so it is NOT the requested calibrated block false-positive probability or list precision:

|Top-block N|Block score threshold|Brownian drift only exceedances / 300|Drift + jump exceedances / 300|Replacement-only exceedances / 300|
|---:|---:|---|---|---|
|100|1.9056|0/300 (0.00%)|6/300 (2.00%)|0/300 (0.00%)|
|500|0.0863|0/300 (0.00%)|35/300 (11.67%)|19/300 (6.33%)|
|1000|-0.4146|1/300 (0.33%)|59/300 (19.67%)|43/300 (14.33%)|

[inference] Finite simulation size matters: even 0/300 has a one-sided binomial 95% upper tail probability about 0.994%, not zero. No observed-block FP count is inferable from these single-locus rates. A linked genotype/path simulation with observed block sizes and missingness would be needed to estimate block maxima under drift.

[observation] Earlier within-node bootstrap top-100 variant membership retention = 43.73% across 30 resamples (agrees with approximately 43.7%). This has not been recomputed at block level or on the merged panel; it warns that ranks near cuts are unstable, not that block retention equals that number.

## Highest-ranked novel blocks: bounded follow-up hypotheses

[observation] The first ten tiles with neither strong nor moderate scored-variant X are listed mechanically (not selected by appealing gene stories). “Interesting” here means rank/proximity under the pre-registered exploratory protocol only. The broader non-strong candidate set additionally includes moderate-only tiles.

|Rank|GRCh37 tile [start,end)|Representative|Score|Best abs X|Scored label|MHC|Genes (overlap unless stated)|Prior panel|479 tag in tile / rep|Cache evidence (POST-REG)|
|---:|---|---|---:|---:|---|---|---|---|---|---|
|58|chr1:[52000000,52250000)|1_52025112_T_TACAC|3.746|3.153|neither|no|CALR4P;OSBPL9|no fixed-panel match|no / no|moderate_only|
|60|chr1:[31000000,31250000)|1_31033851_A_G|3.580|2.930|neither|no|LAPTM5;MATN1;MATN1-AS1;MIR4420|no fixed-panel match|no / no|moderate_only|
|68|chr16:[4750000,5000000)|16_4958961_A_C|2.950|2.538|neither|no|ANKS3;C16orf71;GLYR1;PPL;ROGDI;SEPTIN12;SMIM22;UBN1;ZNF500|no fixed-panel match|no / no|moderate_only|
|69|chr12:[114750000,115000000)|12_114917493_C_T|2.947|1.784|neither|no|TBX5;TBX5-AS1|no fixed-panel match|no / no|moderate_only|
|76|chr7:[155750000,156000000)|7_155996962_C_T|2.712|2.574|neither|no|LOC389602|no fixed-panel match|no / no|moderate_only|
|79|chr3:[138500000,138750000)|3_138651944_G_A|2.468|2.596|neither|no|FOXL2;FOXL2NB;LINC01391;PIK3CB;PRR23A;PRR23B|no fixed-panel match|no / no|moderate_only|
|81|chr4:[1000000,1250000)|4_1213988_C_CAA|2.444|2.217|neither|no|CTBP1;CTBP1-AS;CTBP1-DT;FGFRL1;LOC100130872;LOC105374344;RNF212;SPON2;TMED11P|no fixed-panel match|no / no|moderate_only|
|82|chr5:[161000000,161250000)|5_161229995_G_A|2.426|2.959|neither|no|GABRA6|no fixed-panel match|no / no|moderate_only|
|85|chr9:[105000000,105250000)|9_105230080_C_A|2.175|2.503|neither|no|LINC00587 (nearest, 31918 bp)|no fixed-panel match|no / no|moderate_only|
|87|chr10:[55500000,55750000)|10_55669459_G_T|2.153|2.908|neither|no|PCDH15|no fixed-panel match|no / no|moderate_only|

[speculation] These tiles could motivate a small modern-phenotype follow-up pilot. Their ranks may instead reflect drift/jitter, carriage residuals, LD with a neighboring signal, or dense-panel maxima; **the local evidence base cannot adjudicate them**. Novel means unlabeled in this scored panel, not globally unreported selection. No causal gene or trait association is claimed.

## Validation and bounded recommendation

[observation] `analysis/test_block_survey.py` covers tile-boundary/tie rules, any-variant labels and missing-X unknowns, MHC boundaries, overlap/nearest distance limits, and output provenance/gene traceability. Run `.venv/bin/python -m unittest analysis.test_block_survey -v` after producing the outputs. Run manifest hashes the protocol, script, all shard tables, raw annotation/provenance, PDF, workbook, cache and saved simulations; git metadata intentionally disabled because the operator prohibits all git commands.

[inference] Use this table as an exploratory, LD-reduced queue, not a selection result or validated production discovery ranking. Budget explicitly for the non-Akbari counts above, cross-tile LD and the uncalibrated block drift cost. First downstream step is independent modern-phenotype evidence acquisition and a linked block-null/stability audit before broad follow-up. Earlier authors' conclusions are untouched.

## POST-REGISTRATION CORRECTIONS AND EVIDENCE SENSITIVITY

[observation] Approved after the initial registered run, at approximately 2026-10-09T02:28Z. Earlier text above this appended survey is historical, NOT the current annotation. Initial TSV/spans/summary are retained as `data/derived/block_survey_initial_bounding*`. The scored labels, ranks, scores, tiling and MHC flags are unchanged; this appended annotated survey and current TSV supersede the first annotation attempt.

[observation] **Annotation bug:** a single min-start/max-end span per symbol/chrom/strand bridged disjoint refGene mappings. For example, RNU1-1, RNU1-2, RNU1-3, RNU1-4 and RNVU1-18 acquired an artificial chr1 interval [16840616,149224221), spanning ~132 Mb. The corrected rule is a SET of transcript intervals per symbol/chrom/strand, merging only overlapping or touching intervals (gap tolerance 0 bp). Every overlap must meet an actual transcript interval, not an invented inter-transcript gap.

[observation] Disjoint-interval diagnostics: 299 symbols / 308 symbol-chrom-strand groups have >1 merged interval. Spans >4 Mb are marked and excluded from overlap/nearest annotation: 0 intervals; omitted pathological symbols: none. The interval file retains these flags for audit. Pathology/gap tests in `analysis/test_block_survey.py` check a disjoint synthetic interval and actual refGene union coverage; no merged interval may bridge an uncovered source gap.

[observation] Annotation damage: 1367/10615 blocks change overlap/nearest annotation; 1259 registered neither-labelled blocks change; 3/10 previously listed novel tiles change. MHC flags unchanged: True; rank/score unchanged: False. Detailed lost/gained symbols are in `block_survey_summary.json` annotation_correction; corrected top lists appear in this appended survey.

[observation] **Empty moderate band is a panel CONSTRUCTION ARTIFACT:** cases were selected at |X|>5.45 and controls at |X|<3.16. The 3.16–5.45 band was never scored in these shards. Consequently registered scored-only moderate==strong; registered “neither” must NOT be read as no evidence anywhere in Akbari.

[observation] Supplemental POST-REGISTRATION evidence uses the local QC cache (538,274 variants, 31,939 strong, 356,335 in 3.16<=|X|<5.45; exact |X|=5.45 count 0). Cache moderate denotes the BAND, not moderate-including-strong. TSV fields: `block_max_absX_cache`, `block_has_strong_cache` (>5.45), `block_has_moderate_cache` (3.16<=|X|<5.45), `block_has_any_akbari_cache` (strong OR band), `n_cache_variants_in_block`; matching `_scored` fields distinguish the original panel. Multiple evidence classes can occur in one block. Null maximum if no cache rows; no-cache flags false means absence of evidence in this cache only.

|Frozen top N|Cache strong evidence|Cache moderate-only evidence|No Akbari cache evidence|Registered non-strong novel|Of those: cache moderate evidence|Of those: no cache evidence|
|---:|---:|---:|---:|---:|---:|---:|
|100|86|14|0|14|14|0|
|500|188|299|13|312|299|13|
|1000|244|719|37|756|719|37|

[inference] The genuinely no-Akbari-cache-evidence count is the honest novelty estimate in this available resource, smaller than the registered scored-only novel count wherever moderate evidence appears. This is evidence adjudication, NOT an enrichment analysis of the unscored moderate band: its variants have no measured ranking here. No variant was rescored or block reranked. Novel still does not mean neutral or globally unreported.

### Corrected highest-ranked no-cache-evidence hypotheses

[observation] Additional post-registration display: first ten frozen-ranked blocks with no strong or moderate evidence in the QC cache. This does not change the registered first-ten scored-neither display above; it distinguishes which of those survive cache adjudication. Gene/trait support remains unadjudicated.

|Rank|GRCh37 tile [start,end)|Representative|Score|Best abs X|Scored label|MHC|Genes (overlap unless stated)|Prior panel|479 tag in tile / rep|Cache evidence (POST-REG)|
|---:|---|---|---:|---:|---|---|---|---|---|---|
|119|chr6:[164000000,164250000)|6_164068784_G_C|1.692|1.160|neither|no|LOC102724152|no fixed-panel match|no / no|none_in_QC_cache|
|167|chr3:[66500000,66750000)|3_66745825_G_T|1.201|2.621|neither|no|LRIG1|no fixed-panel match|no / no|none_in_QC_cache|
|176|chr6:[81500000,81750000)|6_81607858_G_A|1.127|2.624|neither|no|null (no gene within 100 kb)|no fixed-panel match|no / no|none_in_QC_cache|
|181|chr13:[75000000,75250000)|13_75094164_A_G|1.090|1.388|neither|no|LINC00347;LINC00381|no fixed-panel match|no / no|none_in_QC_cache|
|275|chr10:[59000000,59250000)|10_59130616_C_T|0.617|1.697|neither|no|MIR3924|no fixed-panel match|no / no|none_in_QC_cache|
|308|chr18:[67500000,67750000)|18_67664150_A_T|0.520|2.224|neither|no|CD226;DOK6;RTTN|no fixed-panel match|no / no|none_in_QC_cache|
|372|chr15:[92250000,92500000)|15_92388014_C_T|0.363|2.333|neither|no|SLCO3A1|no fixed-panel match|no / no|none_in_QC_cache|
|383|chr8:[91750000,92000000)|8_91825458_TATATACACACACACACACACAC_T|0.328|2.243|neither|no|C8orf88;NECAB1|no fixed-panel match|no / no|none_in_QC_cache|
|407|chr11:[39750000,40000000)|11_39909127_A_T|0.252|1.835|neither|no|null (no gene within 100 kb)|no fixed-panel match|no / no|none_in_QC_cache|
|416|chr4:[58250000,58500000)|4_58495595_G_T|0.238|1.551|neither|no|LINC02380|no fixed-panel match|no / no|none_in_QC_cache|


---

## Appended regenerated survey — 2026-10-09T02:36:43.261078+00:00

# Exploratory block-level gene survey — 2026-10-09

[inference] **An exploratory hypothesis list for downstream modern-phenotype follow-up, NOT a selection discovery.** False positives are expected. Gene proximity and Akbari agreement do not establish function, causality, or independent validation.

## Registration and reproducibility

[observation] Protocol persisted first in `docs/block-survey-prereg.md` at 2026-10-09T02:18:36Z; annotation retrieval followed at 2026-10-09T02:19:41.016138+00:00. Inputs are the four scored shard TSVs, not the earlier small panel. Reproduce with `.venv/bin/python analysis/block_survey.py` (add `--fetch-annotation` only if the raw annotation is absent). Explicit TSV null token: `NA`.

[observation] UCSC hg19 refGene URL: https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/refGene.txt.gz; gzip CRC verified; SHA256 `28e7fe22cf1923a6887c7b60e6004c75987f810b8374d4d396021a34ac0410d1`. Provenance: `data/raw/annotation/hg19_refGene.provenance.json`. Parsed symbol/chromosome/interval-union/strand: `data/derived/block_survey_gene_spans.tsv`. All gene names in the candidate/control tables below are generated from annotated rows of `data/derived/block_survey.tsv`, never from PDF gene labels.

[observation] 181,933 finite scored variants, 31,939 strong variants; 10,615 occupied 250 kb tiles, 723 strong tiles and 723 moderate-including-strong tiles. No duplicate loci or missing input X. Annotation statuses: `{'overlap': 8966, 'no_gene_within_100kb': 979, 'nearest_within_100kb': 670}`.

[observation] The task's 181,933 and 31,939 variant counts agree. Approximately 723 strong blocks compares with the measured 723; tiling is explicitly one-based-pos to zero-based-start. Tiles are LD-reduced, not independent: adjacent tiles and MHC can share one signal.

## Top 20 blocks and local evidence

[observation] Strong means any scored variant |X|>5.45; moderate_only means max |X|>=3.16 but not strong. Moderate flags in TSV include strong. Prior panel = task-supplied canonical ancient-DNA support panel, also represented by the existing prioritizer controls; sparse illustrative controls, not an exhaustive literature review. A non-match is not evidence of no prior support. Coordinates/gene spans are zero-based half-open; representative locus strings contain one-based positions.

|Rank|GRCh37 tile [start,end)|Representative|Score|Best abs X|Scored label|MHC|Genes (overlap unless stated)|Prior panel|479 tag in tile / rep|Cache evidence (POST-REG)|
|---:|---|---|---:|---:|---|---|---|---|---|---|
|1|chr6:[31250000,31500000)|6_31359394_C_T|13.419|13.733|strong|yes|ATP6V1G2-DDX39B;DDX39B;HCG26;HCP5;HLA-B;LINC01149;LINC02571;MCCD1;MICA;MICA-AS1;MICB;MICB-DT;MIR6891|no fixed-panel match|yes / no|strong|
|2|chr6:[32500000,32750000)|6_32627830_C_T|12.900|13.704|strong|yes|HLA-DQA1;HLA-DQA2;HLA-DQB1;HLA-DQB1-AS1;HLA-DQB2;HLA-DRB1;HLA-DRB6;MIR3135B|no fixed-panel match|yes / no|strong|
|3|chr6:[28500000,28750000)|6_28689672_C_A|11.405|10.190|strong|yes|GPX5;LINC00533;ZBED9|no fixed-panel match|no / no|strong|
|4|chr6:[30750000,31000000)|6_30789647_G_C|11.125|10.596|strong|yes|DDR1;GTF2H4;HCG20;HCG21;LINC00243;LINC02570;MIR4640;MUC21;MUC22;MUCL3;SFTA2;VARS2|no fixed-panel match|yes / no|strong|
|5|chr6:[31500000,31750000)|6_31577099_A_T|11.010|11.446|strong|yes|ABHD16A;AIF1;APOM;ATP6V1G2;ATP6V1G2-DDX39B;BAG6;C6orf47;CLIC1;CSNK2B;DDAH2;DDX39B;DDX39B-AS1;GPANK1;LOC100287329;LST1;LTA;LTB;LY6G5B;LY6G5C;LY6G6C;LY6G6D;LY6G6E;LY6G6F;LY6G6F-LY6G6D;MIR4646;MIR6832;MPIG6B;MSH5;MSH5-SAPCD1;NCR3;NFKBIL1;PRRC2A;SAPCD1;SAPCD1-AS1;SNORA38;SNORD117;SNORD84;TNF;VARS1;VWA7|no fixed-panel match|yes / no|strong|
|6|chr6:[32000000,32250000)|6_32158319_A_G|10.304|12.460|strong|yes|AGER;AGPAT1;ATF6B;C4B;C4B_2;CYP21A2;EGFL8;FKBPL;GPSM3;LOC100507547;MIR6721;MIR6833;NOTCH4;PBX2;PPT2;PPT2-EGFL8;PRRT1;RNF5;TNXB;TSBP1-AS1|no fixed-panel match|yes / no|strong|
|7|chr6:[32750000,33000000)|6_32778571_A_ATT|10.078|9.699|strong|yes|BRD2;HLA-DMA;HLA-DMB;HLA-DOA;HLA-DOB;LOC100294145;PSMB8;PSMB8-AS1;PSMB9;TAP1;TAP2|no fixed-panel match|yes / no|strong|
|8|chr6:[31750000,32000000)|6_31757881_TA_T|9.797|13.835|strong|yes|C2;C2-AS1;C4A;C4B;C4B_2;CFB;CYP21A1P;DXO;EHMT2;HSPA1A;HSPA1B;HSPA1L;LOC110384692;LSM2;MIR1236;NELFE;NEU1;SKIV2L;SLC44A4;SNHG32;SNORD48;SNORD52;STK19;TNXA;TNXB;VARS1;ZBTB12|no fixed-panel match|yes / no|strong|
|9|chr6:[29000000,29250000)|6_29144532_C_T|9.733|8.248|strong|yes|LOC100129636;LOC101929006;OR2B3;OR2J1;OR2J2;OR2J3;OR2W1|no fixed-panel match|yes / no|strong|
|10|chr2:[136500000,136750000)|2_136608646_G_A|9.682|16.523|strong|no|DARS-AS1;DARS1;LCT;LCT-AS1;LOC107985946;MCM6;UBXN4|LCT;MCM6|yes / yes|strong|
|11|chr2:[136000000,136250000)|2_136161092_C_G|9.450|15.257|strong|no|ZRANB3|no fixed-panel match|no / no|strong|
|12|chr6:[32250000,32500000)|6_32343369_T_C|9.317|11.222|strong|yes|BTNL2;HCG23;HLA-DRA;HLA-DRB5;TSBP1;TSBP1-AS1|no fixed-panel match|yes / no|strong|
|13|chr2:[135250000,135500000)|2_135498924_G_C|9.269|14.895|strong|no|TMEM163|no fixed-panel match|no / no|strong|
|14|chr6:[26250000,26500000)|6_26338128_T_C|8.832|8.353|strong|no|BTN2A1;BTN2A2;BTN2A3P;BTN3A1;BTN3A2;BTN3A3;H2BC10;H2BC9;H3C7;H3C8;H4C8;LOC285819|no fixed-panel match|yes / no|strong|
|15|chr6:[30250000,30500000)|6_30345118_C_T|8.734|11.967|strong|yes|HCG17;HCG18;HLA-E;LINC02569;RPP21;TRIM39;TRIM39-RPP21|no fixed-panel match|yes / no|strong|
|16|chr2:[136250000,136500000)|2_136381348_G_A|8.682|15.815|strong|no|LOC107985946;MIR128-1;R3HDM1;UBXN4;ZRANB3|no fixed-panel match|yes / no|strong|
|17|chr6:[31000000,31250000)|6_31196862_G_A|8.649|12.082|strong|yes|C6orf15;CCHCR1;CDSN;HCG22;HCG27;HLA-C;MUC22;POU5F1;PSORS1C1;PSORS1C2;PSORS1C3;TCF19|no fixed-panel match|yes / no|strong|
|18|chr6:[33000000,33250000)|6_33029081_CCT_C|8.504|11.391|strong|yes|B3GALT4;COL11A2;HCG24;HCG25;HLA-DPA1;HLA-DPB1;HLA-DPB2;HSD17B8;MIR219A1;RING1;RPS18;RXRB;SLC39A7;VPS52;WDR46|no fixed-panel match|yes / no|strong|
|19|chr2:[135500000,135750000)|2_135689188_G_T|8.474|14.365|strong|no|ACMSD;CCNT2;CCNT2-AS1;MAP3K19;MIR5590|no fixed-panel match|no / no|strong|
|20|chr6:[29500000,29750000)|6_29542947_T_C|8.332|10.709|strong|yes|GABBR1;HLA-F;HLA-F-AS1;IFITM4P;LINC01015;MOG;OR2H2;SNORD32B;UBD;ZFP57|no fixed-panel match|yes / no|strong|

[inference] Overlapping genes are positional candidates only (union of actual transcript spans, not coding-exon intersection). Nearest fallback is allowed only when no gene overlaps, at <=100 kb. Nearest context is NOT positive evidence for a block. MHC flag = overlap chr6:[28500000,33500000). Multiple high-ranked MHC tiles remain one extended LD complex, not multiple independent biological confirmations.

## Prior ancient-DNA support panel: what surfaced?

|Control group|Best annotated block rank|Block representative|Annotated panel genes in that row|Label|
|---|---:|---|---|---|
|LCT/MCM6|10|2_136608646_G_A|LCT;MCM6|strong|
|SLC24A5|10510|15_48269895_T_G|SLC24A5|strong|
|SLC45A2|49|5_33951693_C_G|SLC45A2|strong|
|TYR|99|11_88911696_C_A|TYR|strong|
|OCA2/HERC2|2995|15_28225348_G_A|OCA2|neither|
|FADS1/FADS2|89|11_61618012_C_T|FADS1;FADS2|strong|
|TLR1/TLR6|34|4_38798648_C_A|TLR1;TLR6|strong|
|AMY1/AMY2 region|6814|1_104499856_C_T|AMY1A;AMY1B;AMY1C|neither|

[observation] Pre-registered predictions: top-100 strong enrichment is 12.63× (prediction >1: met); 5 top-100 tiles have a fixed-panel gene (prediction >=1: met). Neither success supplies independent validation.

[observation] Inherited variant-level limits remain: SLC24A5 ranked approximately 9,161 and the AMY1 SNP proxy approximately 698 in the earlier 10,641-variant panel. The table above reports NEW block ranks in the much denser merged survey and does not replace those prior variant ranks. The within-lineage residual design is structurally blind to carriage-driven change predominantly BETWEEN demes; it should not be presented as capturing those mechanisms. The amylase region SNP proxy is not a copy-number selection test.

## Akbari reported independent loci and phenotype evidence

[observation] Extracted all 479 unique indexed tag coordinates from PDF pages 52–91 (Figures S5.1–S5.40); all 479 match local GRCh37 cache coordinates and rounded |X| within 0.015. Output: `data/derived/block_survey_akbari_479_tags.tsv`. The 347 pruned-list identities / LD intervals were not extracted.

[inference] “479 tag in tile” checks whether an extracted tag coordinate lies in the tile, and “rep” checks exact representative locus ID against its verified cached tag ID. It is NOT overlap with an LD-defined independent-locus interval. Such intervals and membership of the more strictly pruned 347 set are not available from this extraction; `akbari_347_locus_overlap=NA` throughout. Absence of a tag from a tile does not imply independence or absence of the reported signal. Section 5 describes 479 tags (410 outside HLA +69 inside); Table S5.1 supplies D-prime pruning COUNTS, not the identities of a 347-tag subset. This is a difference in list/definition, not a claim that 347 is incorrect.

[observation] Table S5 sheet `STable5(696.polygenic.tests)` contains 696 aggregate phenotype tests, with beta/sign-beta analyses with and without HLA. It contains no locus or gene mapping. The strongest five absolute z values in the beta-without-HLA column are shown solely to describe available context, not to connect these phenotypes to the ranked tiles:

|Phenotype ID|Phenotype / coding description|z beta without HLA|z beta with HLA|
|---|---|---:|---:|
|Pan-UKBB.4004|Skin colour|-18.195|-18.249|
|Pan-UKBB.137|Hair colour (natural, before greying) / Black|-15.784|-15.783|
|UKB_460K.pigment_SUNBURN|Sunburn|14.070|14.023|
|Pan-UKBB.136|Hair colour (natural, before greying) / Dark brown|-13.005|-12.968|
|Pan-UKBB.4453|Average total household income before tax|9.462|9.514|

[inference] Aggregate Table S5 results cannot adjudicate any candidate gene; all locus-phenotype cells are explicitly null. **GWAS Catalog, eQTL, and functional annotation are NOT available locally and were not checked.** Their evidence statuses are explicit in every TSV row. RefGene transcript spans are positional annotations, not functional evidence.

## False-positive / follow-up cost at block level

[observation] Base prevalence in the 10,615 scored tiles: strong 723/10,615 = 6.81%; moderate including strong 723/10,615 = 6.81%. Denominator is the ascertained scored panel, not all genomic tiles.

|Top N|Strong (count / N)|Strong enrichment|Moderate incl. strong (count / N)|Moderate enrichment|Moderate only|Non-strong|Neither (strict novel)|MHC|479 tag tiles|
|---:|---|---:|---|---:|---:|---:|---:|---:|---:|
|100|86/100 (86.0%)|12.63×|86/100 (86.0%)|12.63×|0|14|14|19|54|
|500|188/500 (37.6%)|5.52×|188/500 (37.6%)|5.52×|0|312|312|20|118|
|1000|244/1000 (24.4%)|3.58×|244/1000 (24.4%)|3.58×|0|756|756|20|148|

[inference] Counts not strong are novel relative to the strong threshold; counts neither are the stricter novel subset, used below. These are follow-up costs, NOT proven false positives: non-significance is not neutrality. Strong and moderate agreement reuses Akbari data (same imputed genotypes), so enrichment is NOT independent validation; original variant ascertainment and tile maxima further favor strong/dense regions.

[observation] A calibrated LD-honest **block drift-null FP rate could not be computed**. Existing `analysis/prioritization_score.py::drift_simulation` / `--simulations` generates independent simulated single-locus genotype columns with no genomic coordinates or linked haplotypes. Arbitrary grouping would invent an LD model. The saved 300-draw drift+jump tail above the EARLIER top-100 variant threshold is 32/300 = 10.67% (agrees with quoted 10.67%). That legacy run scored 10,641 variants, NOT 181,933; all four merged shard runs used zero simulation replicates. Therefore the quoted value is not a merged-survey or block FDR.

[observation] Reusing the existing simulation draws, the following SINGLE-LOCUS tails exceed each observed BLOCK-rank score threshold. This is only a threshold diagnostic; it omits within-tile multiplicity and LD, so it is NOT the requested calibrated block false-positive probability or list precision:

|Top-block N|Block score threshold|Brownian drift only exceedances / 300|Drift + jump exceedances / 300|Replacement-only exceedances / 300|
|---:|---:|---|---|---|
|100|1.9056|0/300 (0.00%)|6/300 (2.00%)|0/300 (0.00%)|
|500|0.0863|0/300 (0.00%)|35/300 (11.67%)|19/300 (6.33%)|
|1000|-0.4146|1/300 (0.33%)|59/300 (19.67%)|43/300 (14.33%)|

[inference] Finite simulation size matters: even 0/300 has a one-sided binomial 95% upper tail probability about 0.994%, not zero. No observed-block FP count is inferable from these single-locus rates. A linked genotype/path simulation with observed block sizes and missingness would be needed to estimate block maxima under drift.

[observation] Earlier within-node bootstrap top-100 variant membership retention = 43.73% across 30 resamples (agrees with approximately 43.7%). This has not been recomputed at block level or on the merged panel; it warns that ranks near cuts are unstable, not that block retention equals that number.

## Highest-ranked novel blocks: bounded follow-up hypotheses

[observation] The first ten tiles with neither strong nor moderate scored-variant X are listed mechanically (not selected by appealing gene stories). “Interesting” here means rank/proximity under the pre-registered exploratory protocol only. The broader non-strong candidate set additionally includes moderate-only tiles.

|Rank|GRCh37 tile [start,end)|Representative|Score|Best abs X|Scored label|MHC|Genes (overlap unless stated)|Prior panel|479 tag in tile / rep|Cache evidence (POST-REG)|
|---:|---|---|---:|---:|---|---|---|---|---|---|
|58|chr1:[52000000,52250000)|1_52025112_T_TACAC|3.746|3.153|neither|no|CALR4P;OSBPL9|no fixed-panel match|no / no|moderate_only|
|60|chr1:[31000000,31250000)|1_31033851_A_G|3.580|2.930|neither|no|LAPTM5;MATN1;MATN1-AS1;MIR4420|no fixed-panel match|no / no|moderate_only|
|68|chr16:[4750000,5000000)|16_4958961_A_C|2.950|2.538|neither|no|ANKS3;C16orf71;GLYR1;PPL;ROGDI;SEPTIN12;SMIM22;UBN1;ZNF500|no fixed-panel match|no / no|moderate_only|
|69|chr12:[114750000,115000000)|12_114917493_C_T|2.947|1.784|neither|no|TBX5;TBX5-AS1|no fixed-panel match|no / no|moderate_only|
|76|chr7:[155750000,156000000)|7_155996962_C_T|2.712|2.574|neither|no|LOC389602|no fixed-panel match|no / no|moderate_only|
|79|chr3:[138500000,138750000)|3_138651944_G_A|2.468|2.596|neither|no|FOXL2;FOXL2NB;LINC01391;PIK3CB;PRR23A;PRR23B|no fixed-panel match|no / no|moderate_only|
|81|chr4:[1000000,1250000)|4_1213988_C_CAA|2.444|2.217|neither|no|CTBP1;CTBP1-AS;CTBP1-DT;FGFRL1;LOC100130872;LOC105374344;RNF212;SPON2;TMED11P|no fixed-panel match|no / no|moderate_only|
|82|chr5:[161000000,161250000)|5_161229995_G_A|2.426|2.959|neither|no|GABRA6|no fixed-panel match|no / no|moderate_only|
|85|chr9:[105000000,105250000)|9_105230080_C_A|2.175|2.503|neither|no|LINC00587 (nearest, 31918 bp)|no fixed-panel match|no / no|moderate_only|
|87|chr10:[55500000,55750000)|10_55669459_G_T|2.153|2.908|neither|no|PCDH15|no fixed-panel match|no / no|moderate_only|

[speculation] These tiles could motivate a small modern-phenotype follow-up pilot. Their ranks may instead reflect drift/jitter, carriage residuals, LD with a neighboring signal, or dense-panel maxima; **the local evidence base cannot adjudicate them**. Novel means unlabeled in this scored panel, not globally unreported selection. No causal gene or trait association is claimed.

## Validation and bounded recommendation

[observation] `analysis/test_block_survey.py` covers tile-boundary/tie rules, any-variant labels and missing-X unknowns, MHC boundaries, overlap/nearest distance limits, and output provenance/gene traceability. Run `.venv/bin/python -m unittest analysis.test_block_survey -v` after producing the outputs. Run manifest hashes the protocol, script, all shard tables, raw annotation/provenance, PDF, workbook, cache and saved simulations; git metadata intentionally disabled because the operator prohibits all git commands.

[inference] Use this table as an exploratory, LD-reduced queue, not a selection result or validated production discovery ranking. Budget explicitly for the non-Akbari counts above, cross-tile LD and the uncalibrated block drift cost. First downstream step is independent modern-phenotype evidence acquisition and a linked block-null/stability audit before broad follow-up. Earlier authors' conclusions are untouched.

## POST-REGISTRATION CORRECTIONS AND EVIDENCE SENSITIVITY

[observation] Approved after the initial registered run, at approximately 2026-10-09T02:28Z. Earlier text above this appended survey is historical, NOT the current annotation. Initial TSV/spans/summary are retained as `data/derived/block_survey_initial_bounding*`. The scored labels, ranks, scores, tiling and MHC flags are unchanged; this appended annotated survey and current TSV supersede the first annotation attempt.

[observation] **Annotation bug:** a single min-start/max-end span per symbol/chrom/strand bridged disjoint refGene mappings. For example, RNU1-1, RNU1-2, RNU1-3, RNU1-4 and RNVU1-18 acquired an artificial chr1 interval [16840616,149224221), spanning ~132 Mb. The corrected rule is a SET of transcript intervals per symbol/chrom/strand, merging only overlapping or touching intervals (gap tolerance 0 bp). Every overlap must meet an actual transcript interval, not an invented inter-transcript gap.

[observation] Disjoint-interval diagnostics: 299 symbols / 308 symbol-chrom-strand groups have >1 merged interval. Spans >4 Mb are marked and excluded from overlap/nearest annotation: 0 intervals; omitted pathological symbols: none. The interval file retains these flags for audit. Pathology/gap tests in `analysis/test_block_survey.py` check a disjoint synthetic interval and actual refGene union coverage; no merged interval may bridge an uncovered source gap.

[observation] Annotation damage: 1367/10615 blocks change overlap/nearest annotation; 1259 registered neither-labelled blocks change; 3/10 previously listed novel tiles change. MHC flags unchanged: True; rank/score unchanged: True. Detailed lost/gained symbols are in `block_survey_summary.json` annotation_correction; corrected top lists appear in this appended survey.

[observation] **Empty moderate band is a panel CONSTRUCTION ARTIFACT:** cases were selected at |X|>5.45 and controls at |X|<3.16. The 3.16–5.45 band was never scored in these shards. Consequently registered scored-only moderate==strong; registered “neither” must NOT be read as no evidence anywhere in Akbari.

[observation] Supplemental POST-REGISTRATION evidence uses the local QC cache (538,274 variants, 31,939 strong, 356,335 in 3.16<=|X|<5.45; exact |X|=5.45 count 0). Cache moderate denotes the BAND, not moderate-including-strong. TSV fields: `block_max_absX_cache`, `block_has_strong_cache` (>5.45), `block_has_moderate_cache` (3.16<=|X|<5.45), `block_has_any_akbari_cache` (strong OR band), `n_cache_variants_in_block`; matching `_scored` fields distinguish the original panel. Multiple evidence classes can occur in one block. Null maximum if no cache rows; no-cache flags false means absence of evidence in this cache only.

|Frozen top N|Cache strong evidence|Cache moderate-only evidence|No Akbari cache evidence|Registered non-strong novel|Of those: cache moderate evidence|Of those: no cache evidence|
|---:|---:|---:|---:|---:|---:|---:|
|100|86|14|0|14|14|0|
|500|188|299|13|312|299|13|
|1000|244|719|37|756|719|37|

[inference] The genuinely no-Akbari-cache-evidence count is the honest novelty estimate in this available resource, smaller than the registered scored-only novel count wherever moderate evidence appears. This is evidence adjudication, NOT an enrichment analysis of the unscored moderate band: its variants have no measured ranking here. No variant was rescored or block reranked. Novel still does not mean neutral or globally unreported.

### Corrected highest-ranked no-cache-evidence hypotheses

[observation] Additional post-registration display: first ten frozen-ranked blocks with no strong or moderate evidence in the QC cache. This does not change the registered first-ten scored-neither display above; it distinguishes which of those survive cache adjudication. Gene/trait support remains unadjudicated.

|Rank|GRCh37 tile [start,end)|Representative|Score|Best abs X|Scored label|MHC|Genes (overlap unless stated)|Prior panel|479 tag in tile / rep|Cache evidence (POST-REG)|
|---:|---|---|---:|---:|---|---|---|---|---|---|
|119|chr6:[164000000,164250000)|6_164068784_G_C|1.692|1.160|neither|no|LOC102724152|no fixed-panel match|no / no|none_in_QC_cache|
|167|chr3:[66500000,66750000)|3_66745825_G_T|1.201|2.621|neither|no|LRIG1|no fixed-panel match|no / no|none_in_QC_cache|
|176|chr6:[81500000,81750000)|6_81607858_G_A|1.127|2.624|neither|no|null (no gene within 100 kb)|no fixed-panel match|no / no|none_in_QC_cache|
|181|chr13:[75000000,75250000)|13_75094164_A_G|1.090|1.388|neither|no|LINC00347;LINC00381|no fixed-panel match|no / no|none_in_QC_cache|
|275|chr10:[59000000,59250000)|10_59130616_C_T|0.617|1.697|neither|no|MIR3924|no fixed-panel match|no / no|none_in_QC_cache|
|308|chr18:[67500000,67750000)|18_67664150_A_T|0.520|2.224|neither|no|CD226;DOK6;RTTN|no fixed-panel match|no / no|none_in_QC_cache|
|372|chr15:[92250000,92500000)|15_92388014_C_T|0.363|2.333|neither|no|SLCO3A1|no fixed-panel match|no / no|none_in_QC_cache|
|383|chr8:[91750000,92000000)|8_91825458_TATATACACACACACACACACAC_T|0.328|2.243|neither|no|C8orf88;NECAB1|no fixed-panel match|no / no|none_in_QC_cache|
|407|chr11:[39750000,40000000)|11_39909127_A_T|0.252|1.835|neither|no|null (no gene within 100 kb)|no fixed-panel match|no / no|none_in_QC_cache|
|416|chr4:[58250000,58500000)|4_58495595_G_T|0.238|1.551|neither|no|LINC02380|no fixed-panel match|no / no|none_in_QC_cache|
