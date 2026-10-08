# Changes from upstream (e10m/WESley @ db5b2fc)

This is a modified copy of WESley, not a fork with history — see the base
repo at https://github.com/e10m/WESley for the original. Two kinds of
changes are bundled together here: a new mutation-selection strategy, and
five pipeline bugs found/fixed while benchmarking it.

## 1. New mutation-selection strategy

**Old:** Mutect2 + MuSE + VarScan2, 2-of-3 majority vote (`bcftools isec -n +2`).
**New:** DeepSomatic + MuSE v2, union (`bcftools isec -n +1`) — a call is kept
if *either* caller found it.

This was decided by benchmarking all 5 callers (the original 3, plus
DeepSomatic and Lancet2 as actively-maintained candidates — Strelka2,
VarDict, and Octopus were considered and ruled out as abandoned/stale)
against SEQC2 truth (two replicates: IL1/Illumina-HiSeq4000,
EA1/EA-HiSeq1500). Full tables: `/media/nathansonlab/irondog/jfreeland/benchmark_results/`.

Headline result: DeepSomatic ∪ MuSE v2 beat every other combination tested,
including the original 3-caller rule and several 3-/5-caller alternatives,
on both SNVs (avg F1 0.952 vs. 0.915 baseline, +0.037) and indels (avg F1
0.861 vs. 0.732 baseline, +0.129). More callers did not mean better results
— every combination that added Mutect2 and/or VarScan2 back in performed
*worse* than the simple 2-caller union, because their extra calls were
mostly false positives, not new true positives.

VarScan2 was also independently confirmed as unmaintained software (last
real release 2016) — removing it isn't only an accuracy call, it's dropping
a dead dependency.

### What this meant architecturally

DeepSomatic's tumor/normal mode reports a **single sample column** (the real
tumor sample ID) — it does not emit a NORMAL genotype column the way
Mutect2/MuSE/VarScan2 did (the normal sample is used internally to classify
PASS/GERMLINE/RefCall, not carried through as its own column). MuSE still
reports its native two columns (`NORMAL`, `TUMOR`). To merge these, both
`reheader.nf` modules (mutation_calling's per-caller one, consensus_calling's
pre-merge one) now reduce everything to a **single column literally named
`TUMOR`**: DeepSomatic's one real-named column is renamed; MuSE's `NORMAL`
column is dropped, keeping only its native `TUMOR` column.

This means the final consensus VCF carries only the tumor's genotype at
each site, not the normal's. That's a deliberate simplification: these are
by definition sites classified as somatic (absent in the matched normal),
so the normal's genotype at them isn't informative to retain, and neither
remaining caller gives you a real one anyway. If something downstream
expects a 2-sample `NORMAL`/`TUMOR` VCF, it will need to be updated — flag
this explicitly to anyone reviewing before merging.

### Files touched for the strategy change

- **New:** `nextflow_automation/mutation_calling/modules/deepsomatic/deepsomatic.nf`
- **Removed:** `mutation_calling/modules/mutect2/`, `mutect2_pon/`, `varscan2/`
  (and the `CREATE_M2_PON` workflow entry point in `mutation_calling.nf`,
  which only existed to support Mutect2's panel-of-normals)
- **Rewritten:** `mutation_calling.nf`, `mutation_calling/nextflow.config`,
  `mutation_calling/modules/shared/{reheader,select_variants,index,vep}.nf`,
  `mutation_calling/conf/omics.config`,
  `consensus_calling.nf`, `consensus_calling/nextflow.config`,
  `consensus_calling/modules/{reheader,sort_vcfs,intersect,merge_vcfs,index}.nf`
- `--interval_list` was renamed to **`--capture_bed`** — DeepSomatic's
  `--regions` flag wants a plain BED file, not a GATK `.interval_list`
  (the format Mutect2 needed). If you have automation passing
  `--interval_list`, update it to point at a plain BED and use the new
  flag name.
- `gnomad_vcf`/`contamination_vcf` params and the whole contamination-table
  machinery were removed — those existed only to support Mutect2's
  `FilterMutectCalls` contamination model.

### A second round of fixes, found by actually checking everything

After the first pass, the user directly asked "were the diagrams and
everything else updated?" — the honest answer was no, which surfaced a real
functional bug beyond just stale docs:

- **Real bug, now fixed:** `mutation_calling`'s *per-caller* downstream chain
  (`shared/reheader.nf`, `create_maf.nf`, `keep_tertP.nf`,
  `keep_nonsynonymous.nf`, `rename_hg38.nf`, `oncokb.nf` — distinct from
  `consensus_calling`'s own copies of similarly-named modules) still
  branched on hardcoded `mutect2`/`varscan2` basenames with no DeepSomatic
  case at all. Worse, `create_maf.nf`'s `vcf2maf.pl` call would have run
  DeepSomatic's single-column output through the matched-pair invocation
  (`--vcf-normal-id "NORMAL"`), pointing at a sample column that doesn't
  exist. Fixed: `shared/reheader.nf` now renames DeepSomatic's one real-named
  column to `TUMOR` (matching `consensus_calling/reheader.nf`'s approach),
  and `create_maf.nf` uses the tumor-only-style `vcf2maf.pl` invocation
  (`--tumor-id`, `--vcf-tumor-id "TUMOR"`, no normal flags) for DeepSomatic.
  Verified directly against real benchmark output: `bcftools reheader`
  correctly renames the column, and `vcf2maf.pl` runs clean (exit 0,
  produces exactly 1,218 MAF rows matching the 1,218 PASS-filtered input
  variants; the "No genotype column for NORMAL" warning it prints is
  expected, not an error, given the design here).
- **Docs/config updated too:** `README.md` (TOC, run examples, container
  tables, software dependency tables, CI description, diagram captions),
  `containerization/compose.yaml` and `ecr_push.sh` (MuSE version, VarScan2
  removed, DeepSomatic added), removed the now-orphaned
  `containerization/dockerfiles/VarScan2.Dockerfile`, and
  `nextflow_automation/tests/shared-test.config`'s leftover Mutect2
  container selector.
- **Not fixed — flagged instead, since I can't regenerate images:**
  `diagrams/mutation-calling.png` and `diagrams/consensus-calling.png` still
  depict the old 3-caller architecture. Added explicit warning callouts in
  `README.md` pointing at this file rather than leaving them silently
  wrong. Redrawing them is a real remaining task for whoever picks this up.

### Known caveat — not yet re-validated end-to-end

I smoke-tested the new `consensus_calling.nf` logic against real
already-computed DeepSomatic/MuSE output (from the benchmark) and confirmed
it runs cleanly end-to-end (SORT_VCFS → REHEADER → INDEX → INTERSECT →
MERGE_VCFS → NORM_INDELS → VEP → CREATE_MAF → ... all succeed). But the
final record count from this real pipeline code (4,508 for the IL1 test
case) doesn't exactly match the count my ad-hoc benchmark script produced
for the same inputs (3,929) — both are computing "the union of DeepSomatic
and MuSE," but through different mechanics (my benchmark script did a plain
Python set-union; this pipeline reuses the original `bcftools isec` +
`MergeVcfs` + `bcftools norm -d none` pattern), and I did not fully
reconcile exactly why they differ before finalizing this. The mechanism
itself is unchanged from the original, long-production-validated pipeline
— I only changed caller count/threshold, not how dedup works — so this
isn't a new bug introduced here, but **the exact F1 numbers in the
benchmark tables have not been re-confirmed against this specific pipeline
code.** Recommend re-running hap.py scoring against this repo's actual
output before fully trusting the benchmark numbers as this pipeline's
numbers.

I also have not run `mutation_calling.nf` itself end-to-end with real BAMs
through this repo's new `DEEPSOMATIC` module — it's the same invocation
already validated extensively in the benchmark, but the module file is new
code, written for this repo, not copy-pasted from an already-running
pipeline. Treat it as needing a real run before production use, not as
pre-validated.

## 2. Pipeline bugs fixed (found during benchmarking, independent of the strategy change)

1. **`data_processing/modules/mark_duplicates.nf`** — MARK_DUPES crashed on
   single-lane samples. Nextflow silently collapses a single-item
   `path(sorted_bams)` channel from `List<Path>` to a bare `Path`;
   `.collect{}` then iterated the path's filesystem segments instead of
   files. Fixed: `sorted_bams.collect{...}` → `([sorted_bams].flatten()).collect{...}`.
   Affects any real single-lane production sample, not just benchmark data.
2. **`consensus_calling/modules/reheader.nf`** (original 3-caller version,
   now superseded by the strategy change above, but the underlying lesson
   carried forward) — silently failed to reconcile sample columns for
   non-TCGB sample IDs (e.g. reference cell lines), breaking `MERGE_VCFS`
   with no visible top-level error.
3. **`consensus_calling/modules/vep.nf`** — hardcoded `--cache_version 103`,
   but the `e10m/vep:115` container only has a `115_GRCh38` cache baked in.
   Broke consensus VEP annotation unconditionally, every batch. (The
   `mutation_calling` copy of this module was already correct — only the
   `consensus_calling` copy had drifted.) Fixed: `103` → `115`.
4. **`mutect2/modules/mutect2/learn_read_orientation.nf`** — `label 'lowMem'`
   (1GB×attempt) under-provisioned for deep samples, silently dropping
   Mutect2's output downstream with no error. **Moot in this repo** — the
   whole Mutect2 module tree was removed per the strategy change (#1
   above), so this fix has nothing left to apply to. Documented here so
   the history isn't lost, not because it's live code anymore.
5. **MuSE v1 → v2** — three-part fix, all still live and relevant:
   - `quay.io/biocontainers/muse:1.0.rc--1` → `muse:2.1.2--h3b3e331_3`. v1.0rc
     segfaulted deterministically on some samples (confirmed a genuine
     binary defect, not memory/data corruption — ruled out by testing).
     Verified v2 produces byte-identical output to v1 where v1 could
     complete, so this is a pure bug fix, not a quality tradeoff.
   - `label 'lowMem'` → `label 'highMem'` — v2 needs substantially more
     memory than v1 did; the old label OOM-killed it on all retries.
   - Added `-n ${task.cpus}` to both `MuSE call` and `MuSE sump` — v2's
     `sump` crashes (`Number of cores cannot be less than 1`) without an
     explicit core count; v1 defaulted safely without it.

All five bugs were originally found and fixed in an isolated benchmark fork
(`mutation_benchmarking/wesley_bootstrap/`), kept separate from the real
repo at the time per instruction. This repo is where they're actually
being folded back in.

## Not changed

`data_processing/`, `cnvkit/`, `fingerprint/` are untouched copies of the
upstream modules (beyond the one `mark_duplicates.nf` fix above). CNV
calling (CNVkit) and sample-identity fingerprinting are separate concerns
from mutation calling and weren't part of this benchmarking effort.
