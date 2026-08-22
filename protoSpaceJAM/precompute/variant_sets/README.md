# Building a variant set

protoSpaceJAM designs against a reference genome. Edits happen in a cell line, whose genome is not
the reference. Reference homology arms carry bases the cells do not have, and guides are chosen and
scored on sequence that may not exist in the target.

A **variant set** fixes that. It is a named, prebuilt asset holding one sample's genotype on one
haplotype, in two interchangeable forms:

| | path | what it is |
|---|---|---|
| substituted genome | `genome_files/fa_pickle/<name>/` | the reference pickles with that sample's substitutions written in. ~2.9 GB for GRCh38 (unchanged contigs are hardlinked, not copied) |
| sidecar | `genome_files/variant_sets/<name>/` | the differences only, plus every variant we could *not* write and why. ~77 MB |

Both are under `genome_files/`, which is gitignored — they are built, not committed.

At design time:

```
protoSpaceJAM --variant_genome <name>    # read the substituted pickles   <- what the web portal uses
protoSpaceJAM --variant_set    <name>    # patch each slice from the sidecar as it is fetched
```

The two are verified to produce byte-identical designs (see *Verification* below), so the choice is
a performance/inspection one, not a scientific one. **The sidecar is loaded either way** — it is
what supplies the `variant_*` columns in `result.csv`, the per-job `variants_report.csv`, the
GenBank `variation` features and the guide penalties. The substituted genome only supplies sequence.

---

## The pipeline

Steps 1–3 need `bcftools`, `samtools`, `bgzip` and `tabix`. Steps 4–6 are stdlib + Biopython only —
no pysam, no cyvcf2, no bcftools at build or run time (BGZF is gzip-readable for a sequential pass).

```
1_make_homopolymer_bed.py     reference FASTA        -> homopolymer runs BED     (per genome, reused)
2_filter_callset.sh           joint VCF + that BED   -> filtered single-sample VCF
3_make_nocall_bed.sh          joint VCF              -> BED of ./. positions     (per sample)
4_make_variant_genome.py      filtered VCF + no-call -> substituted genome + sidecar
5_validate_variant_genome.py  the substituted genome -> PASS/FAIL
6_validate_runtime_patch.py   sidecar vs the genome  -> PASS/FAIL
```

Step 3 must read the **original** VCF: step 2 extracts with `--exclude-uncalled`, which is exactly
what destroys the information step 3 is looking for.

### Worked example — `GRCh38_KOLF2.1J_hap1`

The asset currently shipped. Source data lives outside this repo (see *Where the inputs live*).

```bash
cd protoSpaceJAM                       # the package dir, so genome_files/ resolves
PY=~/miniforge3/envs/protoSpaceJAM-portal/bin/python
VS=precompute/variant_sets
D=/path/to/variant/full_callset_to_use

# 1. reference-only, ~2 min; build once per genome
$PY $VS/1_make_homopolymer_bed.py --fasta $D/../GRCh38.fa --min_len 5 \
    | bgzip > $D/GRCh38.ensembl.homopolymer_ge5.bed.gz
tabix -p bed $D/GRCh38.ensembl.homopolymer_ge5.bed.gz

# 2. ~6 min
$VS/2_filter_callset.sh --vcf $D/czb_kolf.joint.GRCh38.vcf.gz \
    --sample czML1175_KOLF2.1J-Parental --fasta $D/../GRCh38.fa \
    --hp_bed $D/GRCh38.ensembl.homopolymer_ge5.bed.gz \
    --out $D/kolf2.1j.filtered.vcf.gz

# 3. seconds
$VS/3_make_nocall_bed.sh --vcf $D/czb_kolf.joint.GRCh38.vcf.gz \
    --sample czML1175_KOLF2.1J-Parental --out $D/kolf2.1j.nocall.bed.gz

# 4. ~8 min
$PY $VS/4_make_variant_genome.py \
    --vcf $D/kolf2.1j.filtered.vcf.gz --nocall_bed $D/kolf2.1j.nocall.bed.gz \
    --sample czML1175_KOLF2.1J-Parental --genome_ver GRCh38 \
    --out_genome_ver GRCh38_KOLF2.1J_hap1 --checksum --force

# 5, 6. ~15 s together
$PY $VS/5_validate_variant_genome.py --out_genome_ver GRCh38_KOLF2.1J_hap1
$PY $VS/6_validate_runtime_patch.py  --variant_set   GRCh38_KOLF2.1J_hap1 --windows 40000
```

`--chroms 21,22` on steps 4–6 gives a ~5 s iteration loop.

Result, from a 4,612,205-record callset:

```
applied              2,632,109   substitutions written into hap1 (1 per 1,173 bp)
recorded-not-applied 1,980,096   other_haplotype 1,115,389 | indel_selected_hap 552,610
                                 indel_other_hap 269,712 | homopolymer_snv 39,237
                                 unphased_het 2,002 | indel_unphased 1,142 | symbolic_allele 4
sequence_known_wrong   555,758   windows where the emitted sequence is NOT the real haplotype
no-call regions        400,732   690,857 bp where the sample has no genotype at all
```

`applied + recorded-not-applied` must equal the callset's record count exactly. That identity is the
cheapest check that the build read the whole file.

---

## Substitutions only

**Only length-preserving substitutions are ever written.** Indels are recorded and reported, never
applied. This is not a simplification to be tidied up later — it is forced by everything else in
`genome_files/`, none of which is regenerated per variant set:

- `parsed_gff3/<genome>/loc2posType.pickle`, `ENST_info/`, `ENST_codonPhases/`
- the `start`/`end` columns of `precomputed_gRNAs/gRNAs_<PAM>/gRNA_<genome>/`
- four `pos - Lstart` string-index subtractions in `util/hdr.py`
  (`get_trunc_gRNA`, `get_ins2cut_seq`, `make_gRNA_lowercase`) and the phase-list construction in
  `util/utils.py get_HDR_template`

All of it is keyed to reference coordinates. Insert or delete one base and all of it goes out of
register, silently. Handling indels properly means remapping those assets — a separate project.

What you get instead: any window overlapping an indel on the selected haplotype is **flagged**.
`variant_annot.SEQUENCE_UNRELIABLE` is the authoritative set of reasons that make a window
untrustworthy, and `result.csv` gains `donor_sequence_unreliable` / `arms_sequence_unreliable`
warnings. A design there is not silently wrong; it is loudly uncertain.

---

## The filter, and why each threshold

One rule decides every quality question here:

> At a site where the caller emitted a variant, the pipeline must write either the ALT or the REF.
> Emitting REF is wrong precisely when the call is true. So write whichever allele is more likely:
> **apply when P(true) > 0.5.**

P(true) is estimated from Ti/Tv, taking 2.1 for real germline SNVs and 0.5 for random substitution
error (there are twice as many possible transversions as transitions, so pure noise sits at 0.5).
Measured on the KOLF2.1J callset:

| band | Ti/Tv | P(true) | verdict |
|---|---|---|---|
| GQ 0–19 | ≈1.0 | ~0.35 | drop |
| GQ 20–29 | 1.475 | **0.76** | keep — reference would be wrong 3 times in 4 |
| GQ ≥30 | 2.08 | ~0.97 | keep |
| homopolymer run 5–6 bp | 1.135 | **0.58** | apply — so the homopolymer cut belongs at 7, not 5 |
| homopolymer run 7–9 bp | 0.440 | **~0** | emit reference |
| homopolymer run ≥10 bp | 0.430 | **~0** | emit reference |

Which gives `--gq 20` in step 2 and `--homopolymer_ref_min_len 7` in step 4.

The other thresholds are correctness guards rather than quality trade-offs:

| | | |
|---|---|---|
| `--dp_min 25` / `--dp_max 120` | both tails are error-enriched (DP<10 → Ti/Tv 0.798, DP≥150 → 0.931) | **MT is exempt from the cap.** Mitochondrial depth runs 188–434 from copy number; a flat cap deletes the entire contig, and the first build lost all of MT this way |
| drop het calls on Y and X-nonPAR | for a male line these are impossible; 22.9% of raw chrY non-PAR calls were het, pure X-Y mismapping | PAR coordinates are assembly-specific |
| `norm -m -any` (split multiallelics) | at a `1\|2` site *neither* haplotype is reference. Dropping the site emits reference and is wrong on **both** | splitting turns it into a clean `1\|0` plus `0\|1`, and exactly one of the pair lands on the selected haplotype |
| GQ ≥ 40 rejected | Ti/Tv gains 0.03, but het/hom collapses 1.29 → 0.45 | hets inherently carry lower GQ, so a high cut strips exactly the population a haplotype-aware feature exists to serve |

### Homopolymers: substitutions only

Step 4 declines to write a substitution whose `HPLEN` tag is ≥ 7, emitting the reference base and
recording the call as `homopolymer_snv`. 39,237 of them on KOLF2.1J hap1.

Two things about this are deliberate and worth not undoing:

**A quality gate cannot do this job.** 84.4% of the runs-7–9 SNVs and 62.6% of the runs-≥10 SNVs
carry GQ ≥ 30. Homopolymer slippage is a *systematic* error — the caller is confidently wrong — so
no GQ threshold removes it. That is the entire reason a separate rule exists.

**It stops at substitutions. Indels in homopolymers keep their flag.** Runs ≥7 hold 26% of all
indels and it is tempting to discard them as slippage too. Three reasons not to:

1. Ti/Tv is a transition/transversion ratio, defined only for SNVs. It gives no estimate of indel
   FDR. The 0.43 above describes the 27,105 SNVs in those runs, not the 183,903 indels beside them.
2. Long homopolymers are also where real germline length polymorphism concentrates — the same
   instability drives the true variation and the sequencer error. The "probably an artifact" prior
   that is overwhelming for an SNV is genuinely weak for an indel.
3. The costs are asymmetric. Suppressing a bad SNV costs one possibly-wrong base in a 1 kb arm.
   Suppressing an indel *flag* means shipping an arm of the wrong length and saying nothing.
   Indels are never applied anyway, so the only thing on the table was the warning — and the
   warning is the valuable part.

`homopolymer_snv` is also the one skip reason deliberately **not** in `SEQUENCE_UNRELIABLE`. Every
other reason means "we could not represent this"; this one means "we looked at the call and decided
the reference base was the better bet". Flagging the window would contradict the decision.

### No-call regions

`bcftools view --exclude-uncalled` in step 2 deletes every site where the sample is `./.`, which
makes "we do not know" indistinguishable from "it matches the reference". Step 3 recovers them and
step 4 writes them into a `nocall/` shard; `SKIP_NO_CALL_REGION` is in `SEQUENCE_UNRELIABLE`, so a
donor overlapping one is reported unreliable rather than confidently reference. On KOLF2.1J that is
400,732 regions — about 12.6% of 1 kb donors overlap one.

### What the filter still hides

The filter is **hard**: records below the GQ/DP thresholds are gone, not tagged. A removed variant
enters neither the `applied` nor the `recorded-not-applied` bucket, so its window is emitted as
reference and reported reliable. At GQ≥20 that is 830,993 records, ~1 per 3,700 bp, so roughly a
quarter of 1 kb donors contain one. Lowering the gate from 30 to 20 shrank the blind spot by about
a third; it did not remove it.

Fixing it properly means a soft filter —

```bash
bcftools filter -s LOWQUAL -m + -e 'GQ<20 || DP<25 || (CHROM!="MT" && DP>120)' ...
```

— plus a `low_confidence` branch in `classify()` that lands in `recorded-not-applied` and in
`SEQUENCE_UNRELIABLE`. Same emitted sequence, uncertainty stays visible. Not done.

---

## Verification

Three layers, and they are not redundant.

**`5_validate_variant_genome.py`** — per contig: length/id/name preserved, provenance recorded,
every applied position read REF before and reads ALT after with soft-mask case intact, every
skipped variant untouched, and *minimality*: revert the known applied positions in the substituted
sequence and assert byte-identity with the reference. That last check is what proves the 39,237
homopolymer-suppressed sites really do emit reference — nothing outside the applied set differs.
Check 5 additionally fails the asset if the manifest claims a policy that left no trace (an HPLEN
annotation lost in a rebuild, a `--nocall_bed` that produced zero regions).

**`6_validate_runtime_patch.py`** — random windows on both strands, 20–2000 bp, including windows
straddling variants at every offset and clamped to contig edges, fetched both ways and compared.
1,000,000 windows, 602,889 containing variants, 0 mismatches.

**`tests/run_full_test_variant.py`** — the same equivalence at design level: the full OpenCell
design (1,236 rows) run three ways — `--variant_genome`, `--variant_set`, and no flag — and diffed.
The first two must be byte-identical in `result.csv`, `variants_report.csv` and every GenBank file;
the third must reproduce the committed reference ground truth, because "inert unless asked for" is
the constraint the whole feature was built under. ~20 min; the three runs are launched concurrently.

For a fast check, `tests/run_quick_test_variant.py` covers the same chain on six loci in ~10 s.

---

## Where the inputs live

The source VCFs, the reference FASTA and the homopolymer BEDs are **not in this repo and should not
be**: several GB, and human genotype data. On the build machine they sit in the portal checkout
under `variant/`, which that repo gitignores, alongside two handoff documents recording how the
callset was chosen and what was measured.

Nothing is lost by that, because provenance travels with the built asset. Every variant set's
`manifest.json` records the absolute path and sha256 of the source VCF and of the no-call BED, the
sample, the haplotype, every policy flag, and per-contig counts. Any result can be traced back to
the exact genotype data behind it, and step 4 re-runs a REF-concordance preflight on every build so
a wrong-assembly VCF cannot be used by accident.

**Scientific caveat to carry to release:** a variant set named for a cell line describes *the clones
that were sequenced*. Other labs' stocks of the same line have drifted, so the asset is an
approximation for anyone outside the group that produced it.

---

## What to check first if something looks wrong

1. `tests/run_quick_test_variant.py` — covers the whole chain in ~10 s.
2. The pre-existing quick tests (`insertion`, `insertion_dsDNA`, `SNP`) must still pass byte-identically.
   If they do not, something is emitting variant columns or applying a variant weight with no flag given.
3. `5_validate_variant_genome.py` for the materialized genome, `6_validate_runtime_patch.py` for the
   sequence layer.
4. `applied + recorded-not-applied` in the manifest must equal `bcftools index -n <filtered vcf>`.
5. `tabix -l <filtered vcf>` must include `MT`. Missing MT means the DP cap was applied without the
   MT exemption.
6. Ti/Tv on the filtered VCF should be ~2.0 at GQ≥20. Materially lower means the filter did not apply.

## Gotchas

1. **Filter expressions naming `MT`/`X`/`Y` must come after the contig rename.** Tested against a
   `chr`-prefixed source they silently match nothing, which looks like a working filter that removes
   everything. Step 2 orders this correctly; keep it that way.
2. **`!(...)` is not valid in a bcftools filter expression.** It parses as an error, not negation.
   Write the De Morgan form: `(CHROM!="Y" || GT!="het")`.
3. **`GT!~"\."` does not exclude partial no-calls** — it matches every record. Use `GT="mis"` /
   `GT!="mis"`. (`GT="alt"` already excludes them, so the clause is belt-and-braces.)
4. **A hand-written contig rename map will bite you.** `chrEBV -> EBV` is the classic: the reference
   contig is `chrEBV`, so the rename produces a name `norm` cannot find. Step 2 derives the map from
   the FASTA index and only emits a line when the target name actually exists.
5. **`util/utils.py` is CRLF; `hdr.py` and `protoSpaceJAM.py` are LF.** An editor that normalizes
   line endings rewrites all 1,580 lines of utils.py and buries the real diff. Check with
   `file util/utils.py` after editing.
6. **`lru_cache(maxsize=2)` on `_load_contig_seq`** trades ~500 MB RSS for the saved unpickling.
   Raise it only if you also raise the `test_memory(4200)` gate in `protoSpaceJAM.py`.
