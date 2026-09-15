> **About this document.** This is the implementation spec for genotyping primer design (KASP and allele-specific
> PCR), copied verbatim from the approved design so that it lives next to the code. It was the design reference for
> the `genotyping` branch of gramene-swagger, milestones M0–M12; section numbers quoted elsewhere (§2.9, §5.3, …) refer
> to it. It is not updated as the code changes. **`docs/primer_design_api.md` is authoritative for behaviour**: where
> the two differ, the API document describes what is implemented. The known differences are listed in the final
> section, [Implementation deviations (M3–M7)](#implementation-deviations-m3m7).

# Genotyping primer design (KASP and allele-specific PCR): implementation spec

**Status.** Final implementation spec for the `genotyping` feature of gramene-swagger (primers API), written to be implemented without further design decisions. It is self-contained: no other design document needs to be read. Every finding raised in review is resolved, and §13 says where.

**Where this lands.**
- API: new git worktree `/usr/local/gramene/subsites/sorghum/v11/gramene-swagger-genotyping`, branch `genotyping`, created from `primer-design` (stacked PR later).
- UI: §8 is the front-end spec that a separate session implements in `gramene-primers` (and the `gramene-search` Primers tab).
- Nothing in any repository, pm2 process or port was changed while writing this. No check job was submitted.

**Evidence labels** used throughout:

| Label | Meaning |
|---|---|
| **[V]** | Verified while writing this spec, by running code, the dev API on port 50111 (read-only), `primer3_core` 2.6.1, `oligotm`/`ntthal`, or live Ensembl REST. Scripts and outputs: `genotyping-design/scratch-synth/` (Appendix A). |
| **[R]** | Verified in the research reports `genotyping-research/report-{backend,primer3,rules,frontend,pangenome_alleles}.md`. |
| **[I]** | Inferred, or a design decision that was not tested by running it. |
| **[U]** | An expectation that must be confirmed during implementation. |

---

## 0. Decisions at a glance

| # | Decision | Why | Rejected alternative |
|---|---|---|---|
| D1 | **Dedicated endpoints**: `GET /primers/variants`, `GET /primers/variants/{variant_id}`, `POST /primers/genotyping/design`. `POST /primers/check` gains an optional `genotyping` block; `GET /primers/genomes` gains variation flags. `POST /primers/design` is **unchanged**. | Swagger can require `variant` only in its own strict definition. A *set* has a different shape from a *pair*: putting indel ALT pairs into `pairs[]` breaks the two documented invariants `product_size = right.end − left.start + 1` and `len = end − start + 1`. `design.MODES`, the mode/enum sync test and the 22 existing request fixtures stay untouched. | `mode: "variant"` on `/primers/design` (widens two enums, overloads `pairs[]`, and leaks `'variant'` into `CheckRequest.mode`). |
| D2 | **Variants are served through swagger only**; the browser never calls Ensembl. | One origin, one error model. The server verifies REF against the local FASTA, normalizes the EVA/SAP insertion conventions, merges duplicate ids, flags EMS records and caches. gramene-search's own `ensemblRest` still points at release 108 [R]. | Browser-direct REST. |
| D3 | **One canonical identity**: `key = region:position:REF:ALT`, left-aligned VCF with an anchor base for indels (`1:11109:C:A`, `1:11282:CA:C`, `1:11502:C:CGT`). Designs and jobs key on it; Ensembl ids are for lookup and display only. | `tmp_1_11502_C_CGT` and `rs5413863549` are the same event and hash identically [V]. 88 of 3,976 real ids in one 50 kb window contain `,` or `*` or exceed 128 characters [V], so an id is a poor key. | Keying on the Ensembl id, or echoing raw Ensembl start/end/allele_string. |
| D4 | **Sets are first-class.** Every oligo states `role` (`as_ref`/`as_alt`/`common`), the allele it detects, its 3′ base, dye, `target_seq` (what anneals and the only sequence the check ever receives), `matched_seq`, `tail_seq`, `order_seq`. Each set ships order rows and a ready check fragment; the response ships a complete `check.request`. | The UI and scripts never derive roles, alleles, tails, oligo names or check bodies. | Returning pairs and letting each consumer assemble sets. |
| D5 | **Identity of a set is its content.** `set.key` = first 12 hex of `sha256(orientation\|as_ref.target_seq\|as_alt.target_seq\|common.target_seq)`; `set.id` = `S1…Sn` by rank; check pair ids `S1_REF`/`S1_ALT`. Results are matched to sets by the uppercase target-sequence triple, never by id or rank. | Re-designs renumber sets; matching by rank regresses the existing `matchCheckResults` behaviour. Keys verified: `f9df650ad116`, `1accc54c262d`, `7f9af6b1c938` [V]. | Rank-based ids. |
| D6 | **One Primer3 run per orientation per ladder level**, with `SEQUENCE_FORCE_LEFT_END`/`SEQUENCE_FORCE_RIGHT_END` for the allele-specific (AS) primer **and a `SEQUENCE_TARGET` guard** of 10 nt immediately outside the variant zone on the common-primer side. | Verified [V]: with a 50–120 bp product range and no guard, Primer3 returns a common primer at 1:11108–11131 footprint `TCTTTGACTAGCGAGAAATTCAGA` that covers the variant and overlaps the AS primer (the check would drop that product); with the guard the same run returns the clean 93 bp pair. The guard also protects the ALT haplotype because the zone includes both ALT 3′ anchor mappings. | Post-filtering only (misses multi-base shiftable indels), or a second Primer3 step with the AS primer fixed. |
| D7 | **ALT primers are derived, never re-designed**: same 5′ genomic anchor and length, 3′ base at the first haplotype difference. REF and common numbers come from the design run; the ALT and deliberate-mismatch primers are scored by a Primer3 `check_primers` run on their own haplotype (§4.10); only tailed structures and the mismatch duplex Tm come from `ntthal` with Primer3's salts. | Designing on each allele separately gives different primers [R]. `check_primers` on the REF haplotype reproduces the design run exactly [V], so REF, ALT and common numbers stay on one scale; `oligotm` is not used anywhere. | Two independent designs. |
| D8 | **No `PRIMER_PICK_ANYWAY`, ever.** The relaxation ladder stops at level 2 with hard floors: Tm ≥ 52 °C for every allele-specific primer and GC ≥ 15 % for the derived ones (the ALT primer and any deliberate-mismatch primer), as in §3.1 `as_min_tm`/`as_min_gc`. An orientation with nothing above the floor is reported `no_sets` with its `explain` data. | KASP cycles at 61→55 °C then 26 cycles at 55 °C [R]. `PICK_ANYWAY` returns 45.4 °C primers [R]; such a primer silently turns heterozygotes into wrong homozygous calls. At level 2 the same locus yields a real 52.87 °C primer [V]. | A `PICK_ANYWAY` rung with a ranking penalty. |
| D9 | **Allele calls are exact.** A genome is called `ref` or `alt` only when the observed core (the left-aligned VCF span — anchor base included — extended right by the shift, plus one base on each side; inserted bases included; §5.6) equals that haplotype's core exactly; anything else is `other`, and the observed core is always reported. | Lowest-edit-distance naming is wrong in repeats [V]: at rs5413864115 a genome that deleted **two** A's is 1 edit from ALT and 2 from REF, so distance would call it `alt`; `AAAA` would be called `ref`. Both are `other` under the exact rule [V]. | Nearest-haplotype naming with a divergence cap. |
| D10 | **Amplification prediction reads all three primers.** The common primer's status comes from **its own alignment at the called locus** — never from the other pair's product likelihood — taken as the best over the on-locus products of both pairs: a residual 3′ mismatch at distance 1 (or a blocked common site) means no amplification (genotype `no_call`), at distance 2–3 means `weak`. An off-locus amplifying product of the **other** allele's primer adds that dye and makes the genotype `both`; a same-allele paralog only raises `off_locus_products`. | The common primer is shared by both alleles: a private SNP under its 3′ end gives a lab no-call, which must not be reported as a concordant REF/ALT call. `classify.js` already yields `likely_weak` when *either* primer has a terminal mismatch, so the information is present. | Reading only the AS primer's `mm_pos`. |
| D11 | **Deliberate mismatch from the literal Little 1995 Table 9.8.1**, a 32-cell lookup keyed by each primer's own terminal pair and the template base opposite the mismatch; default position −2 (the position the table is defined for), −3 selectable and labelled extrapolated. | A class-rule shortcut differs from the table in 2 of 32 cells; for tmp_1_11193_C_T reverse at −2 the table gives REF A→C and ALT A→T, while the shortcut gives C for both [V]. | Class-rule shortcut, or Liu 2012 at −3 as the default (its pair notation is undefined). |
| D12 | **Ensembl never holds a design semaphore slot.** Resolution and neighbour fetching run before `semaphore.acquire`, inside the 45 s deadline, behind a process-wide limiter (4 concurrent, 5 s wait), single-flight coalescing, a circuit breaker, bounded LRU caches and byte caps. An outage degrades genotyping design to a warning where possible and never affects `/primers/design`. | `design()` acquires the semaphore before `runDesign`; a 10 s Ensembl wait inside a slot starves ordinary PCR designs (4 running, 16 waiting, 10 s wait → `503 BUSY`). | Resolving inside `runDesign`. |
| D13 | **Job ids and `ALGORITHM_VERSION` are untouched for existing requests.** `request.genotyping` is stored exactly as the client sent it (normalized values only, no server-added keys), and the genotyping algorithm version enters the job id through `jobId`'s `algo` argument as `'2+g1'`. | `canonicalJSON` skips undefined values, so plain jobs keep their ids; `job.request` stays a valid `PrimerCheckRequest`; a caller change invalidates only genotyping jobs. | Writing `algorithm`/`orientation` into the stored request (it then fails request validation). |
| D14 | **The server builds the order sheet** (rows, KASP mix, `[X/Y]` submission string with IUPAC-coded neighbours), and the UI state stays `v: 1` with an additive `genotyping` slice. | Every consumer exports identical, auditable rows; older UI builds ignore unknown state keys and fall back to their default mode [R]. | Client-side assembly; a `v: 2` bump (which discards saved views in old builds). |

### 0.1 What was verified while writing this spec

| # | Finding | Evidence [V] |
|---|---|---|
| F1 | Every primer sequence, coordinate and product size in every example of this document matches the reference genome, checked base by base through the dev API's `template_only` sequence for 1:10650–11960, 1:11150–11750 and 1:13600–13820. | `scratch-synth/verify_seqs.py` |
| F2 | Without the `SEQUENCE_TARGET` guard and with a 50 bp product minimum, Primer3 returns `PRIMER_RIGHT_0=423,24` (`TCTTTGACTAGCGAGAAATTCAGA`) against the forced left primer `372,30` — the common primer covers the variant and overlaps the AS primer. With `SEQUENCE_TARGET=402,10` the same record returns the 93 bp pair with common primer `464,24`. | `scratch-synth/p3run.py` runs E and F |
| F3 | The guard does not disturb the forced AS primer in either orientation, nor for a shiftable deletion: forward `375,29` + `SEQUENCE_TARGET=405,10` → common `462,22` (88 bp); reverse `425,25` + `SEQUENCE_TARGET=390,10` → common `297,28` (129 bp). | `p3run.py` runs G and H |
| F4 | Ladder level 2 (`min_tm 52`, `min_gc 20`, `max_tm_diff 6`, sizes to 32) returns `TCCGAATATAGTCATACTCTATTC` at 52.87 °C at the AT-rich forward orientation of 1:11109 with **no** `PICK_ANYWAY`. | `p3run.py` run I |
| F5 | The exact-core caller gives `ref`/`alt`/`other` correctly where distance-based calling fails: rs5413864115 core `TCAAAG`/`TCAAG`; a two-A deletion (`TCAG`) and an `AAAA` insertion (`TCAAAAG`) are `other`, although their edit distances name `alt` and `ref`. Same for the TGG repeat rs5413863413 (core `ATTGGTGGTGGTGGTGGTA` / `ATTGGTGGTGGTGGTA`). | `scratch-synth/caller.py` |
| F6 | Guard zones computed from both discriminating positions and both ALT anchor mappings: rs871475760 `[11109,11109]`, tmp_1_11193_C_T `[11193,11193]`, rs5413864115 `[11282,11286]`, tmp_1_11502_C_CGT `[11502,11503]`, rs5413863413 `[13735,13752]`. The last one is the case a REF-only guard misses. | `scratch-synth/zone.py` |
| F7 | The literal Table 9.8.1 lookup gives: rs871475760 −2 reverse REF and ALT A→G, forward REF and ALT T→C; tmp_1_11193_C_T −2 reverse REF A→C and ALT A→T, forward both T→G; at −3 reverse both T→G and forward both C→A. | `scratch-synth/little.py` |
| F8 | Every Tm and structure value quoted in the examples comes from Primer3 (`check_primers` for derived oligos) or from `ntthal` with `-mv 50 -dv 1.5 -n 0.6 -d 50 -t 37 -r`: Tm 57.10, 56.51, 58.72, 57.28, 56.34, 57.08, 55.21, 56.43, 59.22, 56.09, 59.19, 60.43, 60.02; tailed hairpins FAM+S1 42.62, HEX+S1 51.50, FAM+S2 36.16, HEX+S2 66.22; tails alone 0; cross-dimers `ANY` 6.89 and 5.46, `END1` 8.21 (the maximum of both directions, §4.14); −2 mismatch duplex 54.69 against 57.04 unmodified. Negative values occur (−11.46) and must be clamped to 0. | `scratch-synth` thermo run; re-run by `scratch-fix/score_sets.js` |
| F9 | Cost, computed with the server's own `check/cost.js` and the dev API's real genome sizes (reference 708,735,318; 119 others summing 83,774,241,795): 6 primers full panel 2,666 → **2,690** with genotyping; 12 → 5,356; 13 → **5,800** (largest accepted); 14 → 6,245 (over the 6,000 limit); 6 primers over 3 genomes → 95. | `scratch-synth` cost run |
| F10 | Real Ensembl ids in 1:11080–61079: 3,978 records, 88 contain `,`, 46 contain `*`, 11 exceed 128 characters, longest 255. The pattern `^[A-Za-z0-9][A-Za-z0-9_.:,*-]{0,254}$` accepts all of them; narrower patterns reject 53 and 88. | `scratch-synth` id-pattern run |
| F11 | Live Ensembl neighbours used in the examples: `rs873774986` C/T at 1:11059 (EVA), `tmp_1_11069_G_A` (EMS) at 11069, `rs872438201` A/T at 11161, `rs873026643` A/G at 11182, `rs5413863234` G/A at 11500 (EVA), `rs5413864238` T/C at 11509. Both insertion conventions confirmed: `tmp_1_11502_C_CGT` start 11503 end 11502 (SAP), `rs5413863549` start = end = 11503 (EVA). | live overlap calls |
| F12 | The four `submission_sequence` strings of §4.17 were regenerated from the genome plus live Ensembl neighbours and match this document character for character. | `scratch-synth` submission run |

---

## 1. User flows

The designer gains a fifth mode tab, **Genotyping**, beside Gene, Transcript, Region and Sequence. Each step below names the client call (§8.1) and the response fields it uses (§2).

### 1.1 Flow A — pick a variant from Ensembl (sorghum_bicolor)

| # | The user sees or does | Client call | Response fields used | State written (`state.genotyping`, §8.3) |
|---|---|---|---|---|
| A1 | Opens the Primers tab of SORBI_3001G000200 and chooses **Genotyping**. | `listGenomes('sorghum_bicolor')` (memoized) | `variation.available: true`, `variation.release: "115"`, `genomes[].has_variation` | `mode: 'genotyping'` |
| A2 | The variant picker opens on the gene span ± 2 kb (1:9180–16899), editable up to 50 kb. | `listVariants({system_name, region:'1', start:9180, end:16899})` | `variants[]`: `key`, `label`, `kind`, `ids`, `records[].source`, `ems`, `consequence`, `designable`, `issues`, `multiallelic` | `window` |
| A3 | Filters (SNV only, hide EMS), types `rs871475760`, clicks the row `1:11109 C/A`. | none | `key`, `vcf`, `ids[0]` | `variantKey: '1:11109:C:A'`, `variantId`, `alt` |
| A4 | Chooses **KASP**. Defaults: REF = FAM, ALT = HEX, both orientations, 6 sets, avoid repeats off (§4.5; a default to confirm, §11 item 6). | none | none | `assay` (changed keys only) |
| A5 | Presses **Design**. | `designGenotyping({system_name, variant:{region:'1',position:11109,ref:'C',alt:'A'}, assay:{type:'kasp'}})` | `variant`, `template`, `orientations`, `sets[]`, `check.request`, `warnings` | `designed: true` |
| A6 | Sees the template map (variant marker, neighbour ticks, set lanes) and the sets table: S1 reverse, 65 bp, `usable`, chips `ALT_PRIMER_SUBOPTIMAL` and `TAILED_STRUCTURE`; S2 forward, 92 bp, `poor`, chips `TAILED_STRUCTURE` (high), `ALT_PRIMER_SUBOPTIMAL` and `NEIGHBOUR_IN_PRIMER`. | none | `sets[].primers.*`, `sets[].warnings`, `sets[].quality` | `selectedSetKey` |
| A7 | Ticks S1 and S2. The cost line shows 6 distinct primers, ≈ 2,690 CPU-s for the full 120-genome panel (limit 6,000). | `estimateCheckCpu({…, genotyping: true})` | `sets[].check.pairs` | `checkedSetKeys` |
| A8 | Presses **Check** (specificity + pan-genome). | `runCheck(design.check.request)` or `buildGenotypingCheckRequest(...)` | `job_id`, `progress`, `results.genotyping` | `check.jobId`, `check.submitted` |
| A9 | Opens the **Alleles** tab: the reference control first, then each assembly's allele and each set's predicted signal, with a "disagreements only" filter. | polling (existing) | `results.genotyping.genomes[]`, `.sets[].genomes[]`, `.summary` | `view.tab: 'alleles'` |
| A10 | **Export → Order sheet (TSV)**: 6 rows (REF-FAM, ALT-HEX, common for S1 and S2), the KASP mix note and the `[C/A]` submission string. | none (client export) | `sets[].order[]`, `variant.submission_sequence`, `assay.kasp_mix` | none |

**Why the older tabs look "wrong" for ALT pairs.** On a REF genome the ALT pair has a 3′-terminal mismatch, so the check calls it `likely_weak` / `single_mismatch` [R]. The Specificity and Pan-genome tabs are unchanged; each shows a one-line note pointing at the Alleles tab. Allele questions are answered only by `results.genotyping` (§5.4).

### 1.2 Flow B — manual entry on any genome with sequence

| # | Step | Call | Notes |
|---|---|---|---|
| B1 | On a genome without variation data (e.g. `sorghum_rio`) the picker is hidden and manual inputs appear: region, position, ref, alt. | `listGenomes` → `has_variation: false` | Alleles are A/C/G/T strings, or `-` for an empty allele (§2.1). |
| B2 | **Check reference** | `designGenotyping({…, variant:{region,position,ref,alt}, template_only:true})` | 200 → `variant.ref_verified: true`, `template`, `variant.discriminating`; 400 `REF_MISMATCH` → "the genome has T at 1:…". |
| B3 | **Design** | same without `template_only` | Response warning `NO_VARIATION_DATA`, `neighbours.data: "none"`. |
| B4 | **Check** | `runCheck` with `system_name: 'sorghum_rio'` | The pan-genome check covers the other assemblies; `sorghum_rio` is the reference control. |

The same flow serves an EMS mutation that is not in Ensembl, on sorghum_bicolor; there, known neighbours are still screened.

### 1.3 Flow C — deep link by id

A host passes `state={v:1, mode:'genotyping', genotyping:{variantId:'tmp_1_11502_C_CGT'}}`. The component calls `getVariant('tmp_1_11502_C_CGT', {system_name})` and shows a card reading `1:11502^11503 -/GT`, ids `tmp_1_11502_C_CGT` and `rs5413863549`, with the note "2 Ensembl ids describe the same event". A multi-allelic id returns one entry per alt and the card asks which alt to use. The picker then designs from the entry's `vcf`, not from the id (D3), so ids containing `,` or `*` never reach a design request.

### 1.4 Flow D — gel AS-PCR with a deliberate mismatch

The user picks **Gel AS-PCR**: the deliberate mismatch defaults to on at −2, tails to none, products 150–300 bp.
- The sets table shows each primer's mismatch, e.g. "−2 A→G (weak extra; terminal mismatch maximum) — Little 1995 Table 9.8.1".
- The order sheet has no dyes.
- In the check the wrong-allele primer is predicted `blocked` (`mm_pos [2,1]`, `unlikely`), not `terminal_mismatch` [V].

### 1.5 Flow E — genomes and variants without Ensembl data

If `variation.available` is false for the chosen genome, the picker never appears; if it is true but Ensembl is unreachable, the picker shows a retry banner with a countdown and an **Enter manually** button (§1.6). A manual design still succeeds, with `neighbours.data: "unavailable"` and warning `NEIGHBOURS_UNAVAILABLE`, so the user is told that nothing was screened.

### 1.6 Flow F — degraded and error paths (each is a defined UI state)

| Situation | API answer | UI behaviour |
|---|---|---|
| Ensembl down during the list call | 503 `VARIATION_SOURCE_UNAVAILABLE`, `details.retry_after_s: 30` | Banner with countdown, Retry, and **Enter manually** prefilled from the last selected entry. |
| Ensembl down during a design by id | 503 `VARIATION_SOURCE_UNAVAILABLE` | Same banner; the manual form is prefilled from the variant card. |
| Ensembl slow, outbound queue full | 503 `VARIATION_SOURCE_UNAVAILABLE`, `retry_after_s: 5` | Same banner, shorter countdown. Ordinary PCR designs are unaffected (D12). |
| Ensembl down during a manual design on sorghum_bicolor | 200 + warning `NEIGHBOURS_UNAVAILABLE` | Chip: "known variants near the primers were not screened". |
| REF does not match the genome | 400 `REF_MISMATCH` `{region, position, given, genome}` | Inline error on the ref field. |
| Unknown id | 404 `UNKNOWN_VARIANT` | Inline error on the id search. |
| Genome has no variation data but an id was sent | 422 `NO_VARIATION_DATA` | Switch to manual inputs. |
| One orientation blocked (tmp_1_11502_C_CGT forward) | 200, `orientations.forward.status: "blocked"` with `blockers[]` | "forward blocked by rs5413863234 (G/A) 4 nt from the 3′ end". |
| No sets at all | 200, `sets: []`, warning `NO_SETS`, per-orientation `attempts[].explain` | ExplainPanel per orientation and level, with a hint from the dominant rejection reason. |
| Relaxed constraints | 200, `sets[].relaxation_level ≥ 1`, warning `RELAXED_CONSTRAINTS` | `L1`/`L2` badge listing `settings.ladder[level].changes`. |
| Check too large | client estimate > 6,000; the server would answer 422 `JOB_TOO_LARGE` | Submit disabled, estimate shown, "check fewer sets or genomes". |
| Results expired | 404 `UNKNOWN_JOB` | Existing "Results expired — Re-run". |

### 1.7 Flow G — the order sheet

Rows come from `sets[].order[]` of the checked sets (or the selected set when none is checked), three per set in the order REF, ALT, common:

| name | set | role | allele | dye | sequence (5′→3′) | target-specific part | tail | length | Tm target (°C) | GC target (%) | orientation | product REF/ALT (bp) | variant |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| rs871475760_S1_REF_FAM | S1 | as_ref | C | FAM | GAAGGTGACCAAGTTCATGCTATCTTTGACTAGCGAGAAATTCAG | ATCTTTGACTAGCGAGAAATTCAG | GAAGGTGACCAAGTTCATGCT | 45 | 57.10 | 37.50 | reverse | 65/65 | 1:11109:C:A |
| rs871475760_S1_ALT_HEX | S1 | as_alt | A | HEX | GAAGGTCGGAGTCAACGGATTATCTTTGACTAGCGAGAAATTCAT | ATCTTTGACTAGCGAGAAATTCAT | GAAGGTCGGAGTCAACGGATT | 45 | 56.51 | 33.33 | reverse | 65/65 | 1:11109:C:A |
| rs871475760_S1_COM | S1 | common | | | AGCTTCTCTAAGTGGTTATCCGA | AGCTTCTCTAAGTGGTTATCCGA | | 23 | 58.72 | 43.48 | reverse | 65/65 | 1:11109:C:A |

All values are real [V]. Below the rows, `#`-prefixed note lines carry:
- **KASP primer mix** — 12 µL REF-FAM + 12 µL ALT-HEX + 30 µL common, all at 100 µM, made up to 100 µL with 46 µL water (ratio from Makhoul et al. 2020 [R]; the water volume is inferred [I]).
- **Service submission string** — `YCCTCAAAAAGCTTCTCTAAGTGGTTATCCGAATATAGTCATACTCTATT[C/A]TGAATTTCTCGCTAGTCAAAGATAACAAAAATAGCATATTCTGGATTTCT` [V]; the leading `Y` codes neighbour rs873774986 (C/T) at 1:11059.

Every cell passes through the existing formula-injection guard `tsvCell` (`exporters.ts:24-32`).
---

## 2. API contract

This section is normative. It extends the existing `/primers` conventions (`docs/primer_design_api.md`): JSON bodies with `Content-Type: application/json`, `Cache-Control: no-store`, CORS `*`, handler errors `{message, code, details}`, validator errors `{message, errors[]}`, every 503 carrying `details.retry_after_s` and a `Retry-After` header. Clients omit optional fields rather than sending `null`.

| Endpoint | Purpose | Sync |
|---|---|---|
| `GET /primers/variants` | known variants in a window, normalized, REF-verified, duplicate ids merged | yes, outside the design semaphore |
| `GET /primers/variants/{variant_id}` | resolve one Ensembl variation id | yes, outside the design semaphore |
| `POST /primers/genotyping/design` | design KASP or AS-PCR sets for one variant | yes, 45 s deadline; the semaphore covers only Primer3/thermo work |
| `POST /primers/check` (**changed**) | optional `genotyping` block | job |
| `GET /primers/check/{job_id}` (**changed**) | `results.genotyping` | job |
| `GET /primers/genomes` (**changed**) | `variation` and `genomes[].has_variation` | yes |

The three new paths go into `swagger.yaml` after `/primers/check/{job_id}` (which ends at line 1041) and **before** `/{collection}` (line 1043). Both new GET operations declare `consumes: [application/json]`, because `contract.test.js:107` asserts it for every primers operation and the existing GETs declare it (`swagger.yaml:931-932`, `1014-1015`) [V].

### 2.1 Conventions

**Coordinates.** Genomic coordinates are 1-based inclusive on the reference plus strand. The genotyping template is always the plus strand (`template.strand: 1`); template index `t = g − template.start + 1`. `template.alt_seq` is the ALT haplotype: its indices equal REF indices up to the variant and are shifted by `delta = len(vcf.alt) − len(vcf.ref)` after it.

**Variant representations.** Every variant object carries all of these. The values below are real [V]:

| Field | Meaning | rs871475760 | rs5413864115 | tmp_1_11502_C_CGT |
|---|---|---|---|---|
| `key` | `region:position:REF:ALT`, left-aligned VCF with an anchor base for indels. **Canonical identity**; designs and job hashes use it. | `1:11109:C:A` | `1:11282:CA:C` | `1:11502:C:CGT` |
| `vcf` | `{position, ref, alt}` of the key. This is what the check takes. | `{11109,"C","A"}` | `{11282,"CA","C"}` | `{11502,"C","CGT"}` |
| `minimal` | Ensembl style `{start, end, ref, alt}`; `-` = empty allele. Insertion: `start = end + 1`. | `{11109,11109,"C","A"}` | `{11283,11283,"A","-"}` | `{11503,11502,"-","GT"}` |
| `label` | human text | `1:11109 C/A` | `1:11283-11283 A/-` | `1:11502^11503 -/GT` |
| `shift` | how many bases the event can slide right (0 = not shiftable) | 0 | 2 | 0 |
| `zone` | `{start, end}`: the guarded span — the VCF span, both discriminating positions and both ALT 3′ anchor mappings (§4.4). Nothing but an allele-specific primer may touch it. | `{11109,11109}` | `{11282,11286}` | `{11502,11503}` |
| `discriminating.forward` | `{position, ref_base, alt_base, alt_maps_to}`: the 3′ base of a **forward** AS primer — the first base, left to right, where the haplotypes differ. `alt_maps_to` is the REF coordinate of the ALT base there, `null` when it is an inserted base. | `{11109,"C","A",11109}` | `{11285,"A","G",11286}` | `{11503,"A","G",null}` |
| `discriminating.reverse` | the same, scanning right to left (plus-strand bases) | `{11109,"C","A",11109}` | `{11283,"A","C",11282}` | `{11502,"C","T",null}` |

For the five-TGG-repeat deletion `rs5413863413` (`1:13735:TTGG:T`), `shift` is 13, forward is `{13749,"G","A",13752}`, reverse is `{13738,"G","T",13735}` and `zone` is `{13735,13752}` [V]. A guard built only from the REF span would leave 13752 unprotected, which is exactly how a common primer ends up on top of the ALT allele-specific primer.

**Manual allele input** (`variant.ref`/`variant.alt` in a design request), case-insensitive, at most 50 nt per allele:

| Style | When | `position` means | Example (same event) |
|---|---|---|---|
| VCF | neither allele is `-` | first base of `ref` | `{position: 11282, ref: "CA", alt: "C"}` |
| Ensembl | exactly one allele is `-` | deletion: first deleted base; insertion: the base **after** the insertion point (`minimal.start`) | `{position: 11283, ref: "A", alt: "-"}` |

**Orientation.**

| `orientation` | AS primer | 3′ end at | common primer |
|---|---|---|---|
| `forward` | left primer, plus strand | `discriminating.forward.position` | right primer, downstream of `zone.end` |
| `reverse` | right primer, minus strand | `discriminating.reverse.position` | left primer, upstream of `zone.start` |

**Oligo sequence fields.** Every oligo carries all four, uppercase 5′→3′:

| Field | Meaning | Sent to `/primers/check`? |
|---|---|---|
| `target_seq` | the part that anneals, exactly as ordered, including any deliberate mismatch | **yes** — the only sequence ever sent |
| `matched_seq` | the perfect-match sequence before any deliberate mismatch | no |
| `tail_seq` | 5′ KASP tail (FAM or HEX), or `null` | no |
| `order_seq` | `tail_seq + target_seq`: what the vendor synthesizes | no |

**Genomic footprint.** `genomic = {region, start, end, strand, blocks[]}`; `strand` is the strand whose 5′→3′ sequence is `matched_seq`; `blocks` are ascending reference ranges. An ALT primer across a deletion has two blocks (e.g. `[{11257,11282},{11284,11286}]` [V]); inserted bases have no reference coordinate and are counted in `inserted_bases`.

**Thermodynamics.** °C with 2 decimals, rounded by the rounding rule below (so Primer3's `58.025` is reported as `58.03`). Salts are the effective params, defaulting to Primer3's 50 mM monovalent, 1.5 mM divalent, 0.6 mM dNTP, 50 nM oligo. `*_th` values are structure Tm at 37 °C. `tm_method` is `primer3` for a perfectly matched oligo — the Tm of the design run for `as_ref` and `common`, of the `check_primers` run on the ALT haplotype for `as_alt` (§4.10) — or `ntthal_duplex` for a deliberate-mismatch primer (duplex against its own-allele template). Negative `ntthal` results are clamped to 0 [V]. `gc` is computed from the oligo's own sequence (100 × GC ÷ length, same rounding), not parsed from Primer3, so it never depends on a tool's print precision.

**Rounding — one rule for every reported number.** Every input is either an exact decimal as the tool prints it (Primer3 and `ntthal` print at most 6 decimal places) or an exact ratio of integers (GC %, the product term of §4.16). A derived number — a difference, a maximum, a distance outside a range, a sum such as `score` — is computed exactly from the unrounded inputs in decimal arithmetic (for example as integers in units of 10⁻⁶), never in binary floating point; a ratio is rounded directly from its integer numerator and denominator. The exact value is then rounded **once, half away from zero**: to 2 decimals for every number in a response, except `primer3_penalty` (4 decimals) and the hairpin value quoted in an order-row note (1 decimal, §4.17). Threshold tests (the §4.6 floors, the §2.15 severities, the §4.13 ranges) and the §4.16 ordering compare exact values, never rounded ones.

Worked case [V]: tmp_1_11193_C_T S1 `common_minus_as` = 57.102 − 56.427 = 0.675 exactly, reported **0.68**. In binary floating point the difference is 0.6749999999999972, which `Math.round(x × 100) / 100` turns into the wrong 0.67. Over the 65 candidates scored for the §2.9–§2.10 examples, that float shortcut changes 7 values (`scratch-fix2/score2.js`); §7.3 pins the rule.

**Discrimination.** `mm_pos` lists mismatch distances from the 3′ end (1 = terminal), exactly as in check results. `likelihood` is the check's own class (`likely`, `likely_weak`, `unlikely`) from `classify.js:71-78`, computed with the check's default params. `terminal_mismatch_class` is the Little 1995 class of the primer's 3′ base against the other allele's template base (`max`, `strong`, `medium`, `weak`) — for a KASP set with no deliberate mismatch this is the only in-silico discrimination signal the literature supports [R].

**Identifiers.**

| Id | Rule | Example |
|---|---|---|
| `set.id` | `S` + (rank + 1) | `S1` |
| `set.key` | first 12 hex of `sha256(orientation + "\|" + as_ref.target_seq + "\|" + as_alt.target_seq + "\|" + common.target_seq)` | `f9df650ad116` [V] |
| check pair ids | `{set.id}_REF`, `{set.id}_ALT` | `S1_REF` |
| order row `name` | `{label}_{set.id}_{REF\|ALT\|COM}`, plus `_FAM`/`_HEX` for tailed AS primers | `rs871475760_S1_REF_FAM` |
| `label` | request `label`, else the first id, else `key`; then `[^A-Za-z0-9_.-] → _`, truncated to 40 characters | `rs871475760`, `1_11109_C_A` |

**Enumerations.**

| Name | Values |
|---|---|
| variant `kind` | `snv`, `mnv`, `insertion`, `deletion`, `complex` |
| assay `type` | `kasp`, `as_pcr` |
| request `orientation` | `both`, `forward`, `reverse` |
| `tails` | `none`, `ref_fam_alt_hex`, `ref_hex_alt_fam` |
| `deliberate_mismatch` | `none`, `auto` |
| `mismatch_position` | `2`, `3` |
| `neighbour_policy` | `avoid_3p`, `ignore` |
| oligo `role` | `as_ref`, `as_alt`, `common` |
| orientation `status` | `ok`, `no_sets`, `blocked`, `skipped` |
| set `quality` | `good`, `usable`, `poor` |
| genotype `allele` | `ref`, `alt`, `other`, `ambiguous`, `missing`, `unavailable` |
| primer `status` (check) | `match`, `weak`, `terminal_mismatch`, `uncertain`, `blocked`, `no_product`, `unknown` |
| set `predicted` (check) | `ref`, `alt`, `both`, `none`, `no_call`, `unknown` |

### 2.2 `GET /primers/genomes` (additive)

Two additions, nothing removed: a top-level `variation` object for the query genome and a `has_variation` flag on every genome. The existing entry keys (`system_name`, `display_name`, `taxon_id`, `map_id`, `is_query`, `has_sequence`, `has_blastdb`, `has_cdna_blastdb`, `repeat_masking`, `total_bases`, `warnings`) are unchanged [V].

```jsonc
{ "system_name": "sorghum_bicolor",
  "species": { "taxon_id": 4558, "name": "Sorghum bicolor" },
  "variation": { "available": true, "source": "ensembl", "release": "115" },
  "counts": { "total": 120, "with_blastdb": 120, "with_cdna_blastdb": 120 },
  "genomes": [
    { "system_name": "sorghum_bicolor", "display_name": "Sb bicolor BTx623 v3", "is_query": true, "has_sequence": true,
      "has_blastdb": true, "has_cdna_blastdb": true, "has_variation": true, "repeat_masking": "unmasked_copy",
      "total_bases": 708735318, "taxon_id": 4558006, "map_id": "GCA_000003195.3", "warnings": [] },
    { "system_name": "sorghum_rio", "display_name": "Sb bicolor PI651496 Rio", "has_variation": false, "…": "…" } ] }
```

`variation` for a genome without data is `{"available": false, "source": null, "release": null}` — never `null` itself.

Add to `PrimerGenomesResponse.properties` (`swagger.yaml:2310`):

```yaml
variation:
  type: object
  properties:
    available: { type: boolean }
    source: { type: string, enum: [ensembl], x-nullable: true }
    release: { type: string, x-nullable: true }
```

Add to `PrimerGenomesResponse.properties.genomes.items.properties`:

```yaml
has_variation:
  type: boolean
  description: "known-variant data (Ensembl) is available for this genome; only sorghum_bicolor on sorghum_v11"
```

### 2.3 `GET /primers/variants`

Lists the known variants of a window on a genome that has variation data. One entry per designable biallelic `key`; a multi-allelic record yields one entry per designable alt, and records with the same key merge.

| Query param | Type | Required | Rule |
|---|---|---|---|
| `system_name` | `^[a-z0-9_]+$`, ≤ 128 | yes | in the catalog and has variation data, else `422 NO_VARIATION_DATA` |
| `region` | string ≤ 255 | yes | must be a sequence of the assembly (`404 UNKNOWN_REGION`) |
| `start`, `end` | integer ≥ 1 | yes | `start ≤ end ≤ region length` (`400 REGION_OUT_OF_BOUNDS`); `end − start + 1 ≤ 50,000` (`400 VARIANT_WINDOW_TOO_LONG`) |
| `types` | CSV of `snv,mnv,insertion,deletion,complex` | no | default all |
| `include_ems` | boolean | no | default `true`; `false` drops entries all of whose records are EMS |
| `limit` | integer 1–5000 | no | default 2000; `truncated: true` when more matched |

**Semantics.** An entry is in the window when its `minimal` span overlaps `[start, end]`; an insertion is in when `minimal.end` or `minimal.start` is. Entries are sorted by `vcf.position`, then `key`. `synonyms` is always `[]` here (only the id lookup fills it). Entries with `designable: false` are still returned, with `issues[]`, so the picker can explain why. **No entry is ever dropped because of its id** (§3.4).

**Example [V]** — `GET /primers/variants?system_name=sorghum_bicolor&region=1&start=11180&end=11290`:

```json
{
  "system_name": "sorghum_bicolor",
  "region": "1",
  "start": 11180,
  "end": 11290,
  "source": { "name": "ensembl", "release": "115" },
  "total": 4,
  "returned": 4,
  "truncated": false,
  "variants": [
    { "key": "1:11182:A:G", "ids": ["rs873026643"], "synonyms": [], "label": "1:11182 A/G", "kind": "snv", "region": "1",
      "vcf": { "position": 11182, "ref": "A", "alt": "G" }, "minimal": { "start": 11182, "end": 11182, "ref": "A", "alt": "G" },
      "alleles": ["A", "G"], "multiallelic": null, "shift": 0, "zone": { "start": 11182, "end": 11182 },
      "discriminating": { "forward": { "position": 11182, "ref_base": "A", "alt_base": "G", "alt_maps_to": 11182 },
                          "reverse": { "position": 11182, "ref_base": "A", "alt_base": "G", "alt_maps_to": 11182 } },
      "records": [{ "id": "rs873026643", "source": "EVA", "ems": false }], "ems": false,
      "consequence": "3_prime_UTR_variant", "ref_verified": true, "designable": true, "issues": [] },
    { "key": "1:11193:C:T", "ids": ["tmp_1_11193_C_T"], "synonyms": [], "label": "1:11193 C/T", "kind": "snv", "region": "1",
      "vcf": { "position": 11193, "ref": "C", "alt": "T" }, "minimal": { "start": 11193, "end": 11193, "ref": "C", "alt": "T" },
      "alleles": ["C", "T"], "multiallelic": null, "shift": 0, "zone": { "start": 11193, "end": 11193 },
      "discriminating": { "forward": { "position": 11193, "ref_base": "C", "alt_base": "T", "alt_maps_to": 11193 },
                          "reverse": { "position": 11193, "ref_base": "C", "alt_base": "T", "alt_maps_to": 11193 } },
      "records": [{ "id": "tmp_1_11193_C_T", "source": "EMS_PMID38100514_Jiao", "ems": true }], "ems": true,
      "consequence": "3_prime_UTR_variant", "ref_verified": true, "designable": true, "issues": [] },
    { "key": "1:11203:C:T", "ids": ["tmp_1_11203_C_T"], "synonyms": [], "label": "1:11203 C/T", "kind": "snv", "region": "1",
      "vcf": { "position": 11203, "ref": "C", "alt": "T" }, "minimal": { "start": 11203, "end": 11203, "ref": "C", "alt": "T" },
      "alleles": ["C", "T"], "multiallelic": null, "shift": 0, "zone": { "start": 11203, "end": 11203 },
      "discriminating": { "forward": { "position": 11203, "ref_base": "C", "alt_base": "T", "alt_maps_to": 11203 },
                          "reverse": { "position": 11203, "ref_base": "C", "alt_base": "T", "alt_maps_to": 11203 } },
      "records": [{ "id": "tmp_1_11203_C_T", "source": "EMS_PMID29378822_Addo-Qu", "ems": true }], "ems": true,
      "consequence": "3_prime_UTR_variant", "ref_verified": true, "designable": true, "issues": [] },
    { "key": "1:11282:CA:C", "ids": ["rs5413864115"], "synonyms": [], "label": "1:11283-11283 A/-", "kind": "deletion", "region": "1",
      "vcf": { "position": 11282, "ref": "CA", "alt": "C" }, "minimal": { "start": 11283, "end": 11283, "ref": "A", "alt": "-" },
      "alleles": ["A", "-"], "multiallelic": null, "shift": 2, "zone": { "start": 11282, "end": 11286 },
      "discriminating": { "forward": { "position": 11285, "ref_base": "A", "alt_base": "G", "alt_maps_to": 11286 },
                          "reverse": { "position": 11283, "ref_base": "A", "alt_base": "C", "alt_maps_to": 11282 } },
      "records": [{ "id": "rs5413864115", "source": "EVA", "ems": false }], "ems": false,
      "consequence": "3_prime_UTR_variant", "ref_verified": true, "designable": true, "issues": [] }
  ],
  "warnings": []
}
```

**Multi-allelic [V].** `rs5413863494` (C/T/G) yields `1:10718:C:T` with `multiallelic: {"alleles":["C","T","G"],"other_alts":["G"]}` and `1:10718:C:G` with `other_alts: ["T"]`. `rs5413863901` (C/T/`*`) yields one entry `1:11318:C:T` with `other_alts: ["*"]` and issue `STAR_ALLELE`; the BAP record `tmp_1_11318_C_T` has the same key and merges into it.

**Errors:** validator 400; `400 VARIANT_WINDOW_TOO_LONG {length, max}`; `400 REGION_OUT_OF_BOUNDS`; `404 UNKNOWN_GENOME`; `404 UNKNOWN_REGION`; `422 NO_VARIATION_DATA {system_name}`; `422 NO_SEQUENCE`; `503 VARIATION_SOURCE_UNAVAILABLE`; `503 FEATURE_DISABLED`.
**Warnings:** `VARIATION_RECORDS_SKIPPED {count, reasons}`, `REF_MISMATCHES {count}`.

### 2.4 `GET /primers/variants/{variant_id}?system_name=`

Resolves one Ensembl variation id and returns one entry per designable alt, with every merged alias. The path parameter matches `^[A-Za-z0-9][A-Za-z0-9_.:,*-]{0,254}$` and must be URL-encoded by the client; the server builds the outbound URL with `encodeURIComponent` (§3.3). This pattern accepts all 3,978 ids of the 50 kb research window, including `tmp_1_13549_TTA_T,*` and a 255-character id; the narrower patterns proposed elsewhere reject 53 and 88 of them [V].

**Example [V]** — `GET /primers/variants/tmp_1_11502_C_CGT?system_name=sorghum_bicolor`:

```json
{
  "requested_id": "tmp_1_11502_C_CGT",
  "system_name": "sorghum_bicolor",
  "source": { "name": "ensembl", "release": "115" },
  "variants": [
    { "key": "1:11502:C:CGT", "ids": ["tmp_1_11502_C_CGT", "rs5413863549"], "synonyms": [], "label": "1:11502^11503 -/GT",
      "kind": "insertion", "region": "1",
      "vcf": { "position": 11502, "ref": "C", "alt": "CGT" }, "minimal": { "start": 11503, "end": 11502, "ref": "-", "alt": "GT" },
      "alleles": ["-", "GT"], "multiallelic": null, "shift": 0, "zone": { "start": 11502, "end": 11503 },
      "discriminating": { "forward": { "position": 11503, "ref_base": "A", "alt_base": "G", "alt_maps_to": null },
                          "reverse": { "position": 11502, "ref_base": "C", "alt_base": "T", "alt_maps_to": null } },
      "records": [ { "id": "tmp_1_11502_C_CGT", "source": "SAP_PMID35653240_Boatwri", "ems": false },
                   { "id": "rs5413863549", "source": "EVA", "ems": false } ],
      "ems": false, "consequence": "3_prime_UTR_variant", "ref_verified": true, "designable": true, "issues": [] }
  ],
  "warnings": [
    { "code": "DUPLICATE_VARIANT_IDS", "message": "2 Ensembl ids describe the same event 1:11502:C:CGT; they were merged",
      "details": { "key": "1:11502:C:CGT", "ids": ["tmp_1_11502_C_CGT", "rs5413863549"] } }
  ]
}
```

`GET /primers/variants/rs871475760?system_name=sorghum_bicolor` returns one entry, key `1:11109:C:A`, with `synonyms: ["tmp_1_11109_C_A"]`; Ensembl's `"."` synonym is dropped [V].

**Errors:** validator 400; `404 UNKNOWN_VARIANT {id, system_name}`; `404 UNKNOWN_GENOME`; `422 NO_VARIATION_DATA`; `422 AMBIGUOUS_VARIANT_MAPPING {id, mappings}`; `422 VARIANT_NOT_ON_ASSEMBLY {id}`; `503 VARIATION_SOURCE_UNAVAILABLE`; `503 FEATURE_DISABLED`.

### 2.5 Swagger 2.0: the three new paths

```yaml
  /primers/variants:
    x-swagger-router-controller: primers
    get:
      summary: "Known variants (Ensembl REST) in a genomic window, normalized to left-aligned VCF keys and verified against this site's genome"
      description: "One entry per designable biallelic key; Ensembl records describing the same event are merged. Only genomes with variation data (GET /primers/genomes has_variation) can be queried. Cache-Control: no-store."
      operationId: listPrimerVariants
      tags: [ "Primer design" ]
      consumes: [ application/json ]
      produces: [ application/json ]
      parameters:
        - { in: query, name: system_name, required: true, type: string, pattern: "^[a-z0-9_]+$", maxLength: 128 }
        - { in: query, name: region, required: true, type: string, maxLength: 255 }
        - { in: query, name: start, required: true, type: integer, minimum: 1 }
        - { in: query, name: end, required: true, type: integer, minimum: 1 }
        - { in: query, name: types, required: false, type: array, collectionFormat: csv, items: { type: string, enum: [snv, mnv, insertion, deletion, complex] } }
        - { in: query, name: include_ems, required: false, type: boolean }
        - { in: query, name: limit, required: false, type: integer, minimum: 1, maximum: 5000 }
      responses:
        200: { description: "variants in the window", schema: { $ref: "#/definitions/PrimerVariantList" } }
        400: { description: "validator errors, VARIANT_WINDOW_TOO_LONG, REGION_OUT_OF_BOUNDS", schema: { $ref: "#/definitions/PrimerError" } }
        404: { description: "UNKNOWN_GENOME, UNKNOWN_REGION", schema: { $ref: "#/definitions/PrimerError" } }
        422: { description: "NO_VARIATION_DATA, NO_SEQUENCE, AMBIGUOUS_ASSEMBLY", schema: { $ref: "#/definitions/PrimerError" } }
        503: { description: "VARIATION_SOURCE_UNAVAILABLE, FEATURE_DISABLED, MONGO_UNAVAILABLE; details.retry_after_s", schema: { $ref: "#/definitions/PrimerError" } }

  /primers/variants/{variant_id}:
    x-swagger-router-controller: primers
    get:
      summary: "Resolve one Ensembl variation id to normalized variant entries (one per designable alt)"
      operationId: getPrimerVariant
      tags: [ "Primer design" ]
      consumes: [ application/json ]
      produces: [ application/json ]
      parameters:
        - { in: path, name: variant_id, required: true, type: string, pattern: "^[A-Za-z0-9][A-Za-z0-9_.:,*-]{0,254}$" }
        - { in: query, name: system_name, required: true, type: string, pattern: "^[a-z0-9_]+$", maxLength: 128 }
      responses:
        200: { description: "the variant", schema: { $ref: "#/definitions/PrimerVariantLookup" } }
        400: { description: "validator errors", schema: { $ref: "#/definitions/PrimerError" } }
        404: { description: "UNKNOWN_VARIANT, UNKNOWN_GENOME", schema: { $ref: "#/definitions/PrimerError" } }
        422: { description: "NO_VARIATION_DATA, AMBIGUOUS_VARIANT_MAPPING, VARIANT_NOT_ON_ASSEMBLY", schema: { $ref: "#/definitions/PrimerError" } }
        503: { description: "VARIATION_SOURCE_UNAVAILABLE, FEATURE_DISABLED", schema: { $ref: "#/definitions/PrimerError" } }

  /primers/genotyping/design:
    x-swagger-router-controller: primers
    post:
      summary: "Design KASP or allele-specific PCR primer sets for one variant"
      description: "Synchronous, 45 s deadline. Returns explicit sets (REF-specific, ALT-specific and common primers), order-sheet rows and a ready /primers/check request."
      operationId: designGenotypingPrimers
      tags: [ "Primer design" ]
      consumes: [ application/json ]
      produces: [ application/json ]
      parameters:
        - in: body
          name: body
          required: true
          schema: { $ref: "#/definitions/PrimerGenotypingRequest" }
      responses:
        200: { description: "variant, template, per-orientation attempts, sets (possibly none, warning NO_SETS), check request, warnings", schema: { $ref: "#/definitions/PrimerGenotypingResponse" } }
        400: { description: "validator errors, INVALID_REQUEST, INVALID_PARAMS, INVALID_VARIANT, REF_MISMATCH, ALT_REQUIRED, ALT_NOT_AT_SITE, UNSUPPORTED_ALLELE, VARIANT_TOO_CLOSE_TO_END, REGION_OUT_OF_BOUNDS, PRIMER3_INPUT_ERROR", schema: { $ref: "#/definitions/PrimerError" } }
        404: { description: "UNKNOWN_GENOME, UNKNOWN_REGION, UNKNOWN_VARIANT", schema: { $ref: "#/definitions/PrimerError" } }
        413: { description: "Request body larger than 100 kb" }
        422: { description: "NO_SEQUENCE, AMBIGUOUS_ASSEMBLY, NO_VARIATION_DATA, AMBIGUOUS_VARIANT_MAPPING, VARIANT_NOT_ON_ASSEMBLY", schema: { $ref: "#/definitions/PrimerError" } }
        500: { description: "PRIMER3_FAILED, THERMO_FAILED, INTERNAL", schema: { $ref: "#/definitions/PrimerError" } }
        503: { description: "BUSY, PRIMER3_UNAVAILABLE, THERMO_UNAVAILABLE, VARIATION_SOURCE_UNAVAILABLE, FEATURE_DISABLED; details.retry_after_s", schema: { $ref: "#/definitions/PrimerError" } }
        504: { description: "DEADLINE_EXCEEDED", schema: { $ref: "#/definitions/PrimerError" } }
```

### 2.6 Swagger definitions: variants

Every object below has a `properties` block. `x-nullable: true` marks each place a `null` can appear; the example-validation test of §7.7 is aware of `x-nullable`, because sway 1.0.0 ignores it when validating responses [R of critic evidence; treated as verified behaviour].

```yaml
definitions:
  PrimerVariantRecord:
    type: object
    properties:
      id: { type: string, description: "Ensembl id as reported; may contain ',' or '*' and may be up to 255 characters" }
      source: { type: string, description: "Ensembl overlap source code, e.g. EVA, EMS_PMID38100514_Jiao, SAP_PMID35653240_Boatwri (truncated to 100 characters)" }
      ems: { type: boolean, description: "source matches primers.variation.ems_source_pattern (^EMS_)" }
  PrimerVariantIssue:
    type: object
    properties:
      code: { type: string, enum: [REF_MISMATCH, STAR_ALLELE, ALLELE_TOO_LONG, UNSUPPORTED_ALLELE, REPEAT_TOO_LONG] }
      message: { type: string }
      details: { type: object }
  PrimerAlleleSite:
    type: object
    properties:
      position: { type: integer, description: "plus-strand coordinate of the allele-specific primer's 3' base" }
      ref_base: { type: string }
      alt_base: { type: string }
      alt_maps_to: { type: integer, x-nullable: true, description: "reference coordinate of the ALT haplotype base at this index; null for an inserted base" }
  PrimerVariantZone:
    type: object
    description: "the guarded span: VCF span, both discriminating positions and both ALT 3' anchor mappings"
    properties:
      start: { type: integer }
      end: { type: integer }
  PrimerVariantVcf:
    type: object
    properties:
      position: { type: integer }
      ref: { type: string }
      alt: { type: string }
  PrimerVariantMinimal:
    type: object
    properties:
      start: { type: integer }
      end: { type: integer }
      ref: { type: string, description: "'-' for an insertion" }
      alt: { type: string, description: "'-' for a deletion" }
  PrimerVariantMultiallelic:
    type: object
    x-nullable: true
    properties:
      alleles: { type: array, items: { type: string } }
      other_alts: { type: array, items: { type: string } }
  PrimerVariant:
    type: object
    properties:
      key: { type: string, description: "region:position:REF:ALT, left-aligned VCF with anchor base for indels; canonical identity" }
      ids: { type: array, items: { type: string }, description: "Ensembl ids with this key: requested id first, then rs ids, then the rest alphabetically" }
      synonyms: { type: array, items: { type: string }, description: "id lookup only; '.' dropped" }
      label: { type: string }
      kind: { type: string, enum: [snv, mnv, insertion, deletion, complex] }
      region: { type: string }
      vcf: { $ref: "#/definitions/PrimerVariantVcf" }
      minimal: { $ref: "#/definitions/PrimerVariantMinimal" }
      alleles: { type: array, items: { type: string }, description: "all alleles of the site as reported, reference first ('-' and '*' kept)" }
      multiallelic: { $ref: "#/definitions/PrimerVariantMultiallelic" }
      shift: { type: integer, minimum: 0 }
      zone: { $ref: "#/definitions/PrimerVariantZone" }
      discriminating:
        type: object
        properties:
          forward: { $ref: "#/definitions/PrimerAlleleSite" }
          reverse: { $ref: "#/definitions/PrimerAlleleSite" }
      records: { type: array, items: { $ref: "#/definitions/PrimerVariantRecord" } }
      ems: { type: boolean, description: "every merged record is EMS" }
      consequence: { type: string, x-nullable: true }
      ref_verified: { type: boolean }
      designable: { type: boolean }
      issues: { type: array, items: { $ref: "#/definitions/PrimerVariantIssue" } }
      requested_id: { type: string, x-nullable: true, description: "lookup and design responses only" }
      submission_sequence: { type: string, description: "design response only: 50 bp flanks with [REF/ALT] (minimal alleles, '' for '-'), non-EMS biallelic SNV neighbours as IUPAC codes" }
  PrimerVariantSource:
    type: object
    properties:
      name: { type: string, enum: [ensembl] }
      release: { type: string }
  PrimerVariantList:
    type: object
    properties:
      system_name: { type: string }
      region: { type: string }
      start: { type: integer }
      end: { type: integer }
      source: { $ref: "#/definitions/PrimerVariantSource" }
      total: { type: integer }
      returned: { type: integer }
      truncated: { type: boolean }
      variants: { type: array, items: { $ref: "#/definitions/PrimerVariant" } }
      warnings: { type: array, items: { $ref: "#/definitions/PrimerWarning" } }
  PrimerVariantLookup:
    type: object
    properties:
      requested_id: { type: string }
      system_name: { type: string }
      source: { $ref: "#/definitions/PrimerVariantSource" }
      variants: { type: array, items: { $ref: "#/definitions/PrimerVariant" } }
      warnings: { type: array, items: { $ref: "#/definitions/PrimerWarning" } }
```

`PrimerWarning` (existing, `swagger.yaml:2208`) gains `details: { type: object }`, which every new warning uses.
### 2.7 `POST /primers/genotyping/design`: request

```yaml
  PrimerVariantInput:
    description: "Exactly one of: {id [, alt]} or {region, position, ref, alt}. With id, alt selects the allele at a multi-allelic site (required there). Manual alleles: VCF style (A/C/G/T strings; position = first ref base) or Ensembl style (exactly one allele '-'; deletion: position = first deleted base; insertion: position = the base after the insertion point)."
    type: object
    additionalProperties: false
    properties:
      id: { type: string, pattern: "^[A-Za-z0-9][A-Za-z0-9_.:,*-]{0,254}$" }
      region: { type: string, maxLength: 255 }
      position: { type: integer, minimum: 1 }
      ref: { type: string, pattern: "^([ACGTacgt]{1,50}|-)$" }
      alt: { type: string, pattern: "^([ACGTacgt]{1,50}|-)$" }
  PrimerGenotypingAssay:
    description: "Unset values take the defaults of the assay type: kasp -> tails ref_fam_alt_hex, deliberate_mismatch none; as_pcr -> tails none, deliberate_mismatch auto. Both: orientation both, mismatch_position 2, num_sets 6, max_relaxation 2, neighbour_policy avoid_3p."
    type: object
    additionalProperties: false
    properties:
      type: { type: string, enum: [kasp, as_pcr] }
      orientation: { type: string, enum: [both, forward, reverse] }
      tails: { type: string, enum: [none, ref_fam_alt_hex, ref_hex_alt_fam] }
      deliberate_mismatch: { type: string, enum: [none, auto] }
      mismatch_position: { type: integer, enum: [2, 3] }
      num_sets: { type: integer, minimum: 1, maximum: 10 }
      max_relaxation: { type: integer, minimum: 0, maximum: 2 }
      neighbour_policy: { type: string, enum: [avoid_3p, ignore] }
  PrimerGenotypingParams:
    description: "Closed subset of PrimerDesignParams with the same bounds (design.js PARAM_SPECS, lines 31-54), except product_size_ranges. Unset values come from the assay preset. A param the client sets is pinned: the relaxation ladder never changes it. The effective product-range minimum is raised to 2 x max_size + 1."
    type: object
    additionalProperties: false
    properties:
      opt_size: { type: integer, minimum: 15, maximum: 36 }
      min_size: { type: integer, minimum: 15, maximum: 36 }
      max_size: { type: integer, minimum: 15, maximum: 36 }
      opt_tm: { type: number, minimum: 30, maximum: 90 }
      min_tm: { type: number, minimum: 30, maximum: 90 }
      max_tm: { type: number, minimum: 30, maximum: 90 }
      opt_gc: { type: number, minimum: 0, maximum: 100 }
      min_gc: { type: number, minimum: 0, maximum: 100 }
      max_gc: { type: number, minimum: 0, maximum: 100 }
      max_tm_diff: { type: number, minimum: 0, maximum: 30 }
      max_poly_x: { type: integer, minimum: 0, maximum: 10 }
      gc_clamp: { type: integer, minimum: 0, maximum: 5 }
      max_end_stability: { type: number, minimum: 0, maximum: 100 }
      salt_monovalent: { type: number, minimum: 0, maximum: 1000 }
      salt_divalent: { type: number, minimum: 0, maximum: 100 }
      dntp_conc: { type: number, minimum: 0, maximum: 100 }
      dna_conc: { type: number, minimum: 0, maximum: 10000 }
      product_size_ranges:
        type: array
        minItems: 1
        maxItems: 4
        items:
          type: array
          minItems: 2
          maxItems: 2
          items: { type: integer, minimum: 20, maximum: 1000 }
  PrimerGenotypingRequest:
    description: "Design allele-specific primer sets for one variant. Omit optional fields rather than sending null."
    type: object
    required: [system_name, variant]
    additionalProperties: false
    properties:
      system_name: { type: string, pattern: "^[a-z0-9_]+$", maxLength: 128 }
      variant: { $ref: "#/definitions/PrimerVariantInput" }
      assay: { $ref: "#/definitions/PrimerGenotypingAssay" }
      avoid_repeats: { type: boolean, description: "default false (section 4.5; a default to confirm, section 11 item 6); when true, the allele-specific window is always exempt from the mask" }
      repeat_mask_mode: { type: string, enum: [n_mask, three_prime] }
      template_only: { type: boolean, description: "resolve and verify the variant and return the template; Primer3 is not run" }
      label: { type: string, pattern: "^[A-Za-z0-9_.-]{1,40}$", description: "prefix of order-sheet oligo names; default: first variant id, else key" }
      params: { $ref: "#/definitions/PrimerGenotypingParams" }
```

**Handler rules**, in order. Rules without *I/O* are checked before any I/O:

| # | Rule | Error |
|---|---|---|
| 1 | Unknown top-level or nested keys (defence in depth; swagger rejects them first) | 400 `INVALID_REQUEST {field}` |
| 2 | `variant` has `id` together with `region`, `position` or `ref`; or has neither `id` nor all four manual fields | 400 `INVALID_VARIANT {reason:"id_or_manual"}` |
| 3 | Manual `ref` and `alt` both `-`, equal after uppercasing, or one is `-` while the other is not a base string | 400 `INVALID_VARIANT {reason:"alleles"}` |
| 4 | `repeat_mask_mode` without `avoid_repeats: true` | ignored, as in `/primers/design` |
| 5 | Param cross-field rules, reusing `design.checkParams` (`design.js:265-292`): `min ≤ opt ≤ max`, `gc_clamp ≤ min_size`, each range `a < b` | 400 `INVALID_PARAMS {param}` |
| 6 | `assay.tails ≠ none` with `type: as_pcr`, or `deliberate_mismatch: auto` with `type: kasp` | allowed; both are the user's choice and the response echoes them |
| 7 | *I/O* `system_name` in the catalog; the assembly resolves with a FASTA | 404 `UNKNOWN_GENOME`; 422 `NO_SEQUENCE` / `AMBIGUOUS_ASSEMBLY` |
| 8 | *I/O* `id` given but the genome has no variation data | 422 `NO_VARIATION_DATA` |
| 9 | *I/O* id resolution (§3) | 404 `UNKNOWN_VARIANT`; 422 `AMBIGUOUS_VARIANT_MAPPING` / `VARIANT_NOT_ON_ASSEMBLY`; 503 `VARIATION_SOURCE_UNAVAILABLE` |
| 10 | *I/O* the id's site has more than one designable alt and no `alt` was given | 400 `ALT_REQUIRED {id, alts}` |
| 11 | *I/O* `alt` given with `id` but not an allele of the site | 400 `ALT_NOT_AT_SITE {id, alt, alleles}` |
| 12 | *I/O* an allele is `*` or contains `N` | 400 `UNSUPPORTED_ALLELE {allele}` |
| 13 | *I/O* manual `region` not in the assembly, or `position` beyond it | 404 `UNKNOWN_REGION`; 400 `REGION_OUT_OF_BOUNDS` |
| 14 | *I/O* REF differs from the FASTA | 400 `REF_MISMATCH {region, position, given, genome}` |
| 15 | *I/O* neither orientation has room for the smallest product | 400 `VARIANT_TOO_CLOSE_TO_END {region, position, region_length, needed}` |

Steps 7–14 run **before** `semaphore.acquire` (§4.2), so a slow Ensembl never occupies a Primer3 slot.

### 2.8 `POST /primers/genotyping/design`: response definitions

```yaml
  PrimerNeighbourHit:
    type: object
    properties:
      key: { type: string }
      ids: { type: array, items: { type: string } }
      label: { type: string }
      start: { type: integer }
      end: { type: integer }
      alleles: { type: string, description: "Ensembl allele string, e.g. A/T or -/GT" }
      ems: { type: boolean }
      distance_from_3p:
        type: integer
        x-nullable: true
        description: "1 = the primer's own 3'-terminal base, measured along the primer (minimum over the neighbour's span); null in orientation blockers, where no primer exists yet"
  PrimerLikelihoodAt:
    type: object
    properties:
      mm_pos: { type: array, items: { type: integer }, description: "mismatch distances from the 3' end, 5'->3' (1 = terminal)" }
      likelihood: { type: string, enum: [likely, likely_weak, unlikely] }
  PrimerDeliberateMismatch:
    type: object
    x-nullable: true
    properties:
      position: { type: integer, enum: [2, 3], description: "distance from the 3' end" }
      original_base: { type: string }
      new_base: { type: string }
      template_base: { type: string, description: "template base opposite the mismatch position (complement of original_base)" }
      terminal_pair: { type: string, pattern: "^[ACGT]{2}$", description: "unordered pair used as the table row: the primer's 3' base and the other allele's base on the template strand the primer anneals to, written in alphabetical order, e.g. AG or CT" }
      terminal_mismatch_class: { type: string, enum: [max, strong, medium, weak] }
      added_mismatch_class: { type: string, enum: [max, strong, medium, weak] }
      source: { type: string, description: "Little 1995 Table 9.8.1, or the same table applied at -3 (extrapolated)" }
  PrimerGenotypingDiscrimination:
    type: object
    x-nullable: true
    description: "allele-specific primers only"
    properties:
      own_allele: { $ref: "#/definitions/PrimerLikelihoodAt" }
      other_allele: { $ref: "#/definitions/PrimerLikelihoodAt" }
      terminal_mismatch_class: { type: string, enum: [max, strong, medium, weak], description: "Little 1995 class of the 3' base against the other allele" }
      in_shift_tract: { type: boolean, description: "the primer's 3' base, on its own haplotype, lies inside the variant's shift tract (section 3.5: the reference bases an indel with shift > 0 can slide through, anchor base excluded); for as_alt an inserted base counts as inside. On the other allele's genomes the check reports this primer as uncertain (sections 4.12, 5.7)" }
  PrimerTailedStructures:
    type: object
    x-nullable: true
    properties:
      hairpin_th: { type: number }
      self_any_th: { type: number }
      self_end_th: { type: number }
  PrimerOligoTemplateSpan:
    type: object
    properties:
      start: { type: integer }
      end: { type: integer }
      sequence: { type: string, enum: [ref, alt], description: "which template the span indexes: template.seq or template.alt_seq" }
  PrimerGenotypingOligo:
    type: object
    properties:
      role: { type: string, enum: [as_ref, as_alt, common] }
      allele: { type: string, x-nullable: true, description: "plus-strand allele this primer detects (vcf.ref or vcf.alt); null for the common primer" }
      three_prime_base: { type: string }
      haplotype: { type: string, enum: [ref, alt, both] }
      target_seq: { type: string, description: "anneals; as ordered (includes any deliberate mismatch); the only sequence sent to /primers/check" }
      matched_seq: { type: string, description: "perfect-match sequence before any deliberate mismatch" }
      dye: { type: string, enum: [FAM, HEX], x-nullable: true }
      tail_seq: { type: string, x-nullable: true }
      order_seq: { type: string, description: "tail_seq + target_seq" }
      len: { type: integer, description: "length of target_seq" }
      order_len: { type: integer }
      tm: { type: number, description: "Primer3 Tm of target_seq, or the own-allele duplex Tm for a deliberate-mismatch primer (see tm_method)" }
      tm_method: { type: string, enum: [primer3, ntthal_duplex] }
      matched_tm: { type: number, x-nullable: true, description: "deliberate-mismatch primers only: Primer3 Tm of matched_seq; the Tm floor is applied to this value" }
      gc: { type: number, description: "100 x GC / len of target_seq, rounded half away from zero to 2 decimals from the exact ratio (section 2.1); computed from the sequence, not parsed from Primer3" }
      hairpin_th: { type: number }
      self_any_th: { type: number }
      self_end_th: { type: number }
      end_stability: { type: number }
      primer3_problems: { type: string, x-nullable: true, description: "PRIMER_*_PROBLEMS text from the scoring run, e.g. ' Temperature too low;'; informational, the hard floors of section 4.6 decide acceptance" }
      template: { $ref: "#/definitions/PrimerOligoTemplateSpan" }
      genomic: { $ref: "#/definitions/PrimerGenomicLocation" }
      inserted_bases: { type: integer, description: "bases of this primer that have no reference coordinate (inside an insertion)" }
      deliberate_mismatch: { $ref: "#/definitions/PrimerDeliberateMismatch" }
      discrimination: { $ref: "#/definitions/PrimerGenotypingDiscrimination" }
      tailed: { $ref: "#/definitions/PrimerTailedStructures" }
      neighbours: { type: array, items: { $ref: "#/definitions/PrimerNeighbourHit" }, description: "known variants under this primer, sorted by distance_from_3p ascending" }
  PrimerGenotypingOrderRow:
    type: object
    properties:
      name: { type: string }
      set_id: { type: string }
      set_key: { type: string }
      role: { type: string, enum: [as_ref, as_alt, common] }
      allele: { type: string, x-nullable: true }
      dye: { type: string, enum: [FAM, HEX], x-nullable: true }
      order_seq: { type: string }
      target_seq: { type: string }
      tail_seq: { type: string, x-nullable: true }
      length: { type: integer }
      tm: { type: number }
      gc: { type: number }
      orientation: { type: string, enum: [forward, reverse] }
      product_size_ref: { type: integer }
      product_size_alt: { type: integer }
      variant_key: { type: string }
      notes: { type: string }
  PrimerCheckPairFragment:
    type: object
    description: "one pairs[] item of the check request this set needs; mirrors PrimerCheckRequest.pairs.items plus nothing"
    properties:
      id: { type: string }
      left: { type: string }
      right: { type: string }
      expected:
        type: object
        properties:
          region: { type: string }
          start: { type: integer }
          end: { type: integer }
  PrimerGenotypingSetRef:
    type: object
    properties:
      id: { type: string }
      ref_pair: { type: string }
      alt_pair: { type: string }
  PrimerGenotypingCheckFragment:
    type: object
    properties:
      set: { $ref: "#/definitions/PrimerGenotypingSetRef" }
      pairs: { type: array, items: { $ref: "#/definitions/PrimerCheckPairFragment" } }
  PrimerGenotypingProduct:
    type: object
    properties:
      size: { type: integer }
      template: { $ref: "#/definitions/PrimerOligoTemplateSpan" }
      genomic: { $ref: "#/definitions/PrimerGenomicLocation" }
      inserted_bases: { type: integer, description: "ALT products only: bases of the product with no reference coordinate" }
  PrimerGenotypingPairThermo:
    type: object
    properties:
      compl_any_th: { type: number }
      compl_end_th: { type: number }
  PrimerGenotypingTailedThermo:
    type: object
    x-nullable: true
    properties:
      ref_alt_any_th: { type: number }
      ref_alt_end_th: { type: number }
      ref_common_any_th: { type: number }
      ref_common_end_th: { type: number }
      alt_common_any_th: { type: number }
      alt_common_end_th: { type: number }
  PrimerGenotypingSetIssue:
    type: object
    properties:
      code: { type: string }
      severity: { type: string, enum: [info, warn, high] }
      message: { type: string }
      details: { type: object }
  PrimerGenotypingSet:
    type: object
    properties:
      id: { type: string, pattern: "^S[0-9]+$" }
      key: { type: string, pattern: "^[0-9a-f]{12}$" }
      rank: { type: integer, description: "0-based; id is S(rank+1)" }
      orientation: { type: string, enum: [forward, reverse] }
      relaxation_level: { type: integer, minimum: 0, maximum: 2 }
      quality: { type: string, enum: [good, usable, poor] }
      score: { type: number, description: "ranking score, lower is better: the exact sum of rules S1-S6 of section 4.16, rounded half away from zero to 2 decimals (section 2.1)" }
      primers:
        type: object
        properties:
          as_ref: { $ref: "#/definitions/PrimerGenotypingOligo" }
          as_alt: { $ref: "#/definitions/PrimerGenotypingOligo" }
          common: { $ref: "#/definitions/PrimerGenotypingOligo" }
      products:
        type: object
        properties:
          ref: { $ref: "#/definitions/PrimerGenotypingProduct" }
          alt: { $ref: "#/definitions/PrimerGenotypingProduct" }
      thermo:
        type: object
        properties:
          ref_common: { $ref: "#/definitions/PrimerGenotypingPairThermo" }
          alt_common: { $ref: "#/definitions/PrimerGenotypingPairThermo" }
          tailed: { $ref: "#/definitions/PrimerGenotypingTailedThermo" }
      tm_balance:
        type: object
        properties:
          as_tm_diff: { type: number, description: "|as_ref.tm - as_alt.tm| from unrounded values, rounded to 2 decimals (sections 2.1, 4.13)" }
          common_minus_as: { type: number, description: "common.tm - max(nominal Tm of the two allele-specific primers) from unrounded values, rounded to 2 decimals (sections 2.1, 4.13)" }
      neighbour_sites: { type: integer, description: "distinct non-EMS neighbour keys under the three primers" }
      primer3_penalty: { type: number, description: "PRIMER_PAIR_i_PENALTY of the design-run pair i this set was built from (its as_ref and common primers), rounded half away from zero to 4 decimals; that is PRIMER_PAIR_0_PENALTY only when i = 0 (section 4.16 S1)" }
      warnings: { type: array, items: { $ref: "#/definitions/PrimerWarning" }, description: "the issues of severity warn or high, as {code, message, details}; info issues appear only in issues" }
      issues: { type: array, items: { $ref: "#/definitions/PrimerGenotypingSetIssue" }, description: "every issue of the set, with its severity" }
      check: { $ref: "#/definitions/PrimerGenotypingCheckFragment" }
      order: { type: array, items: { $ref: "#/definitions/PrimerGenotypingOrderRow" } }
  PrimerGenotypingExplainSide:
    type: object
    description: "parsed Primer3 explain line; 'raw' plus one integer per reason, as PrimerDesignResponse.explain"
    properties:
      raw: { type: string }
    additionalProperties: { type: integer }
  PrimerGenotypingRejected:
    type: object
    properties:
      force: { type: integer, description: "Primer3 ignored the forced 3' end (never observed)" }
      overlap: { type: integer, description: "common primer footprint touched the zone or an AS footprint on either haplotype" }
      common_neighbour_3p: { type: integer }
      alt_scoring_failed: { type: integer, description: "the ALT or mismatch check_primers run returned no value for the derived oligo, or PRIMER_ERROR (section 4.7)" }
      below_floor: { type: integer, description: "an AS primer was under the Tm or GC floor" }
      duplicate: { type: integer }
  PrimerGenotypingAttempt:
    type: object
    properties:
      level: { type: integer }
      changes: { type: object, description: "param changes of this level (empty at level 0); pinned params excluded" }
      explain:
        type: object
        properties:
          left: { $ref: "#/definitions/PrimerGenotypingExplainSide" }
          right: { $ref: "#/definitions/PrimerGenotypingExplainSide" }
          pair: { $ref: "#/definitions/PrimerGenotypingExplainSide" }
      pairs_returned: { type: integer, description: "pairs Primer3 returned at this level (PRIMER_NUM_RETURN 20); always the sum of rejected, not_scored and sets" }
      rejected: { $ref: "#/definitions/PrimerGenotypingRejected" }
      not_scored: { type: integer, description: "pairs left unscored because the orientation already held max_scored_per_orientation candidates or the design budget could not cover them (section 4.8)" }
      sets: { type: integer, description: "candidates from this level that were scored and kept" }
  PrimerGenotypingOrientation:
    type: object
    properties:
      status: { type: string, enum: [ok, no_sets, blocked, skipped] }
      reason:
        type: string
        x-nullable: true
        enum: [neighbour_at_3p, too_close_to_end, n_in_primer_window, not_requested, below_floor, budget_exhausted]
      discriminating_position: { type: integer }
      relaxation_level: { type: integer, x-nullable: true }
      sets_found: { type: integer, description: "scored candidates in this orientation over all attempted levels, at most max_scored_per_orientation (8); the response's sets are the best num_sets over both orientations" }
      blockers: { type: array, items: { $ref: "#/definitions/PrimerNeighbourHit" } }
      attempts: { type: array, items: { $ref: "#/definitions/PrimerGenotypingAttempt" } }
  PrimerGenotypingTemplateFeatures:
    type: object
    properties:
      variant: { type: object, properties: { start: { type: integer }, end: { type: integer } }, description: "template span of vcf.ref" }
      zone: { type: object, properties: { start: { type: integer }, end: { type: integer } }, description: "guarded span in template coordinates" }
      discriminating: { type: object, properties: { forward: { type: integer }, reverse: { type: integer } } }
      alt_offset: { type: integer, description: "len(vcf.alt) - len(vcf.ref)" }
      exempt: { type: array, items: { $ref: "#/definitions/PrimerInterval" }, description: "template spans exempted from the repeat mask" }
  PrimerGenotypingTemplate:
    type: object
    properties:
      system_name: { type: string }
      region: { type: string }
      start: { type: integer }
      end: { type: integer }
      strand: { type: integer, enum: [1] }
      length: { type: integer }
      alt_length: { type: integer }
      seq: { type: string }
      alt_seq: { type: string }
      masked: { type: boolean }
      mask_source: { type: string, enum: [softmask, blast_depth], x-nullable: true }
      mask: { type: array, items: { $ref: "#/definitions/PrimerInterval" } }
      masked_fraction: { type: number }
      features: { $ref: "#/definitions/PrimerGenotypingTemplateFeatures" }
  PrimerKaspMix:
    type: object
    x-nullable: true
    properties:
      stock_uM: { type: number }
      as_ref_uL: { type: number }
      as_alt_uL: { type: number }
      common_uL: { type: number }
      water_uL: { type: number }
      total_uL: { type: number }
      source: { type: string }
  PrimerGenotypingAssayEffective:
    type: object
    properties:
      type: { type: string, enum: [kasp, as_pcr] }
      orientation: { type: string, enum: [both, forward, reverse] }
      tails: { type: string, enum: [none, ref_fam_alt_hex, ref_hex_alt_fam] }
      deliberate_mismatch: { type: string, enum: [none, auto] }
      mismatch_position: { type: integer, enum: [2, 3] }
      num_sets: { type: integer }
      max_relaxation: { type: integer }
      neighbour_policy: { type: string, enum: [avoid_3p, ignore] }
      ems_target: { type: boolean, description: "every Ensembl record of the target is EMS; natural neighbours then warn instead of blocking" }
      kasp_mix: { $ref: "#/definitions/PrimerKaspMix" }
  PrimerGenotypingNeighbourSummary:
    type: object
    properties:
      data: { type: string, enum: [ensembl, none, unavailable] }
      window: { type: object, properties: { start: { type: integer }, end: { type: integer } } }
      variants: { type: integer }
      non_ems: { type: integer }
      ems: { type: integer }
      dense_non_ems: { type: integer, description: "non-EMS neighbour keys whose minimal span overlaps [vcf.position - 30, vcf.position + len(vcf.ref) - 1 + 30], the target excluded" }
  PrimerGenotypingLadderStep:
    type: object
    properties:
      level: { type: integer }
      changes: { type: object }
  PrimerGenotypingSettings:
    type: object
    properties:
      preset: { type: string, enum: [kasp, as_pcr] }
      params: { type: object, description: "level-0 effective params, product-range minimum already raised" }
      pinned: { type: array, items: { type: string } }
      ladder: { type: array, items: { $ref: "#/definitions/PrimerGenotypingLadderStep" } }
      floors: { type: object, properties: { as_min_tm: { type: number }, as_min_gc: { type: number } } }
  PrimerGenotypingEngine:
    type: object
    properties:
      primer3: { type: string, x-nullable: true }
      thermo: { type: string, x-nullable: true, description: "ntthal version string" }
      genotyping_design: { type: string, description: "design algorithm version, currently '1'" }
      variation_source: { type: string, x-nullable: true, description: "ensembl 115" }
  PrimerGenotypingCheckBlock:
    type: object
    x-nullable: true
    properties:
      request: { $ref: "#/definitions/PrimerCheckRequest" }
      set_ids: { type: array, items: { type: string } }
      unique_primers: { type: integer }
      omitted_set_ids: { type: array, items: { type: string }, description: "sets left out because of the 5-set, 10-pair or 13-primer caps" }
  PrimerGenotypingResponse:
    type: object
    properties:
      variant: { $ref: "#/definitions/PrimerVariant" }
      template: { $ref: "#/definitions/PrimerGenotypingTemplate" }
      assay: { $ref: "#/definitions/PrimerGenotypingAssayEffective" }
      neighbours: { $ref: "#/definitions/PrimerGenotypingNeighbourSummary" }
      orientations:
        type: object
        x-nullable: true
        description: "null for template_only requests"
        properties:
          forward: { $ref: "#/definitions/PrimerGenotypingOrientation" }
          reverse: { $ref: "#/definitions/PrimerGenotypingOrientation" }
      sets: { type: array, items: { $ref: "#/definitions/PrimerGenotypingSet" } }
      check: { $ref: "#/definitions/PrimerGenotypingCheckBlock" }
      settings: { $ref: "#/definitions/PrimerGenotypingSettings" }
      engine: { $ref: "#/definitions/PrimerGenotypingEngine" }
      warnings: { type: array, items: { $ref: "#/definitions/PrimerWarning" } }
```

Two definitions are reused unchanged: `PrimerGenomicLocation` (`swagger.yaml:2218`) and `PrimerInterval` (`:2076`). Nothing new is defined but unreferenced, so sway reports no `UNUSED_DEFINITION` for the `Primer*` paths, which `contract.test.js:84-90` forbids.
### 2.9 Complete design example: rs871475760, KASP, two sets

Everything below is real [V]: sequences and coordinates come from the genome, Tm/GC/structure/penalty values from `primer3_core` 2.6.1 and `ntthal`, the `explain` strings from the exact records this spec prescribes (forced 3′ end plus the `SEQUENCE_TARGET` guard), the neighbours from live Ensembl, and the set keys from sha256. Only the warning message wording is editorial. `template.seq`/`alt_seq` are 801 nt; they are elided here as `"GCGAGTTCTCAAG…CATATGAT"` and appear in full in the recorded fixtures (§7.4, §6.3).

**Request**

```json
{
  "system_name": "sorghum_bicolor",
  "variant": { "id": "rs871475760", "alt": "A" },
  "assay": { "type": "kasp", "num_sets": 2 }
}
```

**Response `200`**

```json
{
  "variant": {
    "key": "1:11109:C:A",
    "requested_id": "rs871475760",
    "ids": ["rs871475760"],
    "synonyms": ["tmp_1_11109_C_A"],
    "label": "1:11109 C/A",
    "kind": "snv",
    "region": "1",
    "vcf": { "position": 11109, "ref": "C", "alt": "A" },
    "minimal": { "start": 11109, "end": 11109, "ref": "C", "alt": "A" },
    "alleles": ["C", "A"],
    "multiallelic": null,
    "shift": 0,
    "zone": { "start": 11109, "end": 11109 },
    "discriminating": {
      "forward": { "position": 11109, "ref_base": "C", "alt_base": "A", "alt_maps_to": 11109 },
      "reverse": { "position": 11109, "ref_base": "C", "alt_base": "A", "alt_maps_to": 11109 }
    },
    "records": [{ "id": "rs871475760", "source": "EVA", "ems": false }],
    "ems": false,
    "consequence": "downstream_gene_variant",
    "ref_verified": true,
    "designable": true,
    "issues": [],
    "submission_sequence": "YCCTCAAAAAGCTTCTCTAAGTGGTTATCCGAATATAGTCATACTCTATT[C/A]TGAATTTCTCGCTAGTCAAAGATAACAAAAATAGCATATTCTGGATTTCT"
  },
  "template": {
    "system_name": "sorghum_bicolor",
    "region": "1", "start": 10709, "end": 11509, "strand": 1,
    "length": 801, "alt_length": 801,
    "seq": "GCGAGTTCTCAAG…CATATGAT",
    "alt_seq": "GCGAGTTCTCAAG…CATATGAT",
    "masked": false, "mask_source": null, "mask": [], "masked_fraction": 0,
    "features": {
      "variant": { "start": 401, "end": 401 },
      "zone": { "start": 401, "end": 401 },
      "discriminating": { "forward": 401, "reverse": 401 },
      "alt_offset": 0,
      "exempt": [[365, 37], [401, 37]]
    }
  },
  "assay": {
    "type": "kasp", "orientation": "both", "tails": "ref_fam_alt_hex",
    "deliberate_mismatch": "none", "mismatch_position": 2,
    "num_sets": 2, "max_relaxation": 2, "neighbour_policy": "avoid_3p", "ems_target": false,
    "kasp_mix": { "stock_uM": 100, "as_ref_uL": 12, "as_alt_uL": 12, "common_uL": 30, "water_uL": 46, "total_uL": 100,
                  "source": "Makhoul et al. 2020 (12:12:30 at 100 uM); water to 100 uL inferred" }
  },
  "neighbours": { "data": "ensembl", "window": { "start": 10709, "end": 11509 }, "variants": 57, "non_ems": 45, "ems": 12, "dense_non_ems": 0 },
  "orientations": {
    "forward": {
      "status": "ok", "reason": null, "discriminating_position": 11109, "relaxation_level": 0, "sets_found": 8, "blockers": [],
      "attempts": [
        { "level": 0, "changes": {},
          "explain": {
            "left":  { "raw": "considered 13, GC content failed 5, low tm 6, ok 2", "considered": 13, "GC content failed": 5, "low tm": 6, "ok": 2 },
            "right": { "raw": "considered 4771, GC content failed 1718, low tm 1123, high tm 777, ok 1153", "considered": 4771, "GC content failed": 1718, "low tm": 1123, "high tm": 777, "ok": 1153 },
            "pair":  { "raw": "considered 1463, unacceptable product size 1443, ok 20", "considered": 1463, "unacceptable product size": 1443, "ok": 20 } },
          "pairs_returned": 20,
          "rejected": { "force": 0, "overlap": 0, "common_neighbour_3p": 0, "alt_scoring_failed": 0, "below_floor": 0, "duplicate": 0 },
          "not_scored": 12,
          "sets": 8 }
      ]
    },
    "reverse": {
      "status": "ok", "reason": null, "discriminating_position": 11109, "relaxation_level": 0, "sets_found": 8, "blockers": [],
      "attempts": [
        { "level": 0, "changes": {},
          "explain": {
            "left":  { "raw": "considered 4771, GC content failed 598, low tm 1550, high tm 1031, ok 1592", "considered": 4771, "GC content failed": 598, "low tm": 1550, "high tm": 1031, "ok": 1592 },
            "right": { "raw": "considered 13, low tm 6, ok 7", "considered": 13, "low tm": 6, "ok": 7 },
            "pair":  { "raw": "considered 905, unacceptable product size 883, tm diff too large 1, ok 21", "considered": 905, "unacceptable product size": 883, "tm diff too large": 1, "ok": 21 } },
          "pairs_returned": 20,
          "rejected": { "force": 0, "overlap": 0, "common_neighbour_3p": 6, "alt_scoring_failed": 0, "below_floor": 0, "duplicate": 0 },
          "not_scored": 6,
          "sets": 8 }
      ]
    }
  },
  "sets": [
    {
      "id": "S1", "key": "f9df650ad116", "rank": 0, "orientation": "reverse", "relaxation_level": 0,
      "quality": "usable", "score": 9.18,
      "primers": {
        "as_ref": {
          "role": "as_ref", "allele": "C", "three_prime_base": "G", "haplotype": "ref",
          "target_seq": "ATCTTTGACTAGCGAGAAATTCAG", "matched_seq": "ATCTTTGACTAGCGAGAAATTCAG",
          "dye": "FAM", "tail_seq": "GAAGGTGACCAAGTTCATGCT", "order_seq": "GAAGGTGACCAAGTTCATGCTATCTTTGACTAGCGAGAAATTCAG",
          "len": 24, "order_len": 45,
          "tm": 57.1, "tm_method": "primer3", "matched_tm": null, "gc": 37.5,
          "hairpin_th": 0, "self_any_th": 0, "self_end_th": 0, "end_stability": 3.02, "primer3_problems": null,
          "template": { "start": 401, "end": 424, "sequence": "ref" },
          "genomic": { "region": "1", "start": 11109, "end": 11132, "strand": -1, "blocks": [{ "start": 11109, "end": 11132 }] },
          "inserted_bases": 0,
          "deliberate_mismatch": null,
          "discrimination": { "own_allele": { "mm_pos": [], "likelihood": "likely" },
                              "other_allele": { "mm_pos": [1], "likelihood": "likely_weak" },
                              "terminal_mismatch_class": "max", "in_shift_tract": false },
          "tailed": { "hairpin_th": 42.62, "self_any_th": 3.63, "self_end_th": 11.64 },
          "neighbours": []
        },
        "as_alt": {
          "role": "as_alt", "allele": "A", "three_prime_base": "T", "haplotype": "alt",
          "target_seq": "ATCTTTGACTAGCGAGAAATTCAT", "matched_seq": "ATCTTTGACTAGCGAGAAATTCAT",
          "dye": "HEX", "tail_seq": "GAAGGTCGGAGTCAACGGATT", "order_seq": "GAAGGTCGGAGTCAACGGATTATCTTTGACTAGCGAGAAATTCAT",
          "len": 24, "order_len": 45,
          "tm": 56.51, "tm_method": "primer3", "matched_tm": null, "gc": 33.33,
          "hairpin_th": 0, "self_any_th": 0, "self_end_th": 0, "end_stability": 2.57, "primer3_problems": " Temperature too low;",
          "template": { "start": 401, "end": 424, "sequence": "alt" },
          "genomic": { "region": "1", "start": 11109, "end": 11132, "strand": -1, "blocks": [{ "start": 11109, "end": 11132 }] },
          "inserted_bases": 0,
          "deliberate_mismatch": null,
          "discrimination": { "own_allele": { "mm_pos": [], "likelihood": "likely" },
                              "other_allele": { "mm_pos": [1], "likelihood": "likely_weak" },
                              "terminal_mismatch_class": "max", "in_shift_tract": false },
          "tailed": { "hairpin_th": 51.5, "self_any_th": 26.85, "self_end_th": 11.47 },
          "neighbours": []
        },
        "common": {
          "role": "common", "allele": null, "three_prime_base": "A", "haplotype": "both",
          "target_seq": "AGCTTCTCTAAGTGGTTATCCGA", "matched_seq": "AGCTTCTCTAAGTGGTTATCCGA",
          "dye": null, "tail_seq": null, "order_seq": "AGCTTCTCTAAGTGGTTATCCGA",
          "len": 23, "order_len": 23,
          "tm": 58.72, "tm_method": "primer3", "matched_tm": null, "gc": 43.48,
          "hairpin_th": 34.99, "self_any_th": 0, "self_end_th": 0, "end_stability": 4.55, "primer3_problems": null,
          "template": { "start": 360, "end": 382, "sequence": "ref" },
          "genomic": { "region": "1", "start": 11068, "end": 11090, "strand": 1, "blocks": [{ "start": 11068, "end": 11090 }] },
          "inserted_bases": 0,
          "deliberate_mismatch": null,
          "discrimination": null,
          "tailed": null,
          "neighbours": [
            { "key": "1:11069:G:A", "ids": ["tmp_1_11069_G_A"], "label": "1:11069 G/A", "start": 11069, "end": 11069, "alleles": "G/A", "ems": true, "distance_from_3p": 22 }
          ]
        }
      },
      "products": {
        "ref": { "size": 65, "template": { "start": 360, "end": 424, "sequence": "ref" },
                 "genomic": { "region": "1", "start": 11068, "end": 11132, "strand": 1, "blocks": [{ "start": 11068, "end": 11132 }] }, "inserted_bases": 0 },
        "alt": { "size": 65, "template": { "start": 360, "end": 424, "sequence": "alt" },
                 "genomic": { "region": "1", "start": 11068, "end": 11132, "strand": 1, "blocks": [{ "start": 11068, "end": 11132 }] }, "inserted_bases": 0 }
      },
      "thermo": {
        "ref_common": { "compl_any_th": 0, "compl_end_th": 0 },
        "alt_common": { "compl_any_th": 0, "compl_end_th": 0 },
        "tailed": { "ref_alt_any_th": 6.89, "ref_alt_end_th": 8.21, "ref_common_any_th": 0, "ref_common_end_th": 0, "alt_common_any_th": 5.46, "alt_common_end_th": 0 }
      },
      "tm_balance": { "as_tm_diff": 0.59, "common_minus_as": 1.62 },
      "neighbour_sites": 0,
      "primer3_penalty": 7.1797,
      "warnings": [
        { "code": "ALT_PRIMER_SUBOPTIMAL", "message": "the derived as_alt primer would not have been chosen de novo (Primer3: temperature too low); it clears the hard floors",
          "details": { "oligo": "as_alt", "problems": " Temperature too low;" } },
        { "code": "TAILED_STRUCTURE", "message": "HEX-tailed as_alt: hairpin Tm 51.50 °C is above 47 °C",
          "details": { "oligo": "as_alt", "metric": "hairpin_th", "value": 51.5, "severity": "warn" } }
      ],
      "issues": [
        { "code": "ALT_PRIMER_SUBOPTIMAL", "severity": "warn", "message": "as_alt: Primer3 reports temperature too low", "details": { "oligo": "as_alt", "problems": " Temperature too low;" } },
        { "code": "TAILED_STRUCTURE", "severity": "warn", "message": "HEX-tailed as_alt hairpin 51.50 °C", "details": { "oligo": "as_alt", "metric": "hairpin_th", "value": 51.5 } }
      ],
      "check": {
        "set": { "id": "S1", "ref_pair": "S1_REF", "alt_pair": "S1_ALT" },
        "pairs": [
          { "id": "S1_REF", "left": "AGCTTCTCTAAGTGGTTATCCGA", "right": "ATCTTTGACTAGCGAGAAATTCAG", "expected": { "region": "1", "start": 11068, "end": 11132 } },
          { "id": "S1_ALT", "left": "AGCTTCTCTAAGTGGTTATCCGA", "right": "ATCTTTGACTAGCGAGAAATTCAT", "expected": { "region": "1", "start": 11068, "end": 11132 } }
        ]
      },
      "order": [
        { "name": "rs871475760_S1_REF_FAM", "set_id": "S1", "set_key": "f9df650ad116", "role": "as_ref", "allele": "C", "dye": "FAM",
          "order_seq": "GAAGGTGACCAAGTTCATGCTATCTTTGACTAGCGAGAAATTCAG", "target_seq": "ATCTTTGACTAGCGAGAAATTCAG", "tail_seq": "GAAGGTGACCAAGTTCATGCT",
          "length": 45, "tm": 57.1, "gc": 37.5, "orientation": "reverse", "product_size_ref": 65, "product_size_alt": 65, "variant_key": "1:11109:C:A", "notes": "" },
        { "name": "rs871475760_S1_ALT_HEX", "set_id": "S1", "set_key": "f9df650ad116", "role": "as_alt", "allele": "A", "dye": "HEX",
          "order_seq": "GAAGGTCGGAGTCAACGGATTATCTTTGACTAGCGAGAAATTCAT", "target_seq": "ATCTTTGACTAGCGAGAAATTCAT", "tail_seq": "GAAGGTCGGAGTCAACGGATT",
          "length": 45, "tm": 56.51, "gc": 33.33, "orientation": "reverse", "product_size_ref": 65, "product_size_alt": 65, "variant_key": "1:11109:C:A", "notes": "tail hairpin 51.5 °C" },
        { "name": "rs871475760_S1_COM", "set_id": "S1", "set_key": "f9df650ad116", "role": "common", "allele": null, "dye": null,
          "order_seq": "AGCTTCTCTAAGTGGTTATCCGA", "target_seq": "AGCTTCTCTAAGTGGTTATCCGA", "tail_seq": null,
          "length": 23, "tm": 58.72, "gc": 43.48, "orientation": "reverse", "product_size_ref": 65, "product_size_alt": 65, "variant_key": "1:11109:C:A", "notes": "" }
      ]
    },
    {
      "id": "S2", "key": "1accc54c262d", "rank": 1, "orientation": "forward", "relaxation_level": 0,
      "quality": "poor", "score": 19.63,
      "primers": {
        "as_ref": {
          "role": "as_ref", "allele": "C", "three_prime_base": "C", "haplotype": "ref",
          "target_seq": "GGTTATCCGAATATAGTCATACTCTATTC", "matched_seq": "GGTTATCCGAATATAGTCATACTCTATTC",
          "dye": "FAM", "tail_seq": "GAAGGTGACCAAGTTCATGCT", "order_seq": "GAAGGTGACCAAGTTCATGCTGGTTATCCGAATATAGTCATACTCTATTC",
          "len": 29, "order_len": 50,
          "tm": 57.28, "tm_method": "primer3", "matched_tm": null, "gc": 34.48,
          "hairpin_th": 39.2, "self_any_th": 12.01, "self_end_th": 12.01, "end_stability": 1.75, "primer3_problems": null,
          "template": { "start": 373, "end": 401, "sequence": "ref" },
          "genomic": { "region": "1", "start": 11081, "end": 11109, "strand": 1, "blocks": [{ "start": 11081, "end": 11109 }] },
          "inserted_bases": 0,
          "deliberate_mismatch": null,
          "discrimination": { "own_allele": { "mm_pos": [], "likelihood": "likely" },
                              "other_allele": { "mm_pos": [1], "likelihood": "likely_weak" },
                              "terminal_mismatch_class": "max", "in_shift_tract": false },
          "tailed": { "hairpin_th": 36.16, "self_any_th": 3.63, "self_end_th": 12.01 },
          "neighbours": []
        },
        "as_alt": {
          "role": "as_alt", "allele": "A", "three_prime_base": "A", "haplotype": "alt",
          "target_seq": "GGTTATCCGAATATAGTCATACTCTATTA", "matched_seq": "GGTTATCCGAATATAGTCATACTCTATTA",
          "dye": "HEX", "tail_seq": "GAAGGTCGGAGTCAACGGATT", "order_seq": "GAAGGTCGGAGTCAACGGATTGGTTATCCGAATATAGTCATACTCTATTA",
          "len": 29, "order_len": 50,
          "tm": 56.34, "tm_method": "primer3", "matched_tm": null, "gc": 31.03,
          "hairpin_th": 30.63, "self_any_th": 0.58, "self_end_th": 0, "end_stability": 0.98, "primer3_problems": " Temperature too low;",
          "template": { "start": 373, "end": 401, "sequence": "alt" },
          "genomic": { "region": "1", "start": 11081, "end": 11109, "strand": 1, "blocks": [{ "start": 11081, "end": 11109 }] },
          "inserted_bases": 0,
          "deliberate_mismatch": null,
          "discrimination": { "own_allele": { "mm_pos": [], "likelihood": "likely" },
                              "other_allele": { "mm_pos": [1], "likelihood": "likely_weak" },
                              "terminal_mismatch_class": "max", "in_shift_tract": false },
          "tailed": { "hairpin_th": 66.22, "self_any_th": 27.34, "self_end_th": 10.26 },
          "neighbours": []
        },
        "common": {
          "role": "common", "allele": null, "three_prime_base": "A", "haplotype": "both",
          "target_seq": "TCTTTGTCTACTGAGAAATCCAGA", "matched_seq": "TCTTTGTCTACTGAGAAATCCAGA",
          "dye": null, "tail_seq": null, "order_seq": "TCTTTGTCTACTGAGAAATCCAGA",
          "len": 24, "order_len": 24,
          "tm": 57.08, "tm_method": "primer3", "matched_tm": null, "gc": 37.5,
          "hairpin_th": 35.39, "self_any_th": 0, "self_end_th": 0, "end_stability": 3.86, "primer3_problems": null,
          "template": { "start": 441, "end": 464, "sequence": "ref" },
          "genomic": { "region": "1", "start": 11149, "end": 11172, "strand": -1, "blocks": [{ "start": 11149, "end": 11172 }] },
          "inserted_bases": 0,
          "deliberate_mismatch": null,
          "discrimination": null,
          "tailed": null,
          "neighbours": [
            { "key": "1:11161:A:T", "ids": ["rs872438201"], "label": "1:11161 A/T", "start": 11161, "end": 11161, "alleles": "A/T", "ems": false, "distance_from_3p": 13 }
          ]
        }
      },
      "products": {
        "ref": { "size": 92, "template": { "start": 373, "end": 464, "sequence": "ref" },
                 "genomic": { "region": "1", "start": 11081, "end": 11172, "strand": 1, "blocks": [{ "start": 11081, "end": 11172 }] }, "inserted_bases": 0 },
        "alt": { "size": 92, "template": { "start": 373, "end": 464, "sequence": "alt" },
                 "genomic": { "region": "1", "start": 11081, "end": 11172, "strand": 1, "blocks": [{ "start": 11081, "end": 11172 }] }, "inserted_bases": 0 }
      },
      "thermo": {
        "ref_common": { "compl_any_th": 0, "compl_end_th": 0 },
        "alt_common": { "compl_any_th": 0, "compl_end_th": 0 },
        "tailed": { "ref_alt_any_th": 5.46, "ref_alt_end_th": 0, "ref_common_any_th": 0, "ref_common_end_th": 0, "alt_common_any_th": 0, "alt_common_end_th": 0 }
      },
      "tm_balance": { "as_tm_diff": 0.95, "common_minus_as": -0.2 },
      "neighbour_sites": 1,
      "primer3_penalty": 14.6344,
      "warnings": [
        { "code": "ALT_PRIMER_SUBOPTIMAL", "message": "the derived as_alt primer would not have been chosen de novo (Primer3: temperature too low); it clears the hard floors",
          "details": { "oligo": "as_alt", "problems": " Temperature too low;" } },
        { "code": "TAILED_STRUCTURE", "message": "HEX-tailed as_alt: hairpin Tm 66.22 °C is at or above 55 °C, the final KASP annealing temperature",
          "details": { "oligo": "as_alt", "metric": "hairpin_th", "value": 66.22, "severity": "high" } },
        { "code": "NEIGHBOUR_IN_PRIMER", "message": "known variant rs872438201 lies in the common primer, 13 nt from its 3′ end",
          "details": { "role": "common", "ids": ["rs872438201"], "distances": [13] } }
      ],
      "issues": [
        { "code": "ALT_PRIMER_SUBOPTIMAL", "severity": "warn", "message": "as_alt: Primer3 reports temperature too low", "details": { "oligo": "as_alt", "problems": " Temperature too low;" } },
        { "code": "TAILED_STRUCTURE", "severity": "high", "message": "HEX-tailed as_alt hairpin 66.22 °C", "details": { "oligo": "as_alt", "metric": "hairpin_th", "value": 66.22 } },
        { "code": "NEIGHBOUR_IN_PRIMER", "severity": "warn", "message": "rs872438201 in the common primer, 13 nt from its 3′ end", "details": { "role": "common", "ids": ["rs872438201"], "distances": [13] } }
      ],
      "check": {
        "set": { "id": "S2", "ref_pair": "S2_REF", "alt_pair": "S2_ALT" },
        "pairs": [
          { "id": "S2_REF", "left": "GGTTATCCGAATATAGTCATACTCTATTC", "right": "TCTTTGTCTACTGAGAAATCCAGA", "expected": { "region": "1", "start": 11081, "end": 11172 } },
          { "id": "S2_ALT", "left": "GGTTATCCGAATATAGTCATACTCTATTA", "right": "TCTTTGTCTACTGAGAAATCCAGA", "expected": { "region": "1", "start": 11081, "end": 11172 } }
        ]
      },
      "order": [
        { "name": "rs871475760_S2_REF_FAM", "set_id": "S2", "set_key": "1accc54c262d", "role": "as_ref", "allele": "C", "dye": "FAM",
          "order_seq": "GAAGGTGACCAAGTTCATGCTGGTTATCCGAATATAGTCATACTCTATTC", "target_seq": "GGTTATCCGAATATAGTCATACTCTATTC", "tail_seq": "GAAGGTGACCAAGTTCATGCT",
          "length": 50, "tm": 57.28, "gc": 34.48, "orientation": "forward", "product_size_ref": 92, "product_size_alt": 92, "variant_key": "1:11109:C:A", "notes": "" },
        { "name": "rs871475760_S2_ALT_HEX", "set_id": "S2", "set_key": "1accc54c262d", "role": "as_alt", "allele": "A", "dye": "HEX",
          "order_seq": "GAAGGTCGGAGTCAACGGATTGGTTATCCGAATATAGTCATACTCTATTA", "target_seq": "GGTTATCCGAATATAGTCATACTCTATTA", "tail_seq": "GAAGGTCGGAGTCAACGGATT",
          "length": 50, "tm": 56.34, "gc": 31.03, "orientation": "forward", "product_size_ref": 92, "product_size_alt": 92, "variant_key": "1:11109:C:A", "notes": "tail hairpin 66.2 °C" },
        { "name": "rs871475760_S2_COM", "set_id": "S2", "set_key": "1accc54c262d", "role": "common", "allele": null, "dye": null,
          "order_seq": "TCTTTGTCTACTGAGAAATCCAGA", "target_seq": "TCTTTGTCTACTGAGAAATCCAGA", "tail_seq": null,
          "length": 24, "tm": 57.08, "gc": 37.5, "orientation": "forward", "product_size_ref": 92, "product_size_alt": 92, "variant_key": "1:11109:C:A", "notes": "" }
      ]
    }
  ],
  "check": {
    "request": {
      "system_name": "sorghum_bicolor",
      "mode": "region",
      "checks": ["specificity", "pangenome"],
      "pairs": [
        { "id": "S1_REF", "left": "AGCTTCTCTAAGTGGTTATCCGA", "right": "ATCTTTGACTAGCGAGAAATTCAG", "expected": { "region": "1", "start": 11068, "end": 11132 } },
        { "id": "S1_ALT", "left": "AGCTTCTCTAAGTGGTTATCCGA", "right": "ATCTTTGACTAGCGAGAAATTCAT", "expected": { "region": "1", "start": 11068, "end": 11132 } },
        { "id": "S2_REF", "left": "GGTTATCCGAATATAGTCATACTCTATTC", "right": "TCTTTGTCTACTGAGAAATCCAGA", "expected": { "region": "1", "start": 11081, "end": 11172 } },
        { "id": "S2_ALT", "left": "GGTTATCCGAATATAGTCATACTCTATTA", "right": "TCTTTGTCTACTGAGAAATCCAGA", "expected": { "region": "1", "start": 11081, "end": 11172 } }
      ],
      "genotyping": {
        "variant": { "region": "1", "position": 11109, "ref": "C", "alt": "A" },
        "sets": [
          { "id": "S1", "ref_pair": "S1_REF", "alt_pair": "S1_ALT" },
          { "id": "S2", "ref_pair": "S2_REF", "alt_pair": "S2_ALT" }
        ]
      }
    },
    "set_ids": ["S1", "S2"],
    "unique_primers": 6,
    "omitted_set_ids": []
  },
  "settings": {
    "preset": "kasp",
    "params": { "opt_size": 22, "min_size": 18, "max_size": 30, "opt_tm": 60, "min_tm": 57, "max_tm": 63, "min_gc": 30, "max_gc": 70,
                "max_tm_diff": 3, "max_poly_x": 5, "product_size_ranges": [[61, 120]] },
    "pinned": [],
    "ladder": [
      { "level": 1, "changes": { "max_size": 32, "min_tm": 55, "max_tm": 65, "min_gc": 20, "max_gc": 80, "product_size_ranges": [[65, 150]] } },
      { "level": 2, "changes": { "min_tm": 52, "max_tm_diff": 6 } }
    ],
    "floors": { "as_min_tm": 52, "as_min_gc": 15 }
  },
  "engine": { "primer3": "2.6.1", "thermo": "ntthal 2.6.1", "genotyping_design": "1", "variation_source": "ensembl 115" },
  "warnings": []
}
```

The client-side cost of `check.request` is 6 distinct primers over all 119 other assemblies: **2,690 CPU-s** against the 6,000 limit, computed with the server's own `cost.js` and the real genome sizes [V].

### 2.10 Further design examples (excerpts; `"…"` marks elided values)

**(a) tmp_1_11193_C_T — an EMS SNV needing relaxation in both orientations** [V]

```json
{ "system_name": "sorghum_bicolor", "variant": { "id": "tmp_1_11193_C_T" }, "assay": { "num_sets": 2 } }
```

```json
{
  "variant": { "key": "1:11193:C:T", "ids": ["tmp_1_11193_C_T"], "label": "1:11193 C/T", "kind": "snv", "ems": true,
               "zone": { "start": 11193, "end": 11193 },
               "records": [{ "id": "tmp_1_11193_C_T", "source": "EMS_PMID38100514_Jiao", "ems": true }], "…": "…" },
  "template": { "region": "1", "start": 10793, "end": 11593, "length": 801, "…": "…" },
  "assay": { "type": "kasp", "ems_target": true, "…": "…" },
  "neighbours": { "data": "ensembl", "window": { "start": 10793, "end": 11593 }, "variants": 54, "non_ems": 45, "ems": 9, "dense_non_ems": 1 },
  "orientations": {
    "forward": { "status": "ok", "reason": null, "discriminating_position": 11193, "relaxation_level": 1, "sets_found": 8, "blockers": [], "attempts": [
      { "level": 0, "changes": {},
        "explain": { "left":  { "raw": "considered 13, GC content failed 8, low tm 3, ok 2" },
                     "right": { "raw": "considered 4771, GC content failed 1241, low tm 1100, high tm 927, ok 1503" },
                     "pair":  { "raw": "considered 3006, unacceptable product size 3006, ok 0" } },
        "pairs_returned": 0, "not_scored": 0, "sets": 0 },
      { "level": 1, "changes": { "max_size": 32, "min_tm": 55, "max_tm": 65, "min_gc": 20, "max_gc": 80, "product_size_ranges": [[65, 150]] },
        "explain": { "left":  { "raw": "considered 15, low tm 9, ok 6" },
                     "right": { "raw": "considered 5490, GC content failed 293, low tm 1517, high tm 807, ok 2873" },
                     "pair":  { "raw": "considered 6436, unacceptable product size 6413, ok 23" } },
        "pairs_returned": 20, "not_scored": 12, "sets": 8 } ] },
    "reverse": { "status": "ok", "reason": null, "discriminating_position": 11193, "relaxation_level": 1, "sets_found": 8, "blockers": [], "attempts": [
      { "level": 0, "changes": {},
        "explain": { "left":  { "raw": "considered 4771, GC content failed 950, low tm 1792, high tm 507, ok 1522" },
                     "right": { "raw": "considered 13, GC content failed 13, ok 0" },
                     "pair":  { "raw": "considered 0, ok 0" } },
        "pairs_returned": 0, "not_scored": 0, "sets": 0 },
      { "level": 1, "changes": { "max_size": 32, "min_tm": 55, "max_tm": 65, "min_gc": 20, "max_gc": 80, "product_size_ranges": [[65, 150]] },
        "explain": { "left":  { "raw": "considered 5490, GC content failed 104, low tm 1994, high tm 435, ok 2957" },
                     "right": { "raw": "considered 15, GC content failed 1, low tm 12, ok 2" },
                     "pair":  { "raw": "considered 1364, unacceptable product size 1341, tm diff too large 2, ok 21" } },
        "pairs_returned": 20, "not_scored": 12, "sets": 8 } ] }
  },
  "sets": [
    { "id": "S1", "key": "f5a5c6f2ebcf", "rank": 0, "orientation": "forward", "relaxation_level": 1, "quality": "usable", "score": 22.51,
      "primers": {
        "as_ref": { "target_seq": "ACAAAGATAGATAACAAAAATAGCTCTC", "len": 28, "genomic": { "region": "1", "start": 11166, "end": 11193, "strand": 1 },
                    "tm": 56.43, "gc": 28.57, "hairpin_th": 0, "self_any_th": 0, "self_end_th": 0, "end_stability": 3.2,
                    "discrimination": { "own_allele": { "mm_pos": [], "likelihood": "likely" }, "other_allele": { "mm_pos": [1], "likelihood": "likely_weak" },
                                        "terminal_mismatch_class": "weak", "in_shift_tract": false },
                    "neighbours": [{ "key": "1:11182:A:G", "ids": ["rs873026643"], "ems": false, "distance_from_3p": 12 }], "…": "…" },
        "as_alt": { "target_seq": "ACAAAGATAGATAACAAAAATAGCTCTT", "len": 28, "tm": 56.07, "gc": 25, "end_stability": 2.85,
                    "discrimination": { "own_allele": { "mm_pos": [], "likelihood": "likely" }, "other_allele": { "mm_pos": [1], "likelihood": "likely_weak" },
                                        "terminal_mismatch_class": "weak", "in_shift_tract": false }, "…": "…" },
        "common": { "target_seq": "TCACACTTTGAATCATCATTTGGA", "genomic": { "region": "1", "start": 11268, "end": 11291, "strand": -1 },
                    "tm": 57.1, "gc": 33.33, "hairpin_th": 34.96, "end_stability": 3.53,
                    "neighbours": [{ "key": "1:11282:CA:C", "ids": ["rs5413864115"], "ems": false, "distance_from_3p": 16 }], "…": "…" } },
      "products": { "ref": { "size": 126 }, "alt": { "size": 126 } },
      "thermo": { "ref_common": { "compl_any_th": 5.73, "compl_end_th": 0 }, "alt_common": { "compl_any_th": 5.73, "compl_end_th": 0 }, "tailed": { "…": "…" } },
      "tm_balance": { "as_tm_diff": 0.36, "common_minus_as": 0.68 }, "neighbour_sites": 2, "primer3_penalty": 14.4705,
      "warnings": [ { "code": "NEIGHBOUR_IN_PRIMER", "details": { "role": "as_ref", "ids": ["rs873026643"], "distances": [12] }, "message": "…" },
                    { "code": "NEIGHBOUR_IN_PRIMER", "details": { "role": "common", "ids": ["rs5413864115"], "distances": [16] }, "message": "…" } ] },
    { "id": "S2", "key": "c791a4956fd3", "rank": 1, "orientation": "reverse", "relaxation_level": 1, "quality": "usable", "score": 24.84,
      "primers": {
        "as_ref": { "target_seq": "ACAACTTTTTAATATATTGTGTATACTCTAG", "len": 31, "genomic": { "region": "1", "start": 11193, "end": 11223, "strand": -1 },
                    "tm": 55.21, "gc": 22.58, "hairpin_th": 0, "self_any_th": 8.28, "self_end_th": 0, "end_stability": 2.43,
                    "discrimination": { "own_allele": { "mm_pos": [], "likelihood": "likely" }, "other_allele": { "mm_pos": [1], "likelihood": "likely_weak" },
                                        "terminal_mismatch_class": "weak", "in_shift_tract": false },
                    "neighbours": [{ "key": "1:11203:C:T", "ids": ["tmp_1_11203_C_T"], "ems": true, "distance_from_3p": 11 }], "…": "…" },
        "as_alt": { "target_seq": "ACAACTTTTTAATATATTGTGTATACTCTAA", "len": 31, "tm": 54.89, "gc": 19.35, "end_stability": 2.1, "primer3_problems": " GC content too low; Temperature too low;",
                    "discrimination": { "own_allele": { "mm_pos": [], "likelihood": "likely" }, "other_allele": { "mm_pos": [1], "likelihood": "likely_weak" },
                                        "terminal_mismatch_class": "weak", "in_shift_tract": false }, "…": "…" },
        "common": { "target_seq": "TGAATTTCTCGCTAGTCAAAGA", "genomic": { "region": "1", "start": 11110, "end": 11131, "strand": 1 },
                    "tm": 55.51, "gc": 36.36, "hairpin_th": 0, "end_stability": 2.52,
                    "neighbours": [], "…": "…" } },
      "products": { "ref": { "size": 114 }, "alt": { "size": 114 } },
      "thermo": { "ref_common": { "compl_any_th": 0, "compl_end_th": 0 }, "alt_common": { "compl_any_th": 0, "compl_end_th": 0 }, "tailed": { "…": "…" } },
      "tm_balance": { "as_tm_diff": 0.33, "common_minus_as": 0.3 }, "neighbour_sites": 0, "primer3_penalty": 18.2784,
      "warnings": [ { "code": "ALT_PRIMER_SUBOPTIMAL", "details": { "oligo": "as_alt", "problems": " GC content too low; Temperature too low;" }, "message": "…" } ] }
  ],
  "warnings": [
    { "code": "EMS_TARGET", "message": "the target is an EMS mutation, private to BTx623-background mutant lines; natural neighbours are reported as warnings rather than blocking an orientation", "details": { "sources": ["EMS_PMID38100514_Jiao"] } },
    { "code": "RELAXED_CONSTRAINTS", "message": "forward: sets need relaxation level 1", "details": { "orientation": "forward", "level": 1, "changes": { "…": "…" } } },
    { "code": "RELAXED_CONSTRAINTS", "message": "reverse: sets need relaxation level 1", "details": { "orientation": "reverse", "level": 1, "changes": { "…": "…" } } }
  ]
}
```

The ALT primer of S2 is 19.35 % GC. It stays because the floor for a *derived* ALT primer is `as_min_gc` 15 % and its Tm (54.89 °C) clears the 52 °C floor; Primer3's own `problems` string is echoed so the user can see why the ALT oligo would not have been designed de novo (§4.6).

**(b) rs5413864115 — a shiftable deletion entered manually in Ensembl style; ids filled from Ensembl** [V]

```json
{ "system_name": "sorghum_bicolor", "variant": { "region": "1", "position": 11283, "ref": "A", "alt": "-" }, "assay": { "num_sets": 2 } }
```

```json
{
  "variant": { "key": "1:11282:CA:C", "requested_id": null, "ids": ["rs5413864115"], "label": "1:11283-11283 A/-", "kind": "deletion",
               "vcf": { "position": 11282, "ref": "CA", "alt": "C" }, "minimal": { "start": 11283, "end": 11283, "ref": "A", "alt": "-" },
               "shift": 2, "zone": { "start": 11282, "end": 11286 },
               "discriminating": { "forward": { "position": 11285, "ref_base": "A", "alt_base": "G", "alt_maps_to": 11286 },
                                   "reverse": { "position": 11283, "ref_base": "A", "alt_base": "C", "alt_maps_to": 11282 } }, "…": "…" },
  "template": { "start": 10883, "end": 11685, "length": 803, "alt_length": 802,
                "features": { "variant": { "start": 400, "end": 401 }, "zone": { "start": 400, "end": 404 },
                              "discriminating": { "forward": 403, "reverse": 401 }, "alt_offset": -1, "exempt": [[367, 38], [400, 38]] }, "…": "…" },
  "sets": [
    { "id": "S1", "key": "a5277232d8ab", "orientation": "forward", "relaxation_level": 1, "quality": "usable", "score": 15.48,
      "primers": {
        "as_ref": { "target_seq": "ACAGATGATTTTCCAAATGATGATTCAAA", "len": 29, "tm": 59.22, "gc": 27.59, "end_stability": 2.69,
                    "genomic": { "region": "1", "start": 11257, "end": 11285, "strand": 1, "blocks": [{ "start": 11257, "end": 11285 }] },
                    "discrimination": { "own_allele": { "mm_pos": [], "likelihood": "likely" }, "other_allele": { "mm_pos": [1], "likelihood": "likely_weak" },
                                        "terminal_mismatch_class": "weak", "in_shift_tract": true }, "…": "…" },
        "as_alt": { "target_seq": "ACAGATGATTTTCCAAATGATGATTCAAG", "len": 29, "tm": 59.53, "gc": 31.03, "end_stability": 3.02,
                    "genomic": { "region": "1", "start": 11257, "end": 11286, "strand": 1, "blocks": [{ "start": 11257, "end": 11282 }, { "start": 11284, "end": 11286 }] },
                    "inserted_bases": 0, "…": "…" },
        "common": { "target_seq": "CCCCATGTTTTTGTTCCTTCCA", "genomic": { "region": "1", "start": 11323, "end": 11344, "strand": -1 }, "tm": 59.3,
                    "neighbours": [ { "key": "1:11334:A:G", "ids": ["tmp_1_11334_A_G"], "ems": false, "distance_from_3p": 12 },
                                    { "key": "1:11335:A:G", "ids": ["tmp_1_11335_A_G"], "ems": false, "distance_from_3p": 13 },
                                    { "key": "1:11340:TG:T", "ids": ["tmp_1_11340_TG_T"], "ems": false, "distance_from_3p": 19 } ], "…": "…" } },
      "products": { "ref": { "size": 88 }, "alt": { "size": 87 } }, "neighbour_sites": 3,
      "warnings": [ { "code": "NEIGHBOUR_IN_PRIMER", "details": { "role": "common", "ids": ["tmp_1_11334_A_G", "tmp_1_11335_A_G", "tmp_1_11340_TG_T"], "distances": [12, 13, 19] }, "message": "…" },
                    { "code": "SHIFT_TRACT_DISCRIMINATION", "details": { "primer": "as_ref", "shift": 2 }, "message": "…" } ], "…": "…" },
    { "id": "S2", "key": "cb3ef66afd37", "orientation": "reverse", "relaxation_level": 1, "quality": "usable", "score": 23.95,
      "primers": {
        "as_ref": { "target_seq": "AGAGTCTTTTCAAATTTCACACTTT", "len": 25, "tm": 56.09,
                    "genomic": { "region": "1", "start": 11283, "end": 11307, "strand": -1 }, "…": "…" },
        "as_alt": { "target_seq": "AGAGTCTTTTCAAATTTCACACTTG", "len": 25, "tm": 56.71,
                    "genomic": { "region": "1", "start": 11282, "end": 11307, "strand": -1, "blocks": [{ "start": 11282, "end": 11282 }, { "start": 11284, "end": 11307 }] }, "…": "…" },
        "common": { "target_seq": "ACAAAAATAGCTCTCTAGAGTATACACA", "genomic": { "region": "1", "start": 11179, "end": 11206, "strand": 1 }, "tm": 58.12, "…": "…" } },
      "products": { "ref": { "size": 129 }, "alt": { "size": 128 } }, "neighbour_sites": 3,
      "warnings": [ { "code": "NEIGHBOUR_IN_PRIMER", "details": { "role": "as_ref", "ids": ["tmp_1_11298_A_G", "tmp_1_11303_A_T"], "distances": [16, 21] }, "message": "…" },
                    { "code": "NEIGHBOUR_IN_PRIMER", "details": { "role": "common", "ids": ["rs873026643"], "distances": [25] }, "message": "…" },
                    { "code": "SHIFT_TRACT_DISCRIMINATION", "details": { "primer": "as_ref", "shift": 2 }, "message": "…" } ], "…": "…" }
  ],
  "warnings": [
    { "code": "SHIFTABLE_INDEL", "message": "the deletion can slide 2 bases inside AAA; the forward and reverse primers end at different bases and the wrong-allele primer may still prime through a 1-nt bulge",
      "details": { "shift": 2, "forward_position": 11285, "reverse_position": 11283, "zone": { "start": 11282, "end": 11286 } } },
    { "code": "DENSE_NEIGHBOURS", "message": "4 known non-EMS variants lie within 30 bp of the variant", "details": { "count": 4, "window": 30 } },
    { "code": "RELAXED_CONSTRAINTS", "details": { "orientation": "forward", "level": 1 }, "message": "…" },
    { "code": "RELAXED_CONSTRAINTS", "details": { "orientation": "reverse", "level": 1 }, "message": "…" }
  ]
}
```

Both ALT products differ from their REF products by `alt_offset` (−1): 88 → 87 and 129 → 128, confirmed by `check_primers` on the ALT haplotype [V].

**(c) tmp_1_11502_C_CGT — an insertion by id; duplicate ids merged; the forward orientation blocked** [V]

```json
{ "system_name": "sorghum_bicolor", "variant": { "id": "tmp_1_11502_C_CGT" }, "assay": { "num_sets": 1 } }
```

```json
{
  "variant": { "key": "1:11502:C:CGT", "requested_id": "tmp_1_11502_C_CGT", "ids": ["tmp_1_11502_C_CGT", "rs5413863549"],
               "label": "1:11502^11503 -/GT", "kind": "insertion", "shift": 0, "zone": { "start": 11502, "end": 11503 },
               "discriminating": { "forward": { "position": 11503, "ref_base": "A", "alt_base": "G", "alt_maps_to": null },
                                   "reverse": { "position": 11502, "ref_base": "C", "alt_base": "T", "alt_maps_to": null } }, "…": "…" },
  "template": { "start": 11102, "end": 11903, "length": 802, "alt_length": 804, "…": "…" },
  "orientations": {
    "forward": { "status": "blocked", "reason": "neighbour_at_3p", "discriminating_position": 11503, "relaxation_level": null, "sets_found": 0,
                 "blockers": [{ "key": "1:11500:G:A", "ids": ["rs5413863234"], "label": "1:11500 G/A", "start": 11500, "end": 11500, "alleles": "G/A", "ems": false, "distance_from_3p": 4 }],
                 "attempts": [] },
    "reverse": { "status": "ok", "reason": null, "discriminating_position": 11502, "relaxation_level": 0, "sets_found": 8, "blockers": [], "…": "…" }
  },
  "sets": [
    { "id": "S1", "key": "53942cb55348", "orientation": "reverse", "relaxation_level": 0, "quality": "usable", "score": 9.56,
      "primers": {
        "as_ref": { "allele": "C", "three_prime_base": "G", "target_seq": "GCAGGAAAAGAAATCCTAACATCATATG", "len": 28,
                    "genomic": { "region": "1", "start": 11502, "end": 11529, "strand": -1 }, "tm": 59.19, "gc": 35.71, "hairpin_th": 45.59, "end_stability": 1.78,
                    "neighbours": [ { "key": "1:11509:T:C", "ids": ["rs5413864238"], "ems": false, "distance_from_3p": 8 },
                                    { "key": "1:11516:A:T", "ids": ["rs5413863874"], "ems": false, "distance_from_3p": 15 },
                                    { "key": "1:11523:T:A", "ids": ["rs5413863270"], "ems": false, "distance_from_3p": 22 } ], "…": "…" },
        "as_alt": { "allele": "CGT", "three_prime_base": "A", "target_seq": "GCAGGAAAAGAAATCCTAACATCATATA", "len": 28,
                    "genomic": { "region": "1", "start": 11503, "end": 11529, "strand": -1, "blocks": [{ "start": 11503, "end": 11529 }] },
                    "inserted_bases": 1, "tm": 58.03, "gc": 32.14, "hairpin_th": 45.59, "end_stability": 0.86, "…": "…" },
        "common": { "target_seq": "AGGATCTTTGCAACCCTGTGTT", "genomic": { "region": "1", "start": 11458, "end": 11479, "strand": 1 },
                    "tm": 60.43, "gc": 45.45, "end_stability": 3.32,
                    "neighbours": [ { "key": "1:11474:T:A", "ids": ["rs5413863452"], "ems": false, "distance_from_3p": 6 } ], "…": "…" } },
      "products": { "ref": { "size": 72 }, "alt": { "size": 74, "inserted_bases": 2 } },
      "tm_balance": { "as_tm_diff": 1.16, "common_minus_as": 1.24 },
      "neighbour_sites": 4,
      "warnings": [
        { "code": "AS_TM_IMBALANCE", "message": "as_ref and as_alt Tm differ by 1.16 °C (limit 1.0)", "details": { "diff": 1.16 } },
        { "code": "NEIGHBOUR_IN_PRIMER", "details": { "role": "as_ref", "ids": ["rs5413864238", "rs5413863874", "rs5413863270"], "distances": [8, 15, 22] }, "message": "…" },
        { "code": "NEIGHBOUR_IN_PRIMER", "details": { "role": "common", "ids": ["rs5413863452"], "distances": [6] }, "message": "…" }
      ], "…": "…" }
  ],
  "warnings": [
    { "code": "DUPLICATE_VARIANT_IDS", "details": { "key": "1:11502:C:CGT", "ids": ["tmp_1_11502_C_CGT", "rs5413863549"] }, "message": "…" },
    { "code": "ORIENTATION_BLOCKED", "message": "forward: known variant rs5413863234 (G/A) lies 4 nt from the allele-specific primer's 3′ end",
      "details": { "orientation": "forward", "ids": ["rs5413863234"], "distances": [4] } },
    { "code": "DENSE_NEIGHBOURS", "details": { "count": 5, "window": 30 }, "message": "…" }
  ]
}
```

For an insertion the ALT primer's `allele` is `vcf.alt` (`CGT`), because alleles are written in VCF form throughout.

**(d) rs871475760 as gel AS-PCR with a deliberate mismatch at −2** [V]

```json
{ "system_name": "sorghum_bicolor", "variant": { "id": "rs871475760" }, "assay": { "type": "as_pcr", "num_sets": 1 } }
```

```json
{
  "assay": { "type": "as_pcr", "orientation": "both", "tails": "none", "deliberate_mismatch": "auto", "mismatch_position": 2,
             "num_sets": 1, "max_relaxation": 2, "neighbour_policy": "avoid_3p", "ems_target": false, "kasp_mix": null },
  "sets": [
    { "id": "S1", "key": "7f9af6b1c938", "orientation": "reverse", "relaxation_level": 0, "quality": "usable", "score": 7.92,
      "primers": {
        "as_ref": { "role": "as_ref", "allele": "C", "three_prime_base": "G",
                    "target_seq": "ATCTTTGACTAGCGAGAAATTCGG", "matched_seq": "ATCTTTGACTAGCGAGAAATTCAG",
                    "dye": null, "tail_seq": null, "order_seq": "ATCTTTGACTAGCGAGAAATTCGG",
                    "tm": 54.69, "tm_method": "ntthal_duplex", "matched_tm": 57.1, "gc": 41.67,
                    "hairpin_th": 49.2, "self_any_th": 0, "self_end_th": 0, "end_stability": 4.3,
                    "primer3_problems": " Hairpin stability too high;",
                    "deliberate_mismatch": { "position": 2, "original_base": "A", "new_base": "G", "template_base": "T",
                                             "terminal_pair": "AG", "terminal_mismatch_class": "max", "added_mismatch_class": "weak",
                                             "source": "Little 1995 Table 9.8.1" },
                    "discrimination": { "own_allele": { "mm_pos": [2], "likelihood": "likely" },
                                        "other_allele": { "mm_pos": [2, 1], "likelihood": "unlikely" },
                                        "terminal_mismatch_class": "max", "in_shift_tract": false }, "…": "…" },
        "as_alt": { "role": "as_alt", "allele": "A", "three_prime_base": "T",
                    "target_seq": "ATCTTTGACTAGCGAGAAATTCGT", "matched_seq": "ATCTTTGACTAGCGAGAAATTCAT", "order_seq": "ATCTTTGACTAGCGAGAAATTCGT",
                    "tm": 54.58, "tm_method": "ntthal_duplex", "matched_tm": 56.51, "gc": 37.5, "hairpin_th": 50.62, "end_stability": 3.85,
                    "primer3_problems": " Hairpin stability too high;",
                    "deliberate_mismatch": { "position": 2, "original_base": "A", "new_base": "G", "template_base": "T",
                                             "terminal_pair": "CT", "terminal_mismatch_class": "max", "added_mismatch_class": "weak",
                                             "source": "Little 1995 Table 9.8.1" },
                    "discrimination": { "own_allele": { "mm_pos": [2], "likelihood": "likely" },
                                        "other_allele": { "mm_pos": [2, 1], "likelihood": "unlikely" }, "…": "…" }, "…": "…" },
        "common": { "target_seq": "TGCATCAACAAATGTGCTATGTGT", "genomic": { "region": "1", "start": 10880, "end": 10903, "strand": 1 }, "tm": 60.02,
                    "neighbours": [ { "key": "1:10902:G:A", "ids": ["tmp_1_10902_G_A"], "ems": true, "distance_from_3p": 2 },
                                    { "key": "1:10880:T:C", "ids": ["rs5413864006"], "ems": false, "distance_from_3p": 24 } ], "…": "…" } },
      "products": { "ref": { "size": 253 }, "alt": { "size": 253 } },
      "tm_balance": { "as_tm_diff": 0.11, "common_minus_as": 2.92 },
      "warnings": [
        { "code": "ALT_PRIMER_SUBOPTIMAL", "details": { "oligo": "as_alt", "problems": " Hairpin stability too high;" }, "message": "…" },
        { "code": "ALT_PRIMER_SUBOPTIMAL", "details": { "oligo": "as_ref", "problems": " Hairpin stability too high;" }, "message": "…" },
        { "code": "MISMATCH_STRUCTURE", "details": { "oligo": "as_ref", "metric": "hairpin_th", "value": 49.2, "severity": "warn" }, "message": "…" },
        { "code": "MISMATCH_STRUCTURE", "details": { "oligo": "as_alt", "metric": "hairpin_th", "value": 50.62, "severity": "warn" }, "message": "…" },
        { "code": "NEIGHBOUR_IN_PRIMER", "details": { "role": "common", "ids": ["rs5413864006"], "distances": [24] }, "message": "…" }
      ],
      "order": [
        { "name": "rs871475760_S1_REF", "role": "as_ref", "allele": "C", "dye": null, "order_seq": "ATCTTTGACTAGCGAGAAATTCGG", "length": 24, "tm": 54.69,
          "notes": "deliberate mismatch −2 A→G (Little 1995 Table 9.8.1)", "…": "…" },
        { "name": "rs871475760_S1_ALT", "role": "as_alt", "allele": "A", "dye": null, "order_seq": "ATCTTTGACTAGCGAGAAATTCGT", "length": 24, "tm": 54.58,
          "notes": "deliberate mismatch −2 A→G (Little 1995 Table 9.8.1)", "…": "…" },
        { "name": "rs871475760_S1_COM", "role": "common", "allele": null, "dye": null, "order_seq": "TGCATCAACAAATGTGCTATGTGT", "length": 24, "tm": 60.02, "notes": "", "…": "…" }
      ], "…": "…" }
  ],
  "…": "…"
}
```

The EMS neighbour `tmp_1_10902_G_A` sits 2 nt from the common primer's 3′ end and the set is still returned: EMS alleles are private to mutant lines in the BTx623 background, so they do not segregate in the panels a breeder genotypes (§4.15). `common_minus_as` compares the common primer's Tm with the higher **nominal** (perfect-match) Tm of the two allele-specific primers, 57.10 °C, not with their mismatch duplex Tm, so the difference is 2.92 °C and inside the −1…+3 °C window.

**(e) `template_only` with a wrong reference allele**

```json
{ "system_name": "sorghum_bicolor", "variant": { "region": "1", "position": 11109, "ref": "A", "alt": "C" }, "template_only": true }
```

```json
{ "message": "the reference allele A does not match the genome base C at 1:11109", "code": "REF_MISMATCH",
  "details": { "region": "1", "position": 11109, "given": "A", "genome": "C" } }
```

**(f) Other error bodies**

```json
{ "message": "rs5413863494 has more than one alternative allele; choose one with variant.alt", "code": "ALT_REQUIRED",
  "details": { "id": "rs5413863494", "alts": ["T", "G"] } }
```

```json
{ "message": "unknown variant rs0000000001 for sorghum_bicolor", "code": "UNKNOWN_VARIANT",
  "details": { "id": "rs0000000001", "system_name": "sorghum_bicolor" } }
```

```json
{ "message": "sorghum_rio has no known-variant data; enter the variant as region, position, ref and alt", "code": "NO_VARIATION_DATA",
  "details": { "system_name": "sorghum_rio" } }
```

```json
{ "message": "the Ensembl variation service is temporarily unavailable", "code": "VARIATION_SOURCE_UNAVAILABLE",
  "details": { "retry_after_s": 30, "reason": "timeout" } }
```
### 2.11 `POST /primers/check`: the `genotyping` block

Add to `PrimerCheckRequest.properties` (`swagger.yaml:2154`; the definition stays `additionalProperties: false`):

```yaml
genotyping: { $ref: "#/definitions/PrimerCheckGenotyping" }
```

Add to `definitions:`:

```yaml
  PrimerCheckGenotypingVariant:
    description: "the design response's variant.vcf plus its region; ACGT alleles only, left-aligned or not (the server normalizes)"
    type: object
    required: [region, position, ref, alt]
    additionalProperties: false
    properties:
      region: { type: string, maxLength: 255 }
      position: { type: integer, minimum: 1 }
      ref: { type: string, pattern: "^[ACGTacgt]{1,50}$" }
      alt: { type: string, pattern: "^[ACGTacgt]{1,50}$" }
  PrimerCheckGenotypingSet:
    type: object
    required: [id, ref_pair, alt_pair]
    additionalProperties: false
    properties:
      id: { type: string, pattern: "^[A-Za-z0-9_.:-]{1,64}$" }
      ref_pair: { type: string, pattern: "^[A-Za-z0-9_.:-]{1,64}$", description: "id of the pairs[] item whose allele-specific primer carries the REF allele" }
      alt_pair: { type: string, pattern: "^[A-Za-z0-9_.:-]{1,64}$" }
  PrimerCheckGenotyping:
    description: "Allele-specific sets over ordinary pairs. gene and region modes only; every referenced pair needs expected. The server derives each set's orientation and allele side from the shared primer and verifies both allele-specific primers against the reference haplotypes at submit."
    type: object
    required: [variant, sets]
    additionalProperties: false
    properties:
      variant: { $ref: "#/definitions/PrimerCheckGenotypingVariant" }
      sets:
        type: array
        minItems: 1
        maxItems: 5
        items: { $ref: "#/definitions/PrimerCheckGenotypingSet" }
```

All three new definitions join the strict-definition list of `contract.test.js:118-132`.

**Example request** — `design.check.request` of §2.9 narrowed to three genomes:

```json
{
  "system_name": "sorghum_bicolor",
  "mode": "region",
  "checks": ["specificity", "pangenome"],
  "genomes": ["sorghum_bicolorv5", "sorghum_pi180348", "sorghum_pi329250"],
  "pairs": [
    { "id": "S1_REF", "left": "AGCTTCTCTAAGTGGTTATCCGA", "right": "ATCTTTGACTAGCGAGAAATTCAG", "expected": { "region": "1", "start": 11068, "end": 11132 } },
    { "id": "S1_ALT", "left": "AGCTTCTCTAAGTGGTTATCCGA", "right": "ATCTTTGACTAGCGAGAAATTCAT", "expected": { "region": "1", "start": 11068, "end": 11132 } },
    { "id": "S2_REF", "left": "GGTTATCCGAATATAGTCATACTCTATTC", "right": "TCTTTGTCTACTGAGAAATCCAGA", "expected": { "region": "1", "start": 11081, "end": 11172 } },
    { "id": "S2_ALT", "left": "GGTTATCCGAATATAGTCATACTCTATTA", "right": "TCTTTGTCTACTGAGAAATCCAGA", "expected": { "region": "1", "start": 11081, "end": 11172 } }
  ],
  "genotyping": {
    "variant": { "region": "1", "position": 11109, "ref": "C", "alt": "A" },
    "sets": [
      { "id": "S1", "ref_pair": "S1_REF", "alt_pair": "S1_ALT" },
      { "id": "S2", "ref_pair": "S2_REF", "alt_pair": "S2_ALT" }
    ]
  }
}
```

**Response `202`** (the `job_id` and timestamps are illustrative but syntactically valid; the estimate is real: 6 distinct primers over 3 genomes = 95 CPU-s [V]):

```json
{ "job_id": "5d0c3b9e2a7f41c68e1f9a2b7c4d6e80", "status": "queued", "kind": "pangenome", "queue_position": 0,
  "progress": { "done": 0, "total": 4, "stage": "queued", "running": [] },
  "estimate": { "cpu_s": 95 }, "created_at": "2026-09-14T23:30:00.000Z", "warnings": [] }
```

**Submit-time rejection** (one code, many reasons — §5.2):

```json
{ "message": "set S1: S1_REF and S1_ALT must share exactly one primer, the common primer, on the same side",
  "code": "GENOTYPING_SET_INVALID",
  "details": { "set_id": "S1", "reason": "no_shared_common", "pair_ids": ["S1_REF", "S1_ALT"] } }
```

```json
{ "message": "set S1: the REF-specific primer TCTTTGACTAGCGAGAAATTCAGA does not end at the variant (1:11109)",
  "code": "GENOTYPING_SET_INVALID",
  "details": { "set_id": "S1", "reason": "not_at_variant", "pair_ids": ["S1_REF"], "primer": "ref_pair",
               "sequence": "TCTTTGACTAGCGAGAAATTCAGA" } }
```

### 2.12 `GET /primers/check/{job_id}`: `results.genotyping`

```yaml
  PrimerGenotypeVariantInfo:
    type: object
    properties:
      key: { type: string }
      region: { type: string }
      position: { type: integer }
      ref: { type: string }
      alt: { type: string }
      shift: { type: integer, description: "with shift > 0 the indel has the shift tract [position + 1, position + len(ref) - 1 + shift] of the left-aligned variant (section 3.5)" }
      zone: { type: object, properties: { start: { type: integer }, end: { type: integer } } }
      flank: { type: integer, description: "K = max(15, zone length + max(len(ref), len(alt))), the flank on each side of the core" }
      core: { type: object, properties: { ref: { type: string }, alt: { type: string } }, description: "the reference span [position - 1, position + len(ref) - 1 + shift + 1], i.e. the left-aligned VCF span with its anchor base, extended right by shift, plus one base on each side (not the shift tract), read on each haplotype (the ALT core carries the inserted bases or lacks the deleted ones); a genome is called ref or alt only on an exact match (section 5.6)" }
      haplotypes: { type: object, properties: { ref: { type: string }, alt: { type: string } }, description: "core plus K reference bases of flank on each side" }
  PrimerGenotypeCopy:
    type: object
    properties:
      region: { type: string }
      start: { type: integer, description: "envelope of the merged anchor amplicons, or the megablast HSP" }
      end: { type: integer }
      strand: { type: integer, enum: [1, -1] }
      variant_position: { type: integer, x-nullable: true, description: "genome coordinate aligned to vcf.position" }
      identity: { type: number, description: "percent identity over the longest anchor alignment (the anchor amplicon +-50 bp with the most alignment columns): 100 x (columns - edits) / columns, 2 decimals" }
      gap_compressed_identity: { type: number, description: "the same alignment with each gap run counted as one edit and one column; used for the ortholog test" }
      aligned_length: { type: integer, description: "columns of that longest anchor alignment, gap columns included" }
      observed: { type: string, x-nullable: true, description: "the genome's core sequence, always reported, including when the call is other" }
      flank_edits: { type: integer, x-nullable: true, description: "edits outside the core; confidence only, never part of the call" }
      call: { type: string, enum: [ref, alt, other, missing] }
      anchors: { type: integer }
      ortholog: { type: boolean, x-nullable: true, description: "annotation-based in gene mode; null in region mode" }
      source: { type: string, enum: [amplicon, megablast] }
  PrimerGenotypeGenome:
    type: object
    properties:
      system_name: { type: string }
      display_name: { type: string }
      is_reference: { type: boolean }
      allele: { type: string, enum: [ref, alt, other, ambiguous, missing, unavailable] }
      observed: { type: string, x-nullable: true, description: "core shared by all orthologous copies; null when they differ" }
      source: { type: string, enum: [amplicon, megablast], x-nullable: true }
      copies: { type: array, items: { $ref: "#/definitions/PrimerGenotypeCopy" } }
      orthologous_copies: { type: integer }
      paralog_copies: { type: integer }
      reason:
        type: string
        x-nullable: true
        enum: [variant_not_covered, no_orthologous_copy, fallback_failed, fallback_budget, db_unavailable, blast_error, call_failed]
  PrimerGenotypePrimerCall:
    type: object
    description: "status is the best over the primer's on-locus products and copies (section 5.7); likelihood, mm_pos and residual_mm_pos all come from the one product that supplied that status, and are null with no_product"
    properties:
      status: { type: string, enum: [match, weak, terminal_mismatch, uncertain, blocked, no_product, unknown] }
      likelihood: { type: string, enum: [likely, likely_weak, unlikely], x-nullable: true, description: "classify likelihood of that product; for common_primer, the product of either pair whose common site gave status, the better likelihood first when several tie" }
      mm_pos: { type: array, items: { type: integer }, x-nullable: true, description: "this primer's site alignment in that product; null with no_product or an approximate alignment" }
      residual_mm_pos: { type: array, items: { type: integer }, x-nullable: true, description: "mm_pos with the set's declared deliberate-mismatch positions removed; the prediction reads this" }
  PrimerGenotypeSetGenome:
    type: object
    description: "ref_primer and alt_primer come from their own pair's on-locus product; common_primer comes from the common primer's own site at the called locus, the better of its two alignments in the pairs' on-locus products. A deliberate mismatch on an allele-specific primer never changes common_primer. Each call reports the likelihood and mm_pos of the product that supplied its status (section 5.7)."
    properties:
      system_name: { type: string }
      ref_primer: { $ref: "#/definitions/PrimerGenotypePrimerCall" }
      alt_primer: { $ref: "#/definitions/PrimerGenotypePrimerCall" }
      common_primer: { $ref: "#/definitions/PrimerGenotypePrimerCall" }
      predicted: { type: string, enum: [ref, alt, both, none, no_call, unknown] }
      strength: { type: string, enum: [normal, weak], x-nullable: true, description: "weak when the amplifying side, or the common primer, has a residual mismatch at distance 2-3" }
      agrees: { type: boolean, x-nullable: true }
      reasons:
        type: array
        items:
          type: string
          enum: [common_primer_3p_mismatch, common_primer_weak, ref_signal_off_locus, alt_signal_off_locus, shift_tract_uncertain, approx_alignment, no_orthologous_copy, third_allele]
      off_locus_products: { type: integer }
  PrimerGenotypeSetSummary:
    type: object
    properties:
      genomes_total: { type: integer }
      predicted_ref: { type: integer }
      predicted_alt: { type: integer }
      both: { type: integer }
      none: { type: integer }
      no_call: { type: integer }
      unknown: { type: integer }
      weak: { type: integer }
      agree: { type: integer }
      disagree: { type: integer }
      not_comparable: { type: integer }
  PrimerGenotypePairSpecificity:
    type: object
    properties:
      verdict: { type: string }
      consistent_with_allele: { type: boolean }
  PrimerGenotypeReferenceControl:
    type: object
    properties:
      status: { type: string, enum: [pass, warn, fail] }
      allele: { type: string, enum: [ref, alt, other, ambiguous, missing, unavailable] }
      reasons: { type: array, items: { type: string } }
  PrimerGenotypeSetResult:
    type: object
    properties:
      id: { type: string }
      ref_pair: { type: string }
      alt_pair: { type: string }
      orientation: { type: string, enum: [forward, reverse], description: "derived at submit from the shared primer, not taken from the client" }
      deliberate_mismatch_positions: { type: array, items: { type: integer }, description: "positions accepted at submit on each allele-specific primer (2 or 3), removed before prediction" }
      specificity:
        type: object
        properties:
          ref_pair: { $ref: "#/definitions/PrimerGenotypePairSpecificity" }
          alt_pair: { $ref: "#/definitions/PrimerGenotypePairSpecificity" }
          off_target_count: { type: integer }
      control: { $ref: "#/definitions/PrimerGenotypeReferenceControl" }
      reference: { $ref: "#/definitions/PrimerGenotypeSetGenome" }
      summary: { $ref: "#/definitions/PrimerGenotypeSetSummary" }
      genomes: { type: array, items: { $ref: "#/definitions/PrimerGenotypeSetGenome" } }
  PrimerGenotypeSummary:
    type: object
    description: "pan-genome genomes only; the counts add up to genomes_total"
    properties:
      genomes_total: { type: integer }
      ref: { type: integer }
      alt: { type: integer }
      other: { type: integer }
      ambiguous: { type: integer }
      missing: { type: integer }
      unavailable: { type: integer }
  PrimerCheckGenotypingResults:
    type: object
    description: "present only for requests that carried a genotyping block; absent (not null) otherwise"
    properties:
      algorithm_version: { type: string, description: "caller and prediction version, currently 'g1'" }
      variant: { $ref: "#/definitions/PrimerGenotypeVariantInfo" }
      summary: { $ref: "#/definitions/PrimerGenotypeSummary" }
      genomes: { type: array, items: { $ref: "#/definitions/PrimerGenotypeGenome" }, description: "reference control first, then finished pan-genome genomes in request order" }
      sets: { type: array, items: { $ref: "#/definitions/PrimerGenotypeSetResult" } }
```

Add to `PrimerCheckResults.properties` (`swagger.yaml:2395`):

```yaml
genotyping: { $ref: "#/definitions/PrimerCheckGenotypingResults" }
```

### 2.13 Complete `results.genotyping` example

For the §2.11 request. **What is verified:** the copy coordinates, identities, observed cores, calls, mismatch positions and per-primer statuses are from alignments of the designed primers and the reference amplicon ±50 bp against the real assembly FASTAs at the research megablast loci, classified with the worktree's own `realign.js` and `classify.js` [V, research-derived]. Each copy's `identity`, `gap_compressed_identity` and `aligned_length` are those of its longest anchor alignment, set S2's `R_s` 1:11031–11222 (§5.6 step 5; recomputed by `scratch-fix/predict.js` part C). The surrounding `specificity`, `pangenome`, `primers` and `timings_ms` blocks of the same job are **illustrative** [I] — no job was run while writing this spec — and §7.5 requires them to be replaced by a real 50112 run before the fixture is handed to the UI session.

```json
{
  "results": {
    "genotyping": {
      "algorithm_version": "g1",
      "variant": {
        "key": "1:11109:C:A", "region": "1", "position": 11109, "ref": "C", "alt": "A", "shift": 0,
        "zone": { "start": 11109, "end": 11109 }, "flank": 15,
        "core": { "ref": "TCT", "alt": "TAT" },
        "haplotypes": { "ref": "ATAGTCATACTCTATTCTGAATTTCTCGCTAGT", "alt": "ATAGTCATACTCTATTATGAATTTCTCGCTAGT" }
      },
      "summary": { "genomes_total": 3, "ref": 2, "alt": 1, "other": 0, "ambiguous": 0, "missing": 0, "unavailable": 0 },
      "genomes": [
        { "system_name": "sorghum_bicolor", "display_name": "Sb bicolor BTx623 v3", "is_reference": true,
          "allele": "ref", "observed": "TCT", "source": "amplicon",
          "copies": [ { "region": "1", "start": 11068, "end": 11172, "strand": 1, "variant_position": 11109,
                        "identity": 100, "gap_compressed_identity": 100, "aligned_length": 192, "observed": "TCT",
                        "flank_edits": 0, "call": "ref", "anchors": 2, "ortholog": null, "source": "amplicon" } ],
          "orthologous_copies": 1, "paralog_copies": 0, "reason": null },
        { "system_name": "sorghum_bicolorv5", "display_name": "Sb bicolor BTx623 v5", "is_reference": false,
          "allele": "ref", "observed": "TCT", "source": "amplicon",
          "copies": [ { "region": "1", "start": 31523, "end": 31627, "strand": 1, "variant_position": 31564,
                        "identity": 100, "gap_compressed_identity": 100, "aligned_length": 192, "observed": "TCT",
                        "flank_edits": 0, "call": "ref", "anchors": 2, "ortholog": null, "source": "amplicon" } ],
          "orthologous_copies": 1, "paralog_copies": 0, "reason": null },
        { "system_name": "sorghum_pi180348", "display_name": "Sb bicolor PI180348 Juar (IS 12876)", "is_reference": false,
          "allele": "alt", "observed": "TAT", "source": "amplicon",
          "copies": [
            { "region": "1", "start": 15028, "end": 15132, "strand": -1, "variant_position": 15091,
              "identity": 98.44, "gap_compressed_identity": 98.44, "aligned_length": 192, "observed": "TAT",
              "flank_edits": 2, "call": "alt", "anchors": 2, "ortholog": null, "source": "amplicon" },
            { "region": "1", "start": 38342, "end": 38446, "strand": 1, "variant_position": 38383,
              "identity": 98.45, "gap_compressed_identity": 98.45, "aligned_length": 193, "observed": "TAT",
              "flank_edits": 2, "call": "alt", "anchors": 2, "ortholog": null, "source": "amplicon" } ],
          "orthologous_copies": 2, "paralog_copies": 0, "reason": null },
        { "system_name": "sorghum_pi329250", "display_name": "Sb verticilliflorum PI329250 ADAR", "is_reference": false,
          "allele": "ref", "observed": "TCT", "source": "amplicon",
          "copies": [ { "region": "1", "start": 203934, "end": 204038, "strand": 1, "variant_position": 203975,
                        "identity": 100, "gap_compressed_identity": 100, "aligned_length": 192, "observed": "TCT",
                        "flank_edits": 0, "call": "ref", "anchors": 2, "ortholog": null, "source": "amplicon" } ],
          "orthologous_copies": 1, "paralog_copies": 0, "reason": null }
      ],
      "sets": [
        { "id": "S1", "ref_pair": "S1_REF", "alt_pair": "S1_ALT", "orientation": "reverse",
          "deliberate_mismatch_positions": [],
          "specificity": { "ref_pair": { "verdict": "specific", "consistent_with_allele": true },
                           "alt_pair": { "verdict": "specific", "consistent_with_allele": true }, "off_target_count": 0 },
          "control": { "status": "pass", "allele": "ref", "reasons": [] },
          "reference": { "system_name": "sorghum_bicolor",
                         "ref_primer": { "status": "match", "likelihood": "likely", "mm_pos": [], "residual_mm_pos": [] },
                         "alt_primer": { "status": "terminal_mismatch", "likelihood": "likely_weak", "mm_pos": [1], "residual_mm_pos": [1] },
                         "common_primer": { "status": "match", "likelihood": "likely", "mm_pos": [], "residual_mm_pos": [] },
                         "predicted": "ref", "strength": "normal", "agrees": true, "reasons": [], "off_locus_products": 0 },
          "summary": { "genomes_total": 3, "predicted_ref": 2, "predicted_alt": 1, "both": 0, "none": 0, "no_call": 0,
                       "unknown": 0, "weak": 0, "agree": 3, "disagree": 0, "not_comparable": 0 },
          "genomes": [
            { "system_name": "sorghum_bicolorv5",
              "ref_primer": { "status": "match", "likelihood": "likely", "mm_pos": [], "residual_mm_pos": [] },
              "alt_primer": { "status": "terminal_mismatch", "likelihood": "likely_weak", "mm_pos": [1], "residual_mm_pos": [1] },
              "common_primer": { "status": "match", "likelihood": "likely", "mm_pos": [], "residual_mm_pos": [] },
              "predicted": "ref", "strength": "normal", "agrees": true, "reasons": [], "off_locus_products": 0 },
            { "system_name": "sorghum_pi180348",
              "ref_primer": { "status": "terminal_mismatch", "likelihood": "likely_weak", "mm_pos": [1], "residual_mm_pos": [1] },
              "alt_primer": { "status": "match", "likelihood": "likely", "mm_pos": [], "residual_mm_pos": [] },
              "common_primer": { "status": "match", "likelihood": "likely", "mm_pos": [], "residual_mm_pos": [] },
              "predicted": "alt", "strength": "normal", "agrees": true, "reasons": [], "off_locus_products": 0 },
            { "system_name": "sorghum_pi329250",
              "ref_primer": { "status": "match", "likelihood": "likely", "mm_pos": [], "residual_mm_pos": [] },
              "alt_primer": { "status": "terminal_mismatch", "likelihood": "likely_weak", "mm_pos": [1], "residual_mm_pos": [1] },
              "common_primer": { "status": "match", "likelihood": "likely", "mm_pos": [], "residual_mm_pos": [] },
              "predicted": "ref", "strength": "normal", "agrees": true, "reasons": [], "off_locus_products": 0 }
          ] },
        { "id": "S2", "ref_pair": "S2_REF", "alt_pair": "S2_ALT", "orientation": "forward",
          "deliberate_mismatch_positions": [],
          "specificity": { "ref_pair": { "verdict": "specific", "consistent_with_allele": true },
                           "alt_pair": { "verdict": "specific", "consistent_with_allele": true }, "off_target_count": 0 },
          "control": { "status": "pass", "allele": "ref", "reasons": [] },
          "reference": { "system_name": "sorghum_bicolor",
                         "ref_primer": { "status": "match", "likelihood": "likely", "mm_pos": [], "residual_mm_pos": [] },
                         "alt_primer": { "status": "terminal_mismatch", "likelihood": "likely_weak", "mm_pos": [1], "residual_mm_pos": [1] },
                         "common_primer": { "status": "match", "likelihood": "likely", "mm_pos": [], "residual_mm_pos": [] },
                         "predicted": "ref", "strength": "normal", "agrees": true, "reasons": [], "off_locus_products": 0 },
          "summary": { "genomes_total": 3, "predicted_ref": 2, "predicted_alt": 1, "both": 0, "none": 0, "no_call": 0,
                       "unknown": 0, "weak": 0, "agree": 3, "disagree": 0, "not_comparable": 0 },
          "genomes": [
            { "system_name": "sorghum_bicolorv5",
              "ref_primer": { "status": "match", "likelihood": "likely", "mm_pos": [], "residual_mm_pos": [] },
              "alt_primer": { "status": "terminal_mismatch", "likelihood": "likely_weak", "mm_pos": [1], "residual_mm_pos": [1] },
              "common_primer": { "status": "match", "likelihood": "likely", "mm_pos": [], "residual_mm_pos": [] },
              "predicted": "ref", "strength": "normal", "agrees": true, "reasons": [], "off_locus_products": 0 },
            { "system_name": "sorghum_pi180348",
              "ref_primer": { "status": "terminal_mismatch", "likelihood": "likely_weak", "mm_pos": [1], "residual_mm_pos": [1] },
              "alt_primer": { "status": "match", "likelihood": "likely", "mm_pos": [], "residual_mm_pos": [] },
              "common_primer": { "status": "match", "likelihood": "likely", "mm_pos": [], "residual_mm_pos": [] },
              "predicted": "alt", "strength": "normal", "agrees": true, "reasons": [], "off_locus_products": 0 },
            { "system_name": "sorghum_pi329250",
              "ref_primer": { "status": "match", "likelihood": "likely", "mm_pos": [], "residual_mm_pos": [] },
              "alt_primer": { "status": "terminal_mismatch", "likelihood": "likely_weak", "mm_pos": [1], "residual_mm_pos": [1] },
              "common_primer": { "status": "match", "likelihood": "likely", "mm_pos": [], "residual_mm_pos": [] },
              "predicted": "ref", "strength": "normal", "agrees": true, "reasons": [], "off_locus_products": 0 }
          ] }
      ]
    }
  }
}
```

The `pangenome` block of the same job reports `amplifies: 3` for **both** `S1_REF` and `S1_ALT`, while `genotyping` shows the REF primer matching in only 2 of the 3 genomes. That contrast is exactly why consumers must read `results.genotyping` for allele questions (§5.5). pi329250 carries C although Ensembl genotypes that accession `A|A` [R]; the caller reports what the assembly contains.

The stored `job.request.genotyping` is byte-identical to what the client sent, except that alleles are uppercased and the variant is left-aligned; no server-derived key (no `algorithm`, no `orientation`) is added, so `job.request` still validates as a `PrimerCheckRequest` (§5.3).
### 2.14 Error codes

New handler error codes, grouped by HTTP status. No new status class is introduced: the documented set stays 400, 404, 405, 413, 422, 500, 503, 504 (`docs/primer_design_api.md:78-87`), which is also what `gramene-primers` `errors.ts` and `ErrorBanner` already handle.

| Status | Code | Endpoints | When | `details` |
|---|---|---|---|---|
| 400 | `VARIANT_WINDOW_TOO_LONG` | variants list | `end − start + 1 > variation.max_window` (50,000) | `{length, max}` |
| 400 | `INVALID_VARIANT` | genotyping design | wrong combination of variant fields, or invalid alleles (§2.7 rules 2–3) | `{reason: "id_or_manual" \| "alleles" \| "allele_too_long"}` |
| 400 | `REF_MISMATCH` | genotyping design, check | the given reference allele is not the genome's bases | `{region, position, given, genome}` |
| 400 | `ALT_REQUIRED` | genotyping design | an id whose site has more than one designable alt, and no `alt` | `{id, alts}` |
| 400 | `ALT_NOT_AT_SITE` | genotyping design | `alt` is not an allele of the id's site | `{id, alt, alleles}` |
| 400 | `UNSUPPORTED_ALLELE` | genotyping design | allele `*`, or containing `N` | `{allele}` |
| 400 | `VARIANT_TOO_CLOSE_TO_END` | genotyping design | no orientation has room for the smallest product | `{region, position, region_length, needed}` |
| 400 | `VARIANT_TOO_REPETITIVE` | genotyping design, check | the event can slide more than `variation.max_shift` (1,000 bp), so no stable discriminating base exists | `{region, position, shift, max}` |
| 400 | `GENOTYPING_SET_INVALID` | check | any structural or positional problem with a set (§5.2) | `{set_id, reason, pair_ids, primer?, sequence?, mm_pos?}` |
| 400 | `INVALID_REQUEST` (existing code, new cases) | check | `genotyping` in transcript or sequence mode | `{field: "genotyping", reason: "mode"}` |
| 404 | `UNKNOWN_VARIANT` | variant lookup, genotyping design | Ensembl answered 400 with a JSON body whose `error` matches `/not found/i` | `{id, system_name}` |
| 422 | `NO_VARIATION_DATA` | variants list and lookup, genotyping design by id | the genome has no variation data | `{system_name}` |
| 422 | `AMBIGUOUS_VARIANT_MAPPING` | variant lookup, genotyping design | the id maps to more than one sequence of the assembly | `{id, mappings}` |
| 422 | `VARIANT_NOT_ON_ASSEMBLY` | variant lookup, genotyping design | no mapping onto a sequence of this assembly | `{id}` |
| 500 | `THERMO_FAILED` | genotyping design | unparseable `ntthal` output, or a non-zero exit | `{binary}` |
| 503 | `THERMO_UNAVAILABLE` | genotyping design | the `ntthal` binary is missing or not executable (there is **no** startup validation of binaries today [V]; the failure appears on the first spawn) | `{binary, retry_after_s}` |
| 503 | `VARIATION_SOURCE_UNAVAILABLE` | variants endpoints, genotyping design by id | any upstream failure | `{retry_after_s, reason}` with `reason` ∈ `timeout`, `transport`, `http_5xx`, `rate_limited`, `invalid_response`, `breaker_open`, `queue_full` |

`GENOTYPING_SET_INVALID.reason` values: `unknown_pair`, `same_pair`, `pair_reused`, `expected_required`, `expected_differs`, `no_shared_common`, `not_at_variant`, `too_many_edits`, `alleles_swapped`, `common_in_zone`, `expected_mismatch`.

**Reused codes on the new endpoints:**

| Status | Codes |
|---|---|
| 400 | validator errors, `INVALID_REQUEST`, `INVALID_PARAMS`, `REGION_OUT_OF_BOUNDS`, `PRIMER3_INPUT_ERROR` |
| 404 | `UNKNOWN_GENOME`, `UNKNOWN_REGION` |
| 413 | body over 100 kb |
| 422 | `NO_SEQUENCE`, `AMBIGUOUS_ASSEMBLY`, `NO_BLASTDB`, `JOB_TOO_LARGE` |
| 500 | `PRIMER3_FAILED`, `INTERNAL` |
| 503 | `BUSY`, `PRIMER3_UNAVAILABLE`, `FEATURE_DISABLED`, `MONGO_UNAVAILABLE` |
| 504 | `DEADLINE_EXCEEDED` |

`errors.js normalizeStatus` (`errors.js:34-37`) accepts any status in 400–599, so no change is needed there.

### 2.15 Warning codes

Warnings keep `{code, message}` and add `details` (now declared in swagger). One name per condition, shared by all endpoints.

**Variants endpoints**

| Code | Details |
|---|---|
| `DUPLICATE_VARIANT_IDS` | `{key, ids}` |
| `VARIATION_RECORDS_SKIPPED` | `{count, reasons}` — malformed Ensembl records dropped; ids are **never** a reason (§3.4) |
| `REF_MISMATCHES` | `{count}` — entries whose REF disagrees with the local FASTA |
| `VARIANTS_TRUNCATED` | `{returned, total, limit}` |

**Genotyping design, response level**

| Code | When | Details |
|---|---|---|
| `NO_SETS` | no set in any orientation | `{}` |
| `ORIENTATION_BLOCKED` | a non-EMS neighbour inside an AS primer's 3′ window (natural targets) | `{orientation, ids, distances}` |
| `ORIENTATION_NO_SETS` | the ladder was exhausted | `{orientation, levels_tried}` |
| `ORIENTATION_SKIPPED` | too close to the region end, or `N` bases in the AS window | `{orientation, reason}` |
| `RELAXED_CONSTRAINTS` | sets needed level ≥ 1 | `{orientation, level, changes}` |
| `DUPLICATE_VARIANT_IDS` | merged ids | `{key, ids}` |
| `MULTIALLELIC_SITE` | the site has other alternative alleles | `{key, other_alts}` |
| `SHIFTABLE_INDEL` | `shift > 0` | `{shift, forward_position, reverse_position, zone}` |
| `DENSE_NEIGHBOURS` | more than `dense_count` (2) non-EMS neighbour keys whose `minimal` span overlaps `[vcf.position − 30, vcf.position + len(vcf.ref) − 1 + 30]` | `{count, window}` |
| `EMS_TARGET` | every record of the target is EMS; natural neighbours warn instead of blocking | `{sources}` |
| `NO_VARIATION_DATA` | manual variant on a genome without variation data | `{system_name}` |
| `NEIGHBOURS_UNAVAILABLE` | Ensembl unavailable during a manual design | `{reason}` |
| `VARIANT_IN_REPEAT` | the repeat mask overlapped the exempted AS window | `{orientation, masked_bases}` |
| `SUBMISSION_NEIGHBOURS_OMITTED` | indel or multi-allelic neighbours in the 50 bp submission flanks | `{ids}` |
| `DESIGN_BUDGET_EXHAUSTED` | the next candidate's reserved Primer3 runs or thermo calls would exceed `max_primer3_runs` or `max_thermo_calls` (§4.8); the remaining candidates are left unscored and the response is a normal `200` with the sets scored so far, never a `500` | `{orientation, primer3_runs, thermo_calls, max_primer3_runs, max_thermo_calls, not_scored, sets_returned}`, counters as at the end of the request |
| `ASSEMBLY_MISMATCH`, `AMBIGUOUS_ASSEMBLY`, `REPEAT_MASK_FAILED`, `NO_REPEAT_MASK`, `BLAST_DEPTH_MASK`, `MOSTLY_REPEAT`, `PRIMER3_WARNING` | existing meanings | as today |

**Genotyping design, set level** (also emitted as `sets[].issues[]` with a severity)

| Code | When | Severity | Details |
|---|---|---|---|
| `AS_TM_IMBALANCE` | `as_tm_diff > 1.0 °C` | warn (high above 2.0) | `{diff}` |
| `COMMON_TM_OUT_OF_RANGE` | `common_minus_as` outside −1.0 … +3.0 °C | warn | `{value}` |
| `TAILED_STRUCTURE` | a tailed hairpin ≥ 47 °C, or a tailed dimer (`ANY`/`END1`) ≥ 47 °C | warn; **high** when a hairpin ≥ 55 °C | `{oligo, metric, value, severity}` |
| `MISMATCH_STRUCTURE` | a deliberate-mismatch primer's hairpin or self value ≥ 47 °C | warn; high ≥ 55 °C | `{oligo, metric, value}` |
| `MISMATCH_NOT_APPLICABLE` | the base at −k differs between the haplotypes, so one mismatch would change the discriminating chemistry | warn | `{position}` |
| `NEIGHBOUR_IN_PRIMER` | non-EMS neighbours under a primer, outside its 3′ window; at most one for the allele-specific pair (`role: "as_ref"`, the union of both allele-specific primers' neighbours) and one for the common primer (§4.15) | warn | `{role, ids, distances}` |
| `NEIGHBOUR_AT_3P` | a non-EMS neighbour inside a primer's 3′ window that was *not* allowed to block (EMS target, or `neighbour_policy: ignore`); grouped like `NEIGHBOUR_IN_PRIMER` | high | `{role, ids, distances}` |
| `WEAK_DISCRIMINATION` | the other-allele likelihood of an AS primer is `likely` | high | `{primer, mm_pos}` |
| `SHIFT_TRACT_DISCRIMINATION` | an allele-specific primer has `in_shift_tract: true`: its 3′ base, on its own haplotype, lies inside the variant's shift tract (§3.5, §4.12); one issue per such primer | warn | `{primer, shift}` |
| `ALT_PRIMER_SUBOPTIMAL` | an ordered derived oligo echoes a non-empty `primer3_problems`, yet clears the hard floors: `as_alt` without a deliberate mismatch (from its ALT-haplotype run), or either allele-specific primer with one (from its own mismatch run). At most one per oligo; with a deliberate mismatch the matched ALT primer's own ALT run never raises it (§4.16 S6) | warn | `{oligo, problems}` |
| `WEAK_TERMINAL_CLASS` | both allele-specific primers have `terminal_mismatch_class` `weak` (A/G and C/T sites) | info (issues only; not in `warnings`, not scored) | `{}` |

**Check, results level** (`results.warnings`)

| Code | When |
|---|---|
| `REFERENCE_CONTROL_FAILED` | the reference genome is not called `ref`, or a set does not predict `ref` on it |
| `GENOTYPE_FALLBACK_FAILED` | a megablast fallback failed; the message lists the genomes |
| `GENOTYPE_FALLBACK_BUDGET` | the per-job fallback budget (30) was reached; later genomes are `missing` with `reason: "fallback_budget"` |
| `GENOTYPE_FAILED` | the caller threw for a genome; that genome is `unavailable` with `reason: "call_failed"` |
---

## 3. Variant resolution and normalization

### 3.1 Configuration

Added to `config.js` `DEFAULTS` (`config.js:19-99`) and to the `primers:` block of `config/default.yaml`:

```yaml
primers:
  ntthal: /home/olson/primer3-2.6.1/bin/ntthal            # env NTTHAL
  variation:
    enabled: true                                          # env PRIMERS_VARIATION_ENABLED
    base_url: https://data.gramene.org/pansite-ensembl-115 # env PRIMERS_VARIATION_URL
    release: '115'
    species: { sorghum_bicolor: sorghum_bicolor }          # system_name -> Ensembl species segment; nothing else is ever called
    timeout_ms: 8000
    chunk_bp: 10000                # overlap windows are fetched in fixed 10 kb chunks
    max_chunk_bytes: 1000000       # a 10 kb chunk measured 186 KB; a 50 kb window 1.0 MB
    max_lookup_bytes: 262144
    max_window: 50000
    list_limit_default: 2000
    list_limit_max: 5000
    max_concurrent: 4              # process-wide outbound limit
    queue_wait_ms: 5000
    cache_entries: 500
    cache_ttl_ms: 3600000          # release 115 is static
    unknown_cache_ttl_ms: 300000
    unavailable_cache_ttl_ms: 30000
    breaker_failures: 3
    breaker_window_ms: 60000
    breaker_open_ms: 60000
    retry_after_s: 30
    queue_retry_after_s: 5
    max_allele_length: 50
    max_shift: 1000
    ems_source_pattern: '^EMS_'
  genotyping:
    template_flank: 400
    num_return_per_run: 20
    num_sets_default: 6
    max_sets: 10
    max_scored_per_orientation: 8  # candidates scored per orientation, across all levels (§4.8)
    max_primer3_runs: 54           # 2 x 3 design runs + 2 x 8 x 3 scoring runs (worst case: deliberate mismatch)
    max_thermo_calls: 272          # 2 x 8 x 17 ntthal calls (worst case: tails plus deliberate mismatch)
    thermo_concurrency: 4
    thermo_timeout_ms: 5000
    guard_gap: 10                  # SEQUENCE_TARGET length just outside the zone
    mask_exempt_pad: 36            # Primer3's built-in maximum primer length
    as_min_tm: 52                  # hard floor, every allele-specific primer
    as_min_gc: 15                  # hard floor, derived ALT and mismatch primers only
    structure_warn_th: 47
    structure_high_th: 55
    as_tm_diff_warn: 1.0
    common_tm_low: -1.0
    common_tm_high: 3.0
    neighbour_3p_window: 5
    dense_window: 30
    dense_count: 2
    check_max_sets: 5
    check_max_unique_primers: 13
  check:
    genotype_cpu_s_per_genome: 0.2
    genotype_flank_min: 15
    genotype_amplicon_pad: 50
    genotype_ortholog_min_identity: 95
    genotype_ortholog_size_tolerance: 0.2
    genotype_max_copies: 10
    genotype_max_anchors: 50
    genotype_max_megablast: 30
    genotype_megablast_timeout_ms: 20000
    genotype_megablast_min_identity: 95
    genotype_megablast_min_query_cover: 0.8
    genotype_megablast_min_bitscore_frac: 0.9
```

Env overrides go in `ENV_OVERRIDES` (`config.js:103-115`) with `applyEnv` (`:163-185`): `NTTHAL` (path), `PRIMERS_VARIATION_URL` (new `url` type: accepted only when it starts with `https://`, or with `http://127.0.0.1` for the loopback fake server used in tests; otherwise ignored with a warning), `PRIMERS_VARIATION_ENABLED` (bool).

There is **no startup validation of binaries today** — neither `config.js`, `app.js` nor the controllers check `primer3_core` [V] — so `ntthal` is likewise checked only when it is first spawned, and a missing binary becomes `503 THERMO_UNAVAILABLE` (mirroring `primer3.js:30-42`).

### 3.2 Which genomes have variation data

`hasVariation(system_name)` is true when `variation.enabled`, `variation.species` has the key, and the assembly resolves with a FASTA. Only `sorghum_bicolor` qualifies, and its Ensembl coordinates match the local FASTA: 1:11109 = C, 1:11282–11286 = CAAAG, 1:11500–11503 = GCCA, 1:11193 = C [V]. The server never calls Ensembl for any other genome; the REST service answers other species with HTTP 400 and a Perl stack trace [R]. `GET /primers/genomes` exposes the flag so the UI never guesses.

### 3.3 The Ensembl client (`variation/client.js`)

| Purpose | Request |
|---|---|
| Window records | `GET {base_url}/overlap/region/{species}/{region}:{chunkStart}-{chunkEnd}?feature=variation;content-type=application/json` |
| One id | `GET {base_url}/variation/{species}/{id}?content-type=application/json` |

**URL building (security).** `species` comes from config only. `region` must first be a sequence of the resolved assembly (`sequence.regionLength` returns a number) and is then `encodeURIComponent`-ed. Coordinates are safe integers. The id must match `^[A-Za-z0-9][A-Za-z0-9_.:,*-]{0,254}$` (swagger, and re-checked in the handler) and is passed through `encodeURIComponent`; Ensembl answers 200 for `tmp_1_13549_TTA_T,*` both literally and percent-encoded [R]. The URL is built with `new URL(path, base + '/')` and asserted to start with `base_url`. A leading alphanumeric is required, which rules out `.` and `..` path segments.

**Transport.** Node 24 global `fetch`, `redirect: 'error'`, `Accept: application/json`, `User-Agent: gramene-swagger-primers`, `signal: AbortSignal.timeout(timeout_ms)` and nothing else: the shared request never carries a caller's deadline. Each caller awaits `Promise.race([shared, rejectOnAbort(callerDeadline.signal)])`, so a caller whose deadline fires stops waiting while the shared fetch, and the cache entry it fills, continue for everyone else. The body is read through `res.body.getReader()` counting bytes and aborted above the cap (`max_chunk_bytes` for overlaps, `max_lookup_bytes` for lookups).

**Resilience — all four mechanisms are required (§13, finding 6):**
1. **Process-wide limiter.** At most `max_concurrent` (4) outbound calls at a time; a waiter that does not start within `queue_wait_ms` (5 s) fails with `503 VARIATION_SOURCE_UNAVAILABLE {reason: "queue_full", retry_after_s: 5}`.
2. **Single flight.** One in-flight promise per cache key; concurrent callers share it. A cancelled caller never cancels the shared request.
3. **Circuit breaker.** `breaker_failures` (3) transport failures inside `breaker_window_ms` (60 s) open the breaker for `breaker_open_ms` (60 s); while open, calls fail immediately with `reason: "breaker_open"` and `retry_after_s` = the remaining seconds. One success closes it.
4. **Bounded caches.** One LRU of `cache_entries` (500) holding **validated compact records only**, never raw bodies. Keys: `o|{species}|{region}|{chunkStart}` and `v|{species}|{id}`. Positive TTL 1 h (release 115 is static); `UNKNOWN_VARIANT` is cached for 5 min; transport failures are cached negatively for 30 s.

**Chunking.** Overlap windows are always fetched in fixed `chunk_bp` (10 kb) chunks covering `[start − 1, end + 1]`, and filtered locally. A picker listing, a later design over the same locus and the UI's neighbour-track request therefore share the same cache entries. A 50 kb window needs at most 6 chunks, fetched concurrently under the limiter.

**Status mapping.**

| Ensembl answer | Outcome |
|---|---|
| 200 with valid JSON of the expected shape | parsed records, cached |
| 400 on `/variation/…` whose JSON body has `error` matching `/not found/i` [R] | `404 UNKNOWN_VARIANT`, cached 5 min, **not** a breaker failure |
| any other 4xx, a non-JSON body (e.g. an Apache HTML 404 from a mistyped base URL [R]), a body over the cap, or an unexpected shape | `503 VARIATION_SOURCE_UNAVAILABLE {reason: "invalid_response"}`, not cached, counts toward the breaker, logged with the first 200 characters and the host only |
| 429, 5xx, timeout, network error | `503 … {reason: "rate_limited" \| "http_5xx" \| "timeout" \| "transport", retry_after_s: 30}`, negatively cached, counts toward the breaker |

A misconfigured base URL therefore never masquerades as "unknown variant".

**Record validation.** Overlap records are kept when `feature_type === "variation"`, `seq_region_name` equals the requested region, `start`/`end` are safe integers, and `alleles` is an array of 2–10 strings each matching `^([ACGTNacgtn]+|-|\*)$`. **Ids are validated leniently** — any string of at most 255 characters with no control characters — because a record id never enters a URL; 88 of 3,976 real ids in one 50 kb window contain `,` or `*` and 11 exceed 128 characters [V], and dropping them would silently remove real neighbours from screening. `source` and `consequence_type` are truncated to 100 characters. Dropped records are counted in `VARIATION_RECORDS_SKIPPED {count, reasons}`.

For `/variation` lookups the body must be an object with `mappings[]`, each with `seq_region_name`, `start`, `end` and `allele_string`. Mappings onto sequences that are not in this assembly are ignored; zero usable mappings → `422 VARIANT_NOT_ON_ASSEMBLY`, more than one → `422 AMBIGUOUS_VARIANT_MAPPING`. `synonyms` keeps only strings matching the id pattern, which drops Ensembl's `"."` [V]. The batch `POST /variation` endpoint silently omits unknown ids [R] and is not used.

### 3.4 From records to canonical entries (`variation/normalize.js`, pure)

```
for each record r (alleles[0] = REF as reported, alleles[1..] = alts):
  refOk = alleles[0] == '-' ? true : fasta(start..end) == upper(alleles[0])
  for each alt a in alleles[1..]:
    if a == '*': add issue STAR_ALLELE to the sibling entries; continue
    if alleles[0] == '-':          # insertion: BOTH conventions mean "insert before start"
        p = start - 1              #   EVA: end == start      (rs5413863549 1:11503-11503)
        vcf = (p, base(p), base(p) + a)   #   SAP/BAP: end == start - 1 (tmp_1_11502_C_CGT 1:11503-11502)
    elif a == '-':                 # deletion of start..end
        p = start - 1
        vcf = (p, base(p) + alleles[0], base(p))
    else:
        vcf = (start, alleles[0], a)
    v = leftNormalize(vcf)
    key = region:v.pos:v.ref:v.alt
    merge into entry[key]: ids += r.id, records += {id, source, ems}
```

| Record | Convention | `key` | `shift` | forward / reverse discriminating |
|---|---|---|---|---|
| `rs871475760` 1:11109-11109 `C/A` | SNV | `1:11109:C:A` | 0 | 11109 C→A / 11109 C→A |
| `tmp_1_11193_C_T` 1:11193-11193 `C/T` (EMS) | SNV | `1:11193:C:T` | 0 | 11193 / 11193 |
| `rs5413864115` 1:11283-11283 `A/-` (synonym `tmp_1_11282_CA_C`) | deletion inside AAA | `1:11282:CA:C` | 2 | 11285 A vs G / 11283 A vs C |
| `tmp_1_11502_C_CGT` 1:11503-11502 `-/GT` (SAP) | insertion, start = end + 1 | `1:11502:C:CGT` | 0 | 11503 A vs inserted G / 11502 C vs inserted T |
| `rs5413863549` 1:11503-11503 `-/GT` (EVA) | insertion, start = end | `1:11502:C:CGT`, **merged** | 0 | same |
| `rs5413863413` 1:13736-13738 `TGG/-` | deletion in five TGG repeats | `1:13735:TTGG:T` | 13 | 13749 G vs A / 13738 G vs T |

All six rows are verified against the genome and live Ensembl [V].

**Merged entries.** `ids` lists the requested id first, then `rs` ids, then the rest alphabetically. `ems` is true only when every merged record is EMS. `consequence` is the first non-empty `consequence_type`. `ref_verified` is false when any merged record failed `refOk`; such an entry gets issue `REF_MISMATCH` and `designable: false`.

**Multi-allelic sites and `*`.** A record with more than one non-`*` alt yields one entry per alt, each carrying `multiallelic {alleles, other_alts}` [V: `rs5413863494` C/T/G → `1:10718:C:T` and `1:10718:C:G`]. A `*` allele is never an entry; it adds issue `STAR_ALLELE` to its siblings, which stay designable [V: `rs5413863901` C/T/`*`]. Designing a multi-allelic entry adds warning `MULTIALLELIC_SITE`: assemblies carrying a third allele mismatch both AS primers at the 3′ base, so the check calls them `other` and predicts `none`. An allele containing `N` gives `UNSUPPORTED_ALLELE` and `designable: false`; an allele longer than `max_allele_length` (50) gives `ALLELE_TOO_LONG`. Equal-length multi-base substitutions are `mnv`, unequal lengths without a shared anchor are `complex`; both are designable.

### 3.5 Left normalization, shift, discriminating positions and the zone

```
leftNormalize(pos, ref, alt):                       # vt-style
  loop:
    changed = false
    if ref and alt non-empty and last(ref) == last(alt): drop last base of both; changed = true
    if ref or alt is empty: pos -= 1; prepend base(pos) to both; changed = true
    if not changed: break
  while len(ref) >= 2 and len(alt) >= 2 and first(ref) == first(alt): drop first base of both; pos += 1

shift(v):   number of times the event can slide one base right and describe the same haplotype
            (shift > max_shift 1,000 -> 400 VARIANT_TOO_REPETITIVE)

tract(v):   shift(v) > 0 ? [v.pos + 1, v.pos + len(v.ref) - 1 + shift(v)] : null
            # the shift tract: the reference bases an insertion or deletion can slide through, anchor base v.pos excluded;
            # only an anchored indel can slide, so an SNV, MNV, complex event or shift-0 indel has no tract

discriminating(v):
  W = 200                                            # doubled up to 1,600 if the haplotypes do not differ inside W
  a = v.pos - W ; b = v.pos + len(v.ref) - 1 + W     # clamped to the region
  H_ref = fasta(a..b)
  H_alt = fasta(a..v.pos-1) + v.alt + fasta(v.pos+len(v.ref)..b)
  i = first index where H_ref[i] != H_alt[i]                     # left to right  -> forward
  j = first index from the end where they differ                 # right to left  -> reverse
  forward = {position: a+i,              ref_base: H_ref[i],      alt_base: H_alt[i],      alt_maps_to: altToRef(i)}
  reverse = {position: a+len(H_ref)-1-j, ref_base: H_ref[-1-j],   alt_base: H_alt[-1-j],   alt_maps_to: altToRef(len(H_alt)-1-j)}

altToRef(k):   index before the event      -> the same coordinate
               inside v.alt, base shared with v.ref (the anchor, or an MNV base) -> v.pos + offset
               otherwise inside v.alt      -> null (an inserted base)
               after the event             -> v.pos + len(v.ref) + (offset - len(v.alt))

zone(v) = [ min(v.pos, reverse.position, reverse.alt_maps_to ?? +inf),
            max(v.pos + len(v.ref) - 1, forward.position, forward.alt_maps_to ?? -inf) ]
```

The **zone** is the single most important derived value: it is the span that no primer other than an allele-specific primer may touch, on either haplotype. Verified zones [V]: `1:11109:C:A` → [11109, 11109]; `1:11193:C:T` → [11193, 11193]; `1:11282:CA:C` → [11282, 11286]; `1:11502:C:CGT` → [11502, 11503]; `1:13735:TTGG:T` → [13735, 13752]. The last one is the case a REF-only guard misses: a common primer starting at 13751 satisfies "after the REF deletion" yet covers the ALT primer's 3′ base on ALT carriers.

**The shift tract** [V]. `tract(v)` is used in one place only: `in_shift_tract` (§4.12), which raises `SHIFT_TRACT_DISCRIMINATION` and makes the check report a primer `uncertain` (§5.7). It is not the §5.6 core, which is a separate span. For an indel with a tract the zone is always the tract plus one base on each side, and the geometry is fixed:
- a **deletion**'s REF discriminating bases lie inside the tract (forward at `tract.end − len(v.ref) + 2`, reverse at `v.pos + len(v.ref) − 1`) and its ALT anchor mappings just outside it (`tract.end + 1` forward, the anchor `v.pos` reverse);
- an **insertion**'s REF discriminating bases lie just outside the tract (`tract.end + 1` forward, `v.pos` reverse) and its ALT 3′ bases inside it or on inserted bases.

Values: `1:11282:CA:C` has tract [11283, 11285], the AAA run (zone [11282, 11286]); `1:13735:TTGG:T` has [13736, 13751] (zone [13735, 13752]); the other three example variants have none. `scratch-fix2/tract56.py` checks all of these relations for the 432 shiftable indels (235 deletions, 197 insertions) among the 3,978 Ensembl records of 1:11080–61079.

**Why the anchor base is excluded.** The anchor is not a repeat base. A primer ending on it (the reverse ALT primer of every deletion), or on the first base after the tract (the forward ALT primer), meets the other allele with a plain terminal mismatch, as at any non-repeat indel, and 1-bp indels genotype reliably [R]. If the tract included the anchor, the reverse ALT primer of all 235 shiftable deletions would be `uncertain`, so every REF genome would predict `unknown` for a reverse set, and rs5413864115 S2 would gain a second `SHIFT_TRACT_DISCRIMINATION` and score 24.95.

**Why the ALT primer always keeps the REF primer's length.** Forward: the indices before `i` are identical on both haplotypes, so a primer from the same 5′ start to index `i` has the same length. Reverse: the suffixes after `j` are identical, so both 3′ and 5′ ends shift by the same `delta`.

### 3.6 REF verification and FASTA access

One `sequence.fetch` per operation, never one per event (1,000 small fetches cost 1.15 s of CPU and stalled the event loop for up to 685 ms in measurement [R]; one 53.2 kb fetch costs 2.7 ms):

| Operation | Window read once |
|---|---|
| `GET /primers/variants` | `[start − (max_shift + 250), end + (max_shift + 250)]`, clamped — ≤ 53 kb, far below `design.max_fetch_length` |
| `GET /primers/variants/{id}` | `[pos − (max_shift + 250), pos + len(ref) + max_shift + 250]` |
| genotyping design | the template window ± `max_shift` |
| check submit | `[position − 1,700, position + len(ref) + 1,700]` (§5.2) |

Normalization then runs as pure functions over that string, re-fetching only if a shift reaches the window edge. A manual REF mismatch is `400 REF_MISMATCH {region, position, given, genome}`; an entry resolved by id whose `ref_verified` is false is the same error.

### 3.7 EMS flag

`ems` comes from the overlap endpoint's short `source` code matched against `ems_source_pattern` (`^EMS_`). Id lookups take `source` from the overlap records at the same site, because `/variation` returns a long title instead [R]. In listings the flag is informational; in design it changes only the neighbour policy (§4.15) and sets `assay.ems_target` when every record of the *target* is EMS.

### 3.8 When Ensembl is unavailable

| Call | Outcome | Consumer guidance |
|---|---|---|
| `GET /primers/variants`, `GET /primers/variants/{id}` | `503 VARIATION_SOURCE_UNAVAILABLE` with `retry_after_s` and `reason`; repeats inside the negative-cache or breaker window answer immediately | banner with countdown, offer manual entry |
| Design with `variant.id` | the same 503 | fall back to manual entry from the last known `vcf` |
| Design with a manual variant on a genome **with** variation data | `200` with `neighbours.data: "unavailable"` and warning `NEIGHBOURS_UNAVAILABLE` | sets are designed, nothing was screened |
| Design on a genome **without** variation data | `200` with `neighbours.data: "none"` and warning `NO_VARIATION_DATA` | as designed |
| `POST /primers/check`, and the worker | never calls Ensembl | — |
| Ordinary `POST /primers/design` | never calls Ensembl, and never waits behind one (§4.2) | — |

### 3.9 Manual variants

A manual input becomes a pseudo-record: VCF style `{start: position, end: position + len(ref) − 1, alleles: [ref, alt]}`; Ensembl style is passed to §3.4 unchanged with `id: null`. On a genome with variation data the design's template-window overlap also fills `variant.ids` and `records` for the same key, so a manual `1:11109 C/A` is reported as `rs871475760` [V].
---

## 4. Design algorithm (`POST /primers/genotyping/design`)

### 4.1 Pipeline

| Step | What | Module / function | Budget |
|---|---|---|---|
| 1 | Validate the request (§2.7 rules 1–6) and start the 45 s deadline | `genotyping/request.js normalize`, `design.createDeadline` (`design.js:552-584`) | pure |
| 2 | Resolve the assembly and region length | `assemblies.resolve`, `sequence.regionLength` | cached |
| 3 | Resolve the variant (id → Ensembl, or manual) and verify REF | `variation.resolveDesignVariant` | ≤ 2 outbound calls, 8 s each |
| 4 | Normalize; compute `shift`, both discriminating positions and the zone | `variation/normalize.js` | pure over one FASTA read |
| 5 | Fetch the neighbours of the template window (10 kb chunks, cached) | `variation.client.overlap` | ≤ 6 outbound calls, usually 0–1 |
| 6 | **Acquire the design semaphore** | `semaphore.designSemaphore().acquire({signal: dl.signal})` (`semaphore.js:73-100`) | 4 running / 16 waiting |
| 7 | Build the REF and ALT templates and `features` | `genotyping/template.js buildVariantTemplate` | 1 FASTA read |
| 8 | Optional repeat mask with the AS window exempt | `repeat_mask.repeatMask` (`repeat_mask.js:445-476`) | 1 megablast, ~0.15 s |
| 9 | Per orientation: blocker check, then ladder L0…L`max_relaxation`, one Primer3 run per level | `genotyping/design.js runOrientation` | ≤ 6 runs |
| 10 | Pairs → candidate sets (§4.8 filters; at most `max_scored_per_orientation` (8) scored per orientation, across levels) | `genotyping/sets.js candidatesFromPairs` | pure |
| 11 | Derive ALT primers; apply the deliberate mismatch; score each candidate with `check_primers`; apply the hard floors | `sets.js deriveAlt`, `mismatch.js`, `primer3.run` | ≤ 48 runs: 8 per orientation × 1, or × 3 with a deliberate mismatch; `max_primer3_runs` 54 with the design runs |
| 12 | Discrimination on both haplotypes with the check's own DP and classifier | `realign.realignSiteOnWindow` (`realign.js:196-221`), `classify.classifyAmplicon` (`classify.js:71-78`) | pure |
| 13 | Tail structures and tailed cross-dimers | `thermo.js` (`ntthal`) | ≤ 272 calls: 8 per orientation × 17 (15 tailed structures + 2 mismatch duplexes), memoized; pool of 4 |
| 14 | Issues, quality, score, ranking, cut to `num_sets` | `sets.js rankSets` | pure |
| 15 | Order rows, check fragments, the proposed check request | `genotyping/order.js`, `sets.js` | pure |
| 16 | Assemble the response; release the slot | `genotyping/design.js` | pure |

`template_only: true` stops after step 8 and returns `variant`, `template`, `neighbours`, `settings`, `engine`, `warnings`, with `sets: []`, `orientations: null`, `check: null`.

### 4.2 Deadline, semaphore and Ensembl (resolves finding 6)

- `createDeadline(design.deadline_ms)` starts at step 1, **before** the semaphore, and every child process and outbound call receives `min(its own timeout, deadline remaining)`. The end-to-end bound is 45 s, semaphore wait included.
- `semaphore.acquire` is given the **deadline signal**, not only the client signal (`design.js:671` passes only `deps.signal` today); if `acquire` rejects, the deadline is disposed before rethrowing.
- Steps 3–5 — every Ensembl call — run **before** `acquire`. A Primer3 slot is never held across network I/O. This is the difference that keeps an Ensembl brownout from starving ordinary `/primers/design` traffic: with 4 slots, 16 waiters and a 10 s wait, four genotyping designs that each waited 8 s inside a slot would hand `503 BUSY` to every PCR user.
- `template_only` requests skip the semaphore entirely **unless** `avoid_repeats` is true, because only the repeat mask spawns a child.
- Outbound calls also obey the process-wide limiter, single flight and breaker of §3.3.

### 4.3 Template window

- **Flank** `F = max(genotyping.template_flank (400), max product-range upper bound + 40)`.
- `template.start = max(1, min(fwd, rev) − F)`, `template.end = min(region_length, max(fwd, rev) + F)`. Verified windows [V]: rs871475760 1:10709–11509 (801 nt); tmp_1_11193_C_T 1:10793–11593 (801); rs5413864115 1:10883–11685 (803 REF / 802 ALT); tmp_1_11502_C_CGT 1:11102–11903 (802 REF / 804 ALT).
- **Room.** Forward needs `template.end − zone.end + 1 ≥ product_min`; reverse needs `zone.start − template.start + 1 ≥ product_min`. An orientation without room is `skipped {reason: "too_close_to_end"}`; if both lack room the request is `400 VARIANT_TOO_CLOSE_TO_END`.
- `template.seq` is the uppercase plus-strand REF sequence; `template.alt_seq` is the same with `vcf.ref` replaced by `vcf.alt`.
- An `N` inside `[t_fwd − max_size + 1, zone.end]` (forward) or `[zone.start, t_rev + max_size − 1]` (reverse) skips that orientation with `{reason: "n_in_primer_window"}`.
- `template.features` carries `variant`, `zone`, `discriminating`, `alt_offset` and `exempt`, all in template coordinates; intervals use the existing `PrimerInterval` `[start, length]` encoding (`swagger.yaml:2076-2083`) [V].

### 4.4 The variant zone and the common-primer guard (resolves finding 4)

The zone (§3.5) is the span that only an allele-specific primer may touch. It already includes both discriminating positions and both ALT 3′ anchor mappings, so it is correct for multi-base and shiftable indels and for MNVs.

**Primer3-native guard.** Each design run adds one `SEQUENCE_TARGET` of `guard_gap` (10) nt immediately outside the zone, on the side where the common primer must sit:

| Orientation | Forced tag | Guard |
|---|---|---|
| forward | `SEQUENCE_FORCE_LEFT_END = t(discriminating.forward.position)` | `SEQUENCE_TARGET = "{t(zone.end)+1},{gap}"` |
| reverse | `SEQUENCE_FORCE_RIGHT_END = t(discriminating.reverse.position)` | `SEQUENCE_TARGET = "{t(zone.start)−gap},{gap}"` |

Primer3 requires the pair to flank a target and rejects any primer overlapping it, so the common primer is pushed at least 10 nt beyond the zone while the forced allele-specific primer — which lies on the other side of the target — is untouched. The gap is shortened (never below 1) when the template edge is closer.

**Verified** [V]:
- Without the guard and with a 50 bp product minimum, Primer3 returns `PRIMER_RIGHT_0=423,24` = `TCTTTGACTAGCGAGAAATTCAGA`, footprint 1:11108–11131, which covers the variant and overlaps the allele-specific primer; `amplicons.finalize` (`amplicons.js:275-282`) would later discard that product and the set would come back `on_target_missing`.
- With `SEQUENCE_TARGET=402,10` the same record returns the clean 93 bp pair (`PRIMER_RIGHT_0=464,24`).
- The forced primer survives in both orientations and for a shiftable deletion: forward `375,29` + `TARGET=405,10` → common `462,22`, 88 bp; reverse `425,25` + `TARGET=390,10` → common `297,28`, 129 bp.

**Post-filter, on both haplotypes** (belt and braces, and it catches anything Primer3 might do differently):
1. the common primer's REF footprint must not intersect `[zone.start − 1, zone.end + 1]`;
2. mapped onto the ALT haplotype, it must not intersect either allele-specific primer's footprint;
3. the ALT-pair `check_primers` run must succeed — if the common primer does not exist verbatim on the ALT haplotype, Primer3 answers `PRIMER_ERROR=Specified right primer not in sequence` [V] and the set is dropped with `rejected.alt_scoring_failed`.

### 4.5 Repeat masking and the variant exemption

- `avoid_repeats` defaults to **false**, exactly as `/primers/design`. Measured consequence of turning it on for the rs871475760 window: the masker returns `mask [[1, 200]]`, `mask_source: "blast_depth"`, `masked_fraction 0.2497` and warning `BLAST_DEPTH_MASK` — a quarter of the window, including the region that supplies the verified gel AS-PCR common primer at 1:10880–10903 [V]. The default is listed for the user to confirm (§11), with the recommendation to switch it on once the masked-window behaviour has been reviewed.
- When it is on, the existing masker runs on a region-shaped template object and `template.mask`, `mask_source` and `masked_fraction` are reported as in design.
- **Exemption.** Each orientation gets its own copy of the Primer3 template in which the allele-specific window is unmasked: forward `[t_fwd − mask_exempt_pad, t(zone.end)]`, reverse `[t(zone.start), t_rev + mask_exempt_pad]`, with `mask_exempt_pad` = 36, Primer3's built-in maximum primer length, so the exemption covers every ladder level. The exempted spans are echoed in `template.features.exempt`.
- **Why.** A masked base inside the allele-specific footprint kills every forced candidate: the research measured `too many Ns 15, ok 0` with `n_mask` and `lowercase masking of 3' end 15, ok 0` with `three_prime` [R]. The common primer still avoids the mask.
- If any exempted base was masked, warning `VARIANT_IN_REPEAT {orientation, masked_bases}`.

### 4.6 Presets, the relaxation ladder and the hard floors (resolves finding 3)

Level 0 is the preset; each level applies its changes on top of the previous one, **except** params the client set, which are pinned and never relaxed. At every level the effective product-range minimum is raised to `max(range_min, 2 × max_size + 1)` so that a pair can never be shorter than its two primers.

**`kasp`**

| Param | L0 | L1 | L2 | Source |
|---|---|---|---|---|
| `opt_size` / `min_size` / `max_size` | 22 / 18 / 30 | max 32 | — | EasyKASP AS 21–27 nt, common 22–32 [R]; AT-rich sorghum loci needed 28–31 nt [V] |
| `opt_tm` / `min_tm` / `max_tm` | 60 / 57 / 63 | 55 / 65 | min 52 | EasyKASP AS 57.5 ± 1.5 °C, common 60.5 ± 1.5 [R]; LGC touchdown 61→55 °C then 26 cycles at 55 °C [R] |
| `min_gc` / `max_gc` | 30 / 70 | 20 / 80 | — | EasyKASP 30–60 %, 20–70 allowed [R] |
| `max_tm_diff` | 3 | — | 6 | [I] |
| `max_poly_x` | 5 | — | — | EasyKASP: no run of 6 [R] |
| `product_size_ranges` | [[50,120]] → effective [[61,120]] | [[50,150]] → [[65,150]] | — | EasyKASP < 100 bp, PolyMarker 50–150, LGC 50 bp flanks [R] |

**`as_pcr`**

| Param | L0 | L1 | L2 | Source |
|---|---|---|---|---|
| `opt_size` / `min_size` / `max_size` | 24 / 18 / 30 | max 32 | — | Little 1995: common primer ≈ 30 nt [R] |
| `opt_tm` / `min_tm` / `max_tm` | 60 / 57 / 63 | 55 / 65 | min 52 | [I], same targets as kasp |
| `min_gc` / `max_gc` | 30 / 70 | 20 / 80 | — | Little ≈ 50 % GC [R] |
| `max_tm_diff` | 3 | — | 6 | [I] |
| `max_poly_x` | 4 | — | — | design `pcr` preset |
| `product_size_ranges` | [[150,300]] | — | — | Little: 150–250 bp gel products [R] |

**Hard floors, applied after scoring and never relaxed:** every allele-specific primer — the designed REF primer, the derived ALT primer and, when present, their deliberate-mismatch versions — must have `tm ≥ as_min_tm` (52 °C; for a mismatch primer the floor applies to `matched_tm`) and, for derived primers, `gc ≥ as_min_gc` (15 %). A set that fails is dropped and counted in `rejected.below_floor`; an orientation whose candidates all fail is `no_sets {reason: "below_floor"}` with its `attempts[].explain` intact.

**`PRIMER_PICK_ANYWAY` is never used in a design run.** KASP runs 26 cycles at 55 °C, and `PICK_ANYWAY` returns primers at 45.4 °C [R]; such a primer gives no signal for its allele, which shows up as a *wrong homozygous call* rather than a no-call. Level 2 is the floor and it works: at the AT-rich forward orientation of 1:11109, level 2 returns `TCCGAATATAGTCATACTCTATTC` at 52.87 °C with no `PICK_ANYWAY` [V]. (`PICK_ANYWAY` *is* used in the `check_primers` **scoring** runs of §4.10, where it only prevents Primer3 from silently returning nothing for an already-chosen primer; the floors above are then applied by our own code.)

**Reporting.** Each attempt records `{level, changes, explain{left,right,pair}, pairs_returned, rejected{...}, not_scored, sets}`, with `pairs_returned` always equal to the sum of `rejected`, `not_scored` and `sets`; `changes` lists the level's effective changes, product-range minimum already raised; each set records `relaxation_level`; each orientation that needed level ≥ 1 adds `RELAXED_CONSTRAINTS`; `settings.ladder` echoes the levels and `settings.pinned` the client-set params; `settings.floors` echoes the floors.

### 4.7 Primer3 records

The record is the existing `design.buildRecord({id, seq}, {params, target: null, included: null, excluded: [], junctions: null, avoid_repeats, repeat_mask_mode})` (`design.js:436-465`) plus four server-owned tags: `SEQUENCE_FORCE_LEFT_END` **or** `SEQUENCE_FORCE_RIGHT_END`, `SEQUENCE_TARGET` (§4.4) and `PRIMER_NUM_RETURN = genotyping.num_return_per_run` (20).

Semantics relied on, all verified [R, and re-verified here V]:
- `PRIMER_FIRST_BASE_INDEX=1` is always sent (`design.js:78-91`), so forced positions are 1-based.
- `FORCE_LEFT_END = p` puts the left primer's 3′ base at `p`; `FORCE_RIGHT_END = p` puts the right primer's 3′ base at `p`, which is its **lowest** template coordinate, while `PRIMER_RIGHT_i = pos,len` reports the 5′ end; `boulder.extractPairs` (`boulder.js:155-185`) already converts that to the `[pos−len+1, pos]` footprint.
- Length still varies when only the 3′ end is forced ("considered 13" = lengths 18–30 at L0; "considered 15" at L1) [V].
- Forced candidates must still pass every constraint; failures appear in `PRIMER_LEFT_EXPLAIN` / `PRIMER_RIGHT_EXPLAIN`.

**Never sent:** `SEQUENCE_INCLUDED_REGION` (a forced end outside it is silently ignored [R]); `SEQUENCE_TARGET` or `SEQUENCE_EXCLUDED_REGION` covering any base of the zone (either rejects every forced candidate [R]); `PRIMER_PICK_ANYWAY` in a design run; any client-supplied tag. The Boulder guards (`boulder.js:36-53`: `KEY_RE`, the control-character rule, the leading `=` rule, duplicate keys) and `-strict_tags` (`primer3.js:113`) are unchanged, and a `PRIMER_ERROR` in a **design** run becomes `400 PRIMER3_INPUT_ERROR` (`design.js:650-652`). A `PRIMER_ERROR` in a `check_primers` **scoring** run (§4.10) never fails the request: it drops that candidate and increments `rejected.alt_scoring_failed` (for example `Specified right primer not in sequence` when the common primer is absent from the ALT haplotype, §4.4).

**Output guard.** Every returned pair must satisfy `left.end === t_fwd` (forward) or `right.start === t_rev` (reverse); anything else is dropped into `rejected.force`. This has never been observed [V].

### 4.8 From pairs to candidate sets

Orientations run forward, then reverse. Within an orientation the ladder runs one Primer3 design run per level and **stops at the first level that yields at least one scored set**, or when the budget of step 6 is exhausted; later levels are not run. For each pair of a level, in Primer3 rank order:
1. **Force guard** (§4.7).
2. **Zone guard** — the post-filter of §4.4, on both haplotypes (`rejected.overlap`).
3. **Common-primer 3′ neighbour** — under `neighbour_policy: avoid_3p` and a non-EMS target, reject the pair when a non-EMS neighbour lies inside the common primer's last `neighbour_3p_window` (5) bases (`rejected.common_neighbour_3p`). For an EMS target the pair is kept and the set gets the high-severity issue `NEIGHBOUR_AT_3P` instead (§4.15).
4. **Duplicates** — the same `(as_ref.target_seq, common.target_seq)` pair (`rejected.duplicate`).
5. **Scoring cap** — once `max_scored_per_orientation` (8) candidates of this orientation have been scored, counted across all its levels and including candidates that scoring later drops, every further pair that passed steps 1–4 is counted in `not_scored` and skipped.
6. **Budget reservation** — before a candidate is scored its worst-case cost is reserved:
   - **Primer3 runs:** 1 `check_primers` run, or 3 with a deliberate mismatch (the ALT run plus one mismatch run per haplotype).
   - **`ntthal` calls:** 15 when tails are used — 2 tailed hairpins, 4 tailed self `ANY`/`END1`, and 9 tailed cross calls (`ANY` once, `END1` in both directions, for each of the three primer pairings) — plus 2 duplex calls with a deliberate mismatch.

   If `primer3_runs + runs > max_primer3_runs` (54) or `thermo_calls + calls > max_thermo_calls` (272):
   - this pair and every later pair of the level are counted in `not_scored`, and the orientation's ladder stops;
   - the response gains `DESIGN_BUDGET_EXHAUSTED` (§2.15), and the other orientation still makes its design run for `explain` under the same reservation rule;
   - the sets already scored are returned in a normal `200`: the budget never produces a `500`.

   Design runs themselves are not reserved (at most 6). Actual `ntthal` calls are memoized per distinct input (`thermo.js`, §6.1), so real usage is usually far below the reservation.
7. **Scoring** (§4.9–§4.15) — a candidate dropped here counts in `rejected.alt_scoring_failed` or `rejected.below_floor`; a surviving one counts in `sets`.

At every level `pairs_returned` therefore equals the sum of `rejected`, `not_scored` and `sets`, and `sets_found` counts scored candidates only.

**Why these numbers.** The worst case — both orientations, all three levels, tails and a deliberate mismatch — is 2 × 3 design runs + 2 × 8 × 3 scoring runs = **54** Primer3 runs and 2 × 8 × 17 = **272** `ntthal` calls. Per configuration:

| Configuration | Primer3 runs | `ntthal` calls |
|---|---|---|
| `kasp`, no mismatch | ≤ 6 + 16 = 22 | ≤ 240 |
| `as_pcr` defaults (no tails, mismatch `auto`) | ≤ 54 | ≤ 32 |
| any configuration | ≤ 54 | ≤ 272 |

A design with the default caps never exhausts the budget; `DESIGN_BUDGET_EXHAUSTED` exists for lowered caps, which `scratch-fix/score_sets.js` demonstrates (§7.3). The worst case of 54 runs plus 272 calls with a pool of 4 took **1,141 ms and 1,152 ms** with the real `primer3_core` and `ntthal` (`scratch-fix/timing.js`) [V].

The example designs used, as runs/calls:

| Design | runs/calls |
|---|---|
| rs871475760 kasp | 18/123 |
| tmp_1_11193_C_T | 20/150 |
| rs5413864115 | 20/141 |
| tmp_1_11502_C_CGT | 9/66 |
| rs871475760 as_pcr | 29/4 |

`num_return_per_run` is 20 rather than 10 so that steps 2–4 can reject freely without starving an orientation (rs871475760 reverse loses 6 of its 20 pairs to step 3); step 5 then bounds the scoring work. A single-orientation request therefore returns at most 8 sets, even when `num_sets` is 9 or 10.

### 4.9 ALT primer derivation

Let `delta = len(vcf.alt) − len(vcf.ref)` and let the REF allele-specific footprint be `[s, e]` in REF template coordinates.

| Orientation | ALT footprint on `alt_seq` | ALT `matched_seq` |
|---|---|---|
| forward | `[s, e]` | `alt_seq[s..e]` |
| reverse | `[s + delta, e + delta]` | `revcomp(alt_seq[s+delta..e+delta])` |

`genomic.blocks` map each ALT-template index back to the reference with `altToRef` (§3.5); inserted bases have no coordinate and are counted in `inserted_bases`.

**Verified** [V]:

| Variant / set | REF `matched_seq` (footprint) | ALT `matched_seq` | ALT blocks | products REF/ALT |
|---|---|---|---|---|
| rs871475760 reverse | `ATCTTTGACTAGCGAGAAATTCAG` (1:11109–11132, −) | `ATCTTTGACTAGCGAGAAATTCAT` | 11109–11132 | 65 / 65 |
| rs871475760 forward | `GGTTATCCGAATATAGTCATACTCTATTC` (1:11081–11109, +) | `GGTTATCCGAATATAGTCATACTCTATTA` | 11081–11109 | 92 / 92 |
| rs5413864115 forward | `ACAGATGATTTTCCAAATGATGATTCAAA` (1:11257–11285, +) | `ACAGATGATTTTCCAAATGATGATTCAAG` | 11257–11282, 11284–11286 | 88 / 87 |
| rs5413864115 reverse | `AGAGTCTTTTCAAATTTCACACTTT` (1:11283–11307, −) | `AGAGTCTTTTCAAATTTCACACTTG` | 11282–11282, 11284–11307 | 129 / 128 |
| tmp_1_11502_C_CGT reverse | `GCAGGAAAAGAAATCCTAACATCATATG` (1:11502–11529, −) | `GCAGGAAAAGAAATCCTAACATCATATA` | 11503–11529, `inserted_bases: 1` | 72 / 74 |

`products.alt.size = products.ref.size + delta`, confirmed by `check_primers` on the ALT haplotype (87 and 74) [V]. Because sets are returned in their own `sets[]` structure and never as `pairs[]` items, no existing invariant (`product_size = right.end − left.start + 1`, `len = end − start + 1`) is broken by an indel ALT product (finding 11).

### 4.10 Scoring: what comes from where (resolves finding 9 on comparability)

| Oligo / pair | Source of `tm`, `gc`, `hairpin_th`, `self_any_th`, `self_end_th`, `end_stability`, pair `compl_*_th`, product size |
|---|---|
| `as_ref` and `common`, no deliberate mismatch | the design run itself |
| `as_alt` | one `check_primers` run on the **ALT haplotype template** |
| deliberate-mismatch primers | one `check_primers` run per haplotype, on a copy of that haplotype with the template base opposite the mismatch substituted so the primer is present verbatim; `tm` is replaced by the `ntthal` duplex against the **unmodified** own-allele footprint (`tm_method: "ntthal_duplex"`), and `matched_tm` keeps Primer3's Tm of `matched_seq` |
| tailed oligos | `ntthal` structures only (never a Tm) |

The scoring record is `PRIMER_TASK=check_primers`, `SEQUENCE_PRIMER` / `SEQUENCE_PRIMER_REVCOMP` = the set's two primers, `PRIMER_PICK_ANYWAY=1`, and the level's params with the **pair-level** constraints widened (`min_size` 15, `max_size` 36, `max_tm_diff` 30, product range `[36, upper + 100]`) so that a pair constraint cannot suppress the answer, while primer-level constraints keep their values so `PRIMER_*_PROBLEMS` stays meaningful. That string is echoed as `primer3_problems` and, when non-empty on a derived primer that still clears the floors, becomes the warn-level issue `ALT_PRIMER_SUBOPTIMAL`.

**Verified equivalence** [V]: `check_primers` on the REF haplotype reproduces the design run exactly for rs871475760 S1 (`58.721` / `34.99` / `4.5500`; `57.099` / `0.00` / `3.0200`; product 65; penalty `7.179713`). ALT runs give `56.513`, GC `33.333`, end stability `2.5700`, product 65, `PROBLEMS= Temperature too low;`; the deletion ALT gives `59.533` and product 87; the insertion ALT `58.025` (reported as 58.03, §2.1) and product 74 with the §2.10(c) S1 common primer; the −2 mismatch primer gives `59.194` (its perfect-complement Tm, which is why it is *not* reported as `tm`), hairpin `49.20`, end stability `4.3000`.

`ntthal` is always called with `-mv 50 -dv 1.5 -n 0.6 -d 50 -t 37 -r` (its own defaults are `-dv 0 -n 0`, which are not Primer3's [R]); inputs must match `^[ACGT]{1,60}$`; "No secondary structure" and negative values clamp to 0 (a real case: FAM-REF × common `ANY` = −11.46 → 0) [V]. `oligotm` is not used anywhere: every nominal Tm comes from Primer3 itself, which keeps REF, ALT and common numbers on one scale.

**Never reported, with reasons:** the duplex Tm of a primer on the *other* allele and any REF-vs-ALT ΔTm — `ntthal` leaves a 3′-terminal mismatch unpaired as a dangling end, the measured own-allele advantage was only 0.06–0.55 °C, and in one orientation it had the wrong sign [R]; the `oligotm` Tm of a mismatch primer's ordered sequence, because it is the Tm of its perfect complement.

### 4.11 The deliberate mismatch (resolves finding 5)

Applies when `deliberate_mismatch: "auto"` (the `as_pcr` default). Both allele-specific primers get one extra mismatch at `mismatch_position` (2 by default) from the 3′ end; the common primer never does.

**Source.** Little, S. (1995), *Current Protocols in Human Genetics* 9.8, Tables 9.8.1 and 9.8.3 [R]. The implementation is a **literal 32-cell lookup**, not a class shortcut: the shortcut differs from the table in 2 of 32 cells, and one of them changes a real answer for tmp_1_11193_C_T [V].

- **Row** = the primer's own terminal mispair, unordered `{primer 3′ base, template base opposite it on the other allele}`. For a forward primer the template base is the complement of the other allele's plus-strand base; for a reverse primer it is the plus-strand base itself.
- **Column** = the template base opposite the mismatch position, i.e. `complement(original primer base at −k)`.
- Table 9.8.3 classes: maximum `{A,G} {C,T} {T,T}`; strong `{C,C}`; medium `{A,A} {G,G}`; weak `{A,C} {G,T}`.

| Terminal mispair (row) | class | template A | G | C | T |
|---|---|---|---|---|---|
| AA, GG | medium | A | G | A | G |
| AG, TT, CT, CC | max / strong | C | T | A | G |
| AC | weak | G | A | C | T |
| GT | weak | G | A | T | C (Little prints "C or T"; C is chosen — it is maximum in Little and strong in Ye 2001, so it survives their disagreement) |

**Verified outputs** [V] — these are pinned as unit tests:

| Variant | orientation | −2 REF | −2 ALT | −3 REF | −3 ALT |
|---|---|---|---|---|---|
| rs871475760 | reverse | A→G `ATCTTTGACTAGCGAGAAATTCGG` | A→G `ATCTTTGACTAGCGAGAAATTCGT` | — | — |
| rs871475760 | forward | T→C `TGGTTATCCGAATATAGTCATACTCTATCC` | T→C `TGGTTATCCGAATATAGTCATACTCTATCA` | — | — |
| tmp_1_11193_C_T | reverse | A→**C** `ACAACTTTTTAATATATTGTGTATACTCTCG` | A→**T** `ACAACTTTTTAATATATTGTGTATACTCTTA` | T→G `…ACTCGAG` | T→G `…ACTCGAA` |
| tmp_1_11193_C_T | forward | T→G `ACAAAGATAGATAACAAAAATAGCTCGC` | T→G `ACAAAGATAGATAACAAAAATAGCTCGT` | C→A `ACAAAGATAGATAACAAAAATAGCTATC` | C→A `ACAAAGATAGATAACAAAAATAGCTATT` |

The 11193 reverse row is where a class shortcut goes wrong: it would give C for both alleles, while the table gives C for REF and T for ALT. Each primer is looked up on its own row, so the two alleles can legitimately receive different bases.

**Position.** −2 is the default, because Table 9.8.1 is defined for the penultimate base [R]. −3 is selectable and reported as `source: "Little 1995 Table 9.8.1 applied at -3 (extrapolated)"`: Liu et al. 2012 found the third base from the 3′ end best empirically (32.9 % vs 25.1 % polymorphic), but their pair notation is undefined, so their result cannot be turned into a table [R]. The default is listed for the user to confirm (§11). Either position leaves the wrong-allele primer with `mm_pos [2,1]` or `[3,1]`, which `classify.js` calls `unlikely` [V, R].

**Guard.** For an indel, if the base at −k differs between the haplotypes, no mismatch is applied and the set gets warning `MISMATCH_NOT_APPLICABLE {position}` — otherwise the "extra" mismatch would itself be allele-specific.

**Scoring** follows §4.10. Verified for rs871475760 reverse [V]: own-allele duplex 54.69 °C (REF) and 54.58 °C (ALT) against 57.04 / 56.45 unmodified; hairpins 49.20 and 50.62 → `MISMATCH_STRUCTURE`; own allele `mm_pos [2]` → `likely`; other allele `[2,1]` → `unlikely`.

### 4.12 Discrimination analysis

For each allele-specific primer, `target_seq` is aligned with the check's own 3′-anchored DP, `realign.realignSiteOnWindow(site, seq, winStart, three_prime_window, regionLength)` (`realign.js:196-221`), against its own haplotype window `[fp.start − 3, fp.end + 3]` and against the other haplotype (forward: the same indices; reverse: shifted by `∓delta`). The alignment is classified with `classify.classifyAmplicon` against a perfect partner, using the check's default params (`classify.js:21-30`), so design and check speak one vocabulary.

Reported per primer: `discrimination {own_allele {mm_pos, likelihood}, other_allele {mm_pos, likelihood}, terminal_mismatch_class, in_shift_tract}`.

- All KASP allele-specific primers of the four example variants give other-allele `[1]` / `likely_weak`, including the homopolymer deletion — the DP prefers a terminal mismatch over a gap (`realign.js:98-101`: diagonal, then up, then left) [V, R].
- With a −2 mismatch: `[2,1]` / `unlikely` [V].
- `terminal_mismatch_class` is the only in-silico discrimination signal the literature supports for a tail-only KASP assay: A/G and C/T sites are weak+weak, A/C and G/T maximum+maximum, A/T medium+maximum, C/G strong+medium [R]. Both classes weak adds the info-level issue `WEAK_TERMINAL_CLASS`.
- `in_shift_tract` is decided per allele-specific primer from the variant's shift tract (§3.5) and the reference coordinate of the primer's 3′ base **on its own haplotype**:
  - `as_ref`: true when `discriminating.<orientation>.position` lies inside the tract;
  - `as_alt`: true when `discriminating.<orientation>.alt_maps_to` lies inside the tract, or is `null` (an inserted base) while the variant has a tract;
  - false for both when the variant has no tract.

  A primer that ends inside the repeat can meet the other allele through a one-unit bulge with its 3′ base paired, at the same edit cost as the terminal mismatch that `realign.js` reports (its DP prefers the mismatch, `realign.js:98-101`) [I]; TP-ARMS assays of homopolymer indels succeeded in under half of cases [R]. Each such primer adds one `SHIFT_TRACT_DISCRIMINATION`, and on the other allele's genomes the check reports it `uncertain` rather than "does not amplify" (§5.7). For a deletion only `as_ref` can be inside the tract, for an insertion only `as_alt` (§3.5). Verified [V]: rs5413864115 S1 `as_ref` true (11285) and `as_alt` false (11286); S2 `as_ref` true (11283) and `as_alt` false (11282, the anchor).
- `WEAK_DISCRIMINATION` (high severity) fires when the other-allele likelihood is `likely`.

### 4.13 Tm balance

- `tm_balance.as_tm_diff = |as_ref.tm − as_alt.tm|`, computed exactly from the unrounded values and rounded once (§2.1). Above `as_tm_diff_warn` (1.0 °C, EasyKASP [R]) → `AS_TM_IMBALANCE` (high above 2.0 °C). Verified: tmp_1_11502_C_CGT S1 = 1.16 °C [V].
- `tm_balance.common_minus_as = common.tm − max(nominal Tm of the two AS primers)`, where "nominal" is `matched_tm` for a mismatch primer and `tm` otherwise, so the two sides of the comparison are always Primer3 Tm values. It too is computed exactly from the unrounded values (§2.1): tmp_1_11193_C_T S1 is 57.102 − 56.427 = 0.675, reported 0.68 [V]. Outside `[common_tm_low, common_tm_high]` = −1.0 … +3.0 °C → `COMMON_TM_OUT_OF_RANGE` (EasyKASP puts the common primer ≈ 2–3 °C above the allele-specific primers [R]).
- v1 does **not** move a 5′ end to rebalance Tm. Doing so would give the two pairs of one set different `expected` windows, which the check's set validation deliberately forbids (§5.2), and it complicates the order sheet. It is listed as a deferred improvement (§9.3); a prototype that lengthened one 5′ end by 1 nt took tmp_1_11502_C_CGT from 1.16 to 0.13 °C, so the gain is real [I].

### 4.14 Tails and tail-aware checks

- FAM `GAAGGTGACCAAGTTCATGCT`, HEX `GAAGGTCGGAGTCAACGGATT` [R]. `ref_fam_alt_hex` (the kasp default) puts FAM on `as_ref`; `ref_hex_alt_fam` swaps them. The default assignment is listed for the user to confirm (§11).
- A tail appears only in `tail_seq` and `order_seq`. It never reaches Primer3 (built-in 36 nt maximum [R]) and never reaches the check (`^[ACGT]{15,36}$`).
- Structures come from `ntthal` at 37 °C: per oligo `tailed {hairpin_th, self_any_th, self_end_th}`, and per set `thermo.tailed {ref_alt_any_th, ref_alt_end_th, ref_common_any_th, ref_common_end_th, alt_common_any_th, alt_common_end_th}`, where each `*_end_th` is `max(END1(x,y), END1(y,x))`.
- **Thresholds.** A tailed **hairpin** warns at ≥ 47 °C and is **high** severity at ≥ 55 °C, the final KASP annealing temperature [R]; tailed **dimers** (`ANY`, `END1`) warn at ≥ 47 °C, Primer3's own default. Every warn issue adds 1 to the score and every high issue adds 3 (§4.16), so a warn-level tailed hairpin moves a set only slightly. This matters because the HEX tail alone lifts hairpins above 47 °C on ordinary sequences — HEX + REF and HEX + ALT both give 51.50 °C, identical, so the allele base plays no part [V].
- Verified tailed values [V]: FAM+`ATCTTTGACTAGCGAGAAATTCAG` 42.62 / 3.63 / 11.64; HEX+`…TTCAT` 51.50 / 26.85 / 11.47; FAM+`GGTTATCC…TTC` 36.16 / 3.63 / 12.01; HEX+`…TTA` 66.22 / 27.34 / 10.26; both tails alone 0.00; cross-dimers `ANY` 6.89 (S1 REF×ALT), 5.46 (S1 ALT×common), 5.46 (S2 REF×ALT), 0 (S2 ALT×common); `END1` for S1 REF×ALT is −12.20 in one direction and 8.21 in the other, which is why each `*_end_th` is the maximum of both directions (reported 8.21).

### 4.15 Neighbouring known variants (resolves the EMS half of finding 5's neighbour rule)

"Neighbours" are variant keys inside the template window other than the target key and its merged ids. They come from the same cached 10 kb chunks as §3.3, so no per-variant fetch ever happens (finding 8).

**`distance_from_3p` is measured along the primer, from that primer's own 3′ end** — 1 is the base under the 3′-terminal base — taking the minimum over the neighbour's `minimal` span. For a forward set the common primer is the right primer, so its 3′ end is `genomic.start`; for a reverse set it is the left primer and its 3′ end is `genomic.end`. (A formula keyed on the *set* orientation measures from the 5′ end and reports, for example, 2 instead of 23; that error is explicitly not made here.)

| Neighbour | in an AS primer's last 5 nt | in the common primer's last 5 nt | elsewhere under a primer | inside the dense window |
|---|---|---|---|---|
| non-EMS, natural target | orientation `blocked`, no Primer3 run | pair rejected (`rejected.common_neighbour_3p`) | reported; warn issue `NEIGHBOUR_IN_PRIMER` | counted |
| non-EMS, **EMS target** | kept, high-severity issue `NEIGHBOUR_AT_3P` (+3 to the score, quality `poor`) | kept, high-severity issue | warn issue | counted |
| EMS | listed in `neighbours[]` with `ems: true`; never blocks; no issue | never blocks; no issue | listed; no issue | not counted |

EMS mutations are private to mutant lines in the BTx623 background, so a neighbouring EMS site almost never segregates in the panel a breeder genotypes. But an EMS *target* is itself routinely crossed to non-BTx623 parents, where natural neighbours do segregate — so for an EMS target natural neighbours are still detected and reported, they merely stop blocking the only orientation that can work. `neighbour_policy: "ignore"` reports everything and blocks nothing.

**One issue per primer site.** `NEIGHBOUR_IN_PRIMER` and `NEIGHBOUR_AT_3P` are raised at most once per code for the **allele-specific site** and at most once per code for the **common primer**. The two allele-specific primers share every base except the discriminating end, so a neighbour under either of them sits under the same genomic site for both alleles. The allele-specific issue therefore carries `role: "as_ref"` (standing for the pair) and lists the union, by key, of both primers' non-EMS neighbours, each with its smaller `distance_from_3p`, sorted by distance. Each primer's own `neighbours[]` still lists every neighbour under it, and `neighbour_sites` counts distinct non-EMS keys over the three primers, so a set carries at most four neighbour issues. Verified [V]: tmp_1_11193_C_T S1 has one allele-specific issue (`rs873026643` at 12) and one common issue (`rs5413864115` at 16); rs5413864115 S1 has a single common issue listing three neighbours.

`DENSE_NEIGHBOURS {count, window: 30}` fires above `dense_count` (2) non-EMS neighbour keys whose `minimal` span overlaps the dense window `[vcf.position − 30, vcf.position + len(vcf.ref) − 1 + 30]` (the target excluded) — LGC's FAQ limit, kept as a warning because 105 of 172 SNPs in the research window exceed it [R]. Verified counts [V]: rs871475760 0, tmp_1_11193_C_T 1, rs5413864115 4, tmp_1_11502_C_CGT 5. Verified blocking [V]: tmp_1_11502_C_CGT forward is blocked by `rs5413863234` (1:11500 G/A, EVA) at distance 4.

The full neighbour list is not repeated in the design response; the UI draws its track from `GET /primers/variants` over the template window, which is a cache hit. Degenerate/IUPAC bases inside primers are not attempted: Primer3 never places a primer over an IUPAC base, even with `PRIMER_MAX_NS_ACCEPTED=1` [R] (§9.3).

### 4.16 Quality, score, ranking and limits

Every set carries `issues[]` with severities, a `quality` and a `score` (lower is better). The score is the exact sum of the six terms S1–S6, computed and rounded by S7:

| Rule | Term | Definition |
|---|---|---|
| S1 | penalty | `PRIMER_PAIR_i_PENALTY` of the candidate's **own** design-run pair `i` — the Primer3 pair that supplied its `as_ref` and `common` primers — as printed (6 decimals). It is `PRIMER_PAIR_0_PENALTY` only when `i` is 0: rs871475760 S2 is Primer3's pair 2 (14.634369), while pair 0 of that run has 14.376148. |
| S2 | relaxation | `5 × relaxation_level` |
| S3 | allele-specific Tm balance | `2 × max(0, as_tm_diff − as_tm_diff_warn)`, with `as_tm_diff` from §4.13 and `as_tm_diff_warn` 1.0 °C |
| S4 | common-primer Tm | `2 × d`, where `d` is how far `common_minus_as` (§4.13) lies outside `[common_tm_low, common_tm_high]` = [−1.0, +3.0], and 0 inside |
| S5 | product size | `kasp` only, at every level: `max(0, products.ref.size − 100) / 25` |
| S6 | issues | `3 ×` the counted high-severity issues `+ 1 ×` the counted warn-severity issues |
| S7 | arithmetic | the exact decimal sum of S1–S6 (§2.1: integers in units of 10⁻⁶, never binary floating point); `score` is that sum rounded half away from zero to 2 decimals, and ordering uses the exact sum |

**S6 — which issues exist and which count.** Each issue counts once, by its severity, with two exclusions: `AS_TM_IMBALANCE` and `COMMON_TM_OUT_OF_RANGE`, which S3 and S4 already price, and every `info` issue (`WEAK_TERMINAL_CLASS`). How many issues of one code a set can carry is fixed:

| Code | At most one per | Rule |
|---|---|---|
| `NEIGHBOUR_IN_PRIMER`, `NEIGHBOUR_AT_3P` | primer site: the allele-specific pair, and the common primer | §4.15. The allele-specific issue has `role: "as_ref"` and lists the union, by key, of both allele-specific primers' non-EMS neighbours, each with its smaller distance. EMS neighbours raise no issue. |
| `ALT_PRIMER_SUBOPTIMAL` | ordered derived oligo | From that oligo's own echoed `primer3_problems` (§4.10). `as_alt` without a deliberate mismatch uses its ALT-haplotype run; with a deliberate mismatch, `as_ref` and `as_alt` each use their own mismatch run. The matched ALT primer's ALT run then never raises the issue: that sequence is not ordered, and its Tm is already floored through `matched_tm`. |
| `SHIFT_TRACT_DISCRIMINATION`, `WEAK_DISCRIMINATION` | allele-specific primer | §4.12 |
| `TAILED_STRUCTURE` | tailed oligo × metric (`hairpin_th`, `self_any_th`, `self_end_th`), and tailed pairing × metric (the six `thermo.tailed` values) | §4.14 thresholds |
| `MISMATCH_STRUCTURE` | deliberate-mismatch primer × metric (`hairpin_th`, `self_any_th`, `self_end_th`) | §2.15 thresholds |
| `MISMATCH_NOT_APPLICABLE` | set | §4.11 |

So `ALT_PRIMER_SUBOPTIMAL`, `NEIGHBOUR_IN_PRIMER`, `SHIFT_TRACT_DISCRIMINATION`, `MISMATCH_STRUCTURE` and a warn-level `TAILED_STRUCTURE` each add 1, and a high `TAILED_STRUCTURE`, `NEIGHBOUR_AT_3P` or `WEAK_DISCRIMINATION` adds 3. The re-verifier's reading switches show what each rule decides [V] (`scratch-reverify/rescore_oligo.log`, `rescore_suboptall.log`, `rescore_pair0.log`):
- one neighbour issue per oligo would make tmp_1_11193_C_T S1 23.51 and tmp_1_11502_C_CGT S1 10.56;
- counting the matched ALT run would make rs871475760 AS-PCR S1 8.92 and promote the forward set `f2adeabf4265` (8.78);
- a pair-0 penalty would tie six reverse rs871475760 candidates at 9.18.

`quality` is `poor` when any issue (the two uncounted Tm issues included) is high severity, `usable` when there is a warn issue or `relaxation_level ≥ 1`, and `good` otherwise. Quality does not change the ordering — a `poor` set outranks a `usable` one whose score is higher — and no set is ever hidden.

**Ordering.** Within an orientation, candidates sort by the exact score (S7), then by `key` as a stable tie-break. Across orientations: the best set of each orientation comes first (ordered between themselves the same way), then all remaining sets by exact score and key. The list is cut to `assay.num_sets`, and `id` = `S{rank+1}` is assigned afterwards.

**Worked examples** [V]. Every score, rank and key of §2.9 and §2.10 is recomputed from S1–S7 by `scratch-fix2/score2.js`, on real Primer3 2.6.1, `ntthal` and live-Ensembl neighbours. The re-verifier's independent `rescore.js` gives the same values with the §3.5 tract:

| Variant | Set | Key | S1 penalty | S2 5 × level | S3 Tm diff | S4 common | S6 high × 3 | S6 warn × 1 | S5 product | **score** | exact sum |
|---|---|---|---|---|---|---|---|---|---|---|---|
| rs871475760 kasp | S1 reverse | `f9df650ad116` | 7.1797 | 0 | 0 | 0 | 0 | 2: `ALT_PRIMER_SUBOPTIMAL` (as_alt), `TAILED_STRUCTURE` (as_alt hairpin 51.50 °C) | 0 (65 bp) | **9.18** | 9.179713 |
| | S2 forward | `1accc54c262d` | 14.6344 | 0 | 0 | 0 | 1: `TAILED_STRUCTURE` (as_alt hairpin 66.22 °C) | 2: `ALT_PRIMER_SUBOPTIMAL` (as_alt), `NEIGHBOUR_IN_PRIMER` (common) | 0 (92 bp) | **19.63** | 19.634369 |
| tmp_1_11193_C_T | S1 forward | `f5a5c6f2ebcf` | 14.4705 | 5 | 0 | 0 | 0 | 2: `NEIGHBOUR_IN_PRIMER` (allele-specific; common) | 1.04 (126 bp) | **22.51** | 22.510497 |
| | S2 reverse | `c791a4956fd3` | 18.2784 | 5 | 0 | 0 | 0 | 1: `ALT_PRIMER_SUBOPTIMAL` (as_alt) | 0.56 (114 bp) | **24.84** | 24.838417 |
| rs5413864115 | S1 forward | `a5277232d8ab` | 8.4793 | 5 | 0 | 0 | 0 | 2: `NEIGHBOUR_IN_PRIMER` (common), `SHIFT_TRACT_DISCRIMINATION` (as_ref) | 0 (88 bp) | **15.48** | 15.479320 |
| | S2 reverse | `cb3ef66afd37` | 14.7932 | 5 | 0 | 0 | 0 | 3: `NEIGHBOUR_IN_PRIMER` (allele-specific; common), `SHIFT_TRACT_DISCRIMINATION` (as_ref only: as_alt ends on the anchor 11282, outside the tract) | 1.16 (129 bp) | **23.95** | 23.953151 |
| tmp_1_11502_C_CGT | S1 reverse | `53942cb55348` | 7.2375 | 0 | 0.32 (1.16 °C) | 0 | 0 | 2: `NEIGHBOUR_IN_PRIMER` (allele-specific; common); `AS_TM_IMBALANCE` not counted | 0 (72 bp) | **9.56** | 9.561494 |
| rs871475760 as_pcr | S1 reverse | `7f9af6b1c938` | 2.9227 | 0 | 0 | 0 | 0 | 5: `ALT_PRIMER_SUBOPTIMAL` (as_ref, as_alt: mismatch runs), `MISMATCH_STRUCTURE` × 2, `NEIGHBOUR_IN_PRIMER` (common) | — | **7.92** | 7.922723 |

`WEAK_TERMINAL_CLASS` is present on tmp_1_11193_C_T S1 and S2, rs5413864115 S1 and tmp_1_11502_C_CGT S1, and counts 0.

The eight scored forward candidates of rs871475760 score, in order:
- 19.63 — `1accc54c262d`, Primer3's pair 2, the top forward set
- 19.79, 19.93, 21.53, 21.67, 21.76 and 21.90
- 22.30 — Primer3's pair 0, the 30-mer `TGGTTATCCGAATATAGTCATACTCTATTC` with two ≥ 55 °C tailed hairpins

So the formula decides the order, not the Primer3 penalty.

| Limit | Value |
|---|---|
| Sets returned | `num_sets` ≤ 10 (default 6); at most 8 from one orientation |
| Pairs per design run | `num_return_per_run` 20 |
| Candidates scored | ≤ `max_scored_per_orientation` (8) per orientation, across levels (§4.8) |
| Primer3 runs | ≤ 6 design runs + ≤ 48 scoring runs = `max_primer3_runs` (54); a candidate whose reservation does not fit stays unscored (`not_scored`) with warning `DESIGN_BUDGET_EXHAUSTED`, never a `500` |
| `ntthal` calls | ≤ `max_thermo_calls` (272), pool of 4, 5 s each |
| Ensembl calls | ≤ 2 lookups + ≤ 6 overlap chunks, 8 s each, all before the semaphore |
| FASTA reads | 2 |
| Deadline | 45 s end to end |

Measured cost per call: Primer3 forced design 16–24 ms, `check_primers` 9 ms, `ntthal` 1.9–5.1 ms [R]. The worst case the caps allow — 54 Primer3 runs (6 design + 48 `check_primers`) and 272 `ntthal` calls with a pool of 4 — took 1,141 ms and 1,152 ms with the real binaries (`scratch-fix/timing.js`) [V]. A design therefore stays well under 1.5 s inside its slot.

### 4.17 Order rows, labels and the submission string

Three order rows per set, in the order REF, ALT, common, with the fields of `PrimerGenotypingOrderRow`. `notes` carries `deliberate mismatch −2 A→G (Little 1995 Table 9.8.1)` for mismatch primers and `tail hairpin 66.2 °C` (the hairpin Tm rounded to 1 decimal, §2.1) whenever the primer has a `TAILED_STRUCTURE` hairpin issue of either severity (§2.9: 51.5 °C, warn, on the S1 ALT row; 66.2 °C, high, on the S2 ALT row).

`assay.kasp_mix` (kasp only) is `{stock_uM: 100, as_ref_uL: 12, as_alt_uL: 12, common_uL: 30, water_uL: 46, total_uL: 100}` — the 12:12:30 ratio at 100 µM is from Makhoul et al. 2020 [R]; the 46 µL of water to reach 100 µL is arithmetic [I].

`variant.submission_sequence` is 50 bp of reference each side of the `minimal` alleles, written `[REF/ALT]` with `-` rendered as an empty string (LGC's `[/T]` notation [R]). Biallelic non-EMS SNV neighbours inside the flanks become IUPAC codes (AG→R, CT→Y, CG→S, AT→W, GT→K, AC→M); indel and multi-allelic neighbours stay as reference bases and add `SUBMISSION_NEIGHBOURS_OMITTED {ids}`. A `*` allele (a spanning deletion, §3.4) is ignored when deciding whether a neighbour is biallelic, so `rs5413863901` C/T/`*` at 1:11318 is written `Y` in the rs5413864115 string below and adds no warning. All four strings were regenerated from the genome plus live Ensembl while writing this spec [V]:

| Variant | `submission_sequence` |
|---|---|
| rs871475760 | `YCCTCAAAAAGCTTCTCTAAGTGGTTATCCGAATATAGTCATACTCTATT[C/A]TGAATTTCTCGCTAGTCAAAGATAACAAAAATAGCATATTCTGGATTTCT` |
| tmp_1_11193_C_T | `GCATATTCTGGATTTCTCWGTAGACAAAGATAGATAACARAAATAGCTCT[C/T]TAGAGTATACACAATATATTAAAAAGTTGTTAGAGAGTGAAAATATATAG` |
| rs5413864115 | `AAAATATATAGAAAACAATTTTATACAGATGATTTTCCAAATGATGATTC[A/]AAGTGTGAAATTTGRAAAGWCTCTTRGASATGMTYTAAGTGGAAGGAACA` |
| tmp_1_11502_C_CGT | `ATGTTAGGATCTTTGCAACCCWGTGTTGCGTGCAATCTCGGTATCTCRCC[/GT]ATATGAYGTTAGGWTTTCTTWTCCTGCAACTGCCAAGAGAAYAAATATAT` |

### 4.18 The proposed check request

`check.request` is a ready-to-POST `PrimerCheckRequest`: `system_name` from the design, `mode: "region"`, `checks: ["specificity","pangenome"]`, no `genomes`, no `params`. Sets are added in rank order while the totals stay inside `check_max_sets` (5), 10 pairs and `check_max_unique_primers` (13, the largest count whose full-panel estimate stays under 6,000 CPU-s [V]); the rest are listed in `check.omitted_set_ids`. `genotyping` is `{variant: variant.vcf plus region, sets: the included sets' check.set}`. The client may post it verbatim or rebuild it from any subset with `buildGenotypingCheckRequest` (§8.6).
---

## 5. Check extension

### 5.1 Request additions and validation order

The schema is §2.11. `check/normalize.js` gains these steps (★ = new); line numbers are from the `primer-design` worktree [V]:

| # | Step | Where | Error |
|---|---|---|---|
| 1 | existing `validateShape`, plus ★ `genotype.validateShape(body.genotyping)` mirroring the swagger definition (keys, types, patterns, 1–5 sets) | `normalize.js:83-143` | 400 `INVALID_REQUEST {field}` |
| 2 | existing: unique pair ids, uppercase primers, ≤ 20 distinct primers | `normalize.js:310-329` | existing |
| 3 | ★ mode: `genotyping` requires `mode` ∈ {`gene`, `region`}; other modes do not honour `expected` (`specificity.js:80-82`) | after `normalize.js:300` | 400 `INVALID_REQUEST {field: "genotyping", reason: "mode"}` |
| 4 | ★ structural link checks (pure, §5.2 part A) | after the pairs map (`normalize.js:310-323`) | 400 `GENOTYPING_SET_INVALID` |
| 5 | existing expected/params rules, catalog, gene lookup, reference resolve | `normalize.js:335-411` | existing |
| 6 | ★ variant and positional validation, one FASTA read (§5.2 part B) | after `normalize.js:404-411` | 400 `REF_MISMATCH`, `VARIANT_TOO_REPETITIVE`, `GENOTYPING_SET_INVALID`; 404 `UNKNOWN_REGION`; 400 `REGION_OUT_OF_BOUNDS`; 422 `NO_SEQUENCE` |
| 7 | existing pan-genome genome resolution | `normalize.js:413-453` | existing |
| 8 | ★ normalized block and the cost estimate with `genotyping: true` | `normalize.js:455-461` | 422 `JOB_TOO_LARGE` |

Steps 1–4 are pure and run before any catalog lookup, so the contract test's "stop at the catalog stub" pass still works.

### 5.2 Submit-time set validation (resolves the check half of finding 1 and of reuse's set-validation gap)

**Part A — structure (pure).** For each set, in order: `ref_pair` and `alt_pair` exist (`unknown_pair`) and differ (`same_pair`); neither pair is used by another set (`pair_reused`); both carry `expected` (`expected_required`) with identical values (`expected_differs`); the two pairs share **exactly one** primer sequence, on the same side (`no_shared_common`). The shared side is the common primer; the other side holds the allele-specific primers. **Orientation and the allele side are derived here, never taken from the client**: a shared left primer means the allele-specific primers are the right primers, i.e. `reverse`; a shared right primer means `forward`.

**Part B — position (one FASTA read of ≤ 3.5 kb).**

```
fasta = resolved.assemblies[system_name].fasta.dna
lo = max(1, position − 1,700); hi = min(regionLength, position + len(ref) + 1,700)
seq = sequence.fetch(fasta, region, lo, hi, 1)
ref check: seq at [position, position+len(ref)−1] == upper(ref)     else 400 REF_MISMATCH
v = leftNormalize(variant); d = discriminating(v); zone = zone(v)    # the same pure code as design
H_ref, H_alt = the two haplotype windows over [lo, hi]
for each set:
  as_ref = the REF pair's allele-specific primer; as_alt = the ALT pair's
  face   = orientation === 'forward' ? 'F' : 'R'
  a = realign.realignSiteOnWindow({key, primer: as_ref, face, p3: index of d.<orientation> in H_ref}, H_ref, lo, 5, regionLength)
  b = realign.realignSiteOnWindow({key, primer: as_alt, face, p3: matching index in H_alt},          H_alt, lo, 5, regionLength)
  each alignment must have its 3′ base at the anchor        else reason not_at_variant
  each may carry at most one extra mismatch, at distance 2 or 3, and no indel ops   else too_many_edits
  a must fit the REF haplotype better than the ALT one, and b the reverse           else alleles_swapped
  the common primer's footprint must not intersect [zone.start − 1, zone.end + 1]   else common_in_zone
  the 5′ genomic end of as_ref must equal expected.start (forward) or expected.end (reverse)  else expected_mismatch
```

A single extra mismatch at distance 2 or 3 is accepted **without being declared** — that is a deliberate-mismatch primer — and its position is recorded in `results.genotyping.sets[].deliberate_mismatch_positions` so the prediction can subtract it (§5.7). Every failure is `400 GENOTYPING_SET_INVALID {set_id, reason, pair_ids, primer?, sequence?, mm_pos?}`.

Rejecting at submit matters because a mislabelled or misplaced set would otherwise queue a job of up to 6,000 CPU-s and return confident, meaningless allele predictions for 120 genomes. The cost is one ≤ 3.5 kb FASTA read plus two DP alignments per set, under 5 ms [I].

### 5.3 Job ids and `ALGORITHM_VERSION` (resolves finding 10)

- **`request.genotyping` is stored exactly as the client sent it**, with only value normalization (uppercased alleles, left-aligned variant). No server-derived key — no `algorithm`, no `orientation` — is added, so `job.request` still validates as a `PrimerCheckRequest`. Orientation is re-derived by the worker from the shared primer (a pure function of the stored pairs).
- **Requests without `genotyping` keep their exact job ids**: `canonicalJSON` skips undefined values (`jobs/index.js:23-42`), and the key is simply absent.
- **The genotyping algorithm version enters the job id through `jobId`'s `algo` argument.** `check.algorithmVersionFor(request)` returns `'2'` for plain requests and `'2+g1'` when `request.genotyping` is present. The call site (`jobs/index.js:164`) becomes
  `const algo = typeof check.algorithmVersionFor === 'function' ? check.algorithmVersionFor(norm.request) : check.ALGORITHM_VERSION;`
  The `typeof` guard is required: `job_id.test.js:183` injects a check-module stub that has only `ALGORITHM_VERSION` and `normalize`, and an unguarded call would throw there [V].
- **`ALGORITHM_VERSION` stays `'2'`** (`check/index.js:12`, with the fallback copy at `run.js:51`). Nothing changes for existing requests: no new results key, and BLAST, sites, amplicons, classification and verdicts are untouched. `job_id.test.js:95-102`, the docs and the swagger descriptions stay as they are.
- Results report `results.engine.algorithm_version: "2"` and `results.genotyping.algorithm_version: "g1"`. Changing the caller or the prediction rules bumps only the `g` part, which invalidates genotyping jobs alone.

### 5.4 Cost model

```
cpu_s = ceil( existing terms + (genotyping ? genome_tasks × genotype_cpu_s_per_genome : 0) )
genome_tasks = 1 + number of pan-genome genomes        # genotyping requires a genome target
```

`breakdown.genotyping` is added **only when genotyping is requested**, so `cost.test.js:52`, which deep-equals the four existing breakdown keys, still passes [V]. The coefficient 0.2 CPU-s per genome covers the anchor fetch plus DP (< 0.01 CPU-s) and, when needed, one megablast of a ≈ 430 bp segment (0.10–0.15 CPU-s) [R].

Computed with the server's own `cost.js` and the real genome sizes (reference 708,735,318; 119 others summing 83,774,241,795) [V]:

| Check | Without genotyping | With genotyping |
|---|---|---|
| 1 set (3 primers), full panel | 1,333 | 1,357 |
| 2 sets (6 primers), full panel | 2,666 | **2,690** |
| 4 sets (12 primers), full panel | 5,332 | 5,356 |
| 13 primers, full panel | 5,776 | **5,800** (largest accepted) |
| 14 primers, full panel | 6,221 | 6,245 → `422 JOB_TOO_LARGE` |
| 2 sets (6 primers), 3 genomes | 95 | **95** (the genotyping term does not change the rounded value) |

Sets that share an allele-specific primer but differ in the common primer cost 4 primers for 2 sets, so `check_max_unique_primers` is 13, not a set count.

### 5.5 Allele-aware interpretation of the existing blocks

`specificity` and `pangenome` are unchanged, and on their own they do not answer allele questions:

| Pair | on a REF genome | on an ALT genome |
|---|---|---|
| REF pair, no deliberate mismatch | `likely`; verdict `specific`; `single_perfect` | `likely_weak` (AS `mm_pos [1]`); `single_mismatch`, `terminal_mismatch: true` |
| ALT pair, no deliberate mismatch | `likely_weak`; `single_mismatch` | `likely`; `single_perfect` |
| REF pair with a −2 mismatch | `likely` (`[2]`); `single_mismatch`, not terminal | `unlikely` (`[2,1]`); `no_amplicon` with `nearest` |
| ALT pair with a −2 mismatch | `unlikely`; reference verdict `on_target_missing` | `likely` (`[2]`); `single_mismatch` |

Rules for consumers: allele and amplification questions are answered **only** by `results.genotyping`. `sets[].specificity.{ref_pair,alt_pair}.consistent_with_allele` states explicitly whether an apparently poor verdict is the expected consequence of allele specificity (`match` → an amplifying verdict; `blocked` → `on_target_missing`), and `off_target_count` is still reported.

### 5.6 The per-assembly allele caller (`check/genotype.js`) — resolves finding 1

**Prepared once per job**, inside `referenceStage`:
- re-normalize the variant against the reference; compute `shift`, the zone and the shift tract (§3.5), and for each set the `in_shift_tract` flag of both allele-specific primers from the set's derived orientation (§4.12);
- `K = max(genotype_flank_min (15), (zone.end − zone.start + 1) + max(len(ref), len(alt)))`;
- **core** `C_ref = ref[vcf.position − 1 .. vcf.position + len(vcf.ref) − 1 + shift + 1]` — the left-aligned VCF span **with its anchor base**, extended right by `shift`, plus one flanking base on each side — and `C_alt` = the same reference span read on the ALT haplotype, i.e. with the inserted bases added or the deleted bases removed. The core is a span of its own, not the shift tract of §3.5, which excludes the anchor and serves only `in_shift_tract`. For an indel with a tract the core is `[tract.start − 2, tract.end + 1]` = `[zone.start − 1, zone.end]`, so every base that differs between the alleles, or can slide, lies strictly inside it. The zone ±1 is not used because it runs one base further (`TCAAAGT` at rs5413864115);
- **haplotypes** `H_ref` / `H_alt` = the reference span `[core.start − K, core.end + K]` read on each haplotype, so the two differ only inside the core;
- per set, the reference segment `R_s = ref[expected.start − 50 .. expected.end + 50]`.

Verified values [V], computed by `scratch-fix/cores.py` and, with the shift tract, by `scratch-fix2/tract56.py` on the genome read from the 50111 API:

| Variant | shift | shift tract (§3.5) | core span | `C_ref` / `C_alt` | K | haplotype span (`H_ref` / `H_alt` length) |
|---|---|---|---|---|---|---|
| rs871475760 | 0 | none | 1:11108–11110 | `TCT` / `TAT` | 15 | 1:11093–11125 (33 / 33 nt): `ATAGTCATACTCTATTCTGAATTTCTCGCTAGT` / `ATAGTCATACTCTATTATGAATTTCTCGCTAGT` |
| tmp_1_11193_C_T | 0 | none | 1:11192–11194 | `TCT` / `TTT` | 15 | 1:11177–11209 (33 / 33 nt) |
| rs5413864115 | 2 | 1:11283–11285 | 1:11281–11286 | `TCAAAG` / `TCAAG` | 15 | 1:11266–11301 (36 / 35 nt) |
| tmp_1_11502_C_CGT | 0 | none | 1:11501–11503 | `CCA` / `CCGTA` | 15 | 1:11486–11518 (33 / 35 nt) |
| rs5413863413 | 13 | 1:13736–13751 | 1:13734–13752 | `ATTGGTGGTGGTGGTGGTA` / `ATTGGTGGTGGTGGTA` | 22 | 1:13712–13774 (63 / 60 nt) |

**Per genome** (the reference first, as a control, then each pan-genome genome; genome targets only):

1. **Anchors.** Every product of the set's two pairs from `amplicons ∪ unlikely` with orientation `LR` or `RL`, de-duplicated by `(region, start, end)`, at most `genotype_max_anchors` (50). `unlikely` products are included deliberately: the wrong-allele primer can be blocked while its amplicon is still a valid place to read the allele.
2. **Align.** Fetch `[min − 50, max + 50]` with `this.fetchWindow` (clamped), reverse-complement for `RL`, and align `R_s` semi-globally (query fully consumed, genome ends free, unit costs, traceback diagonal-first).
3. **Read the core.** Map the two end bases of the core span through the alignment.
   - **An end base that aligns to a gap** (the genome lacks it) moves outward, left for the start and right for the end, to the nearest reference position that aligns to a genome base.
   - If an end falls outside the aligned span, or no aligned base exists outward, the anchor is `missing` (`variant_not_covered`).
   - Otherwise `S` = the genome bases between the two mapped ends, gaps removed.

   Examples [V]:
   - At rs871475760, a genome lacking 1:11110 reads `TCG` (the end moves to 11111) → `other`.
   - A genome lacking 1:11108 reads `TCT` (the start moves to 11107, the same base) → `ref`. That is correct: a 1-bp deletion inside `TT` cannot be placed, and the observed core equals `C_ref`.
   - At rs5413864115, a genome lacking 1:11281 reads `TCAAAG` → `ref` for the same reason.
   - At rs5413863413, a genome lacking 1:13752 reads `ATTGGTGGTGGTGGTGGTG` → `other`.
   - An anchor that aligns over only part of `R_s` (43 % identity in the test case) is removed by the identity test of step 6 and is never read as an allele.
4. **Call — exact core match only.**
   - `S === C_ref` → `ref`; `S === C_alt` → `alt`; **anything else → `other`**, with `observed: S` always reported.
   - Edits outside the core are counted in `flank_edits` and are **confidence only**; they never change the call.
   - This is the rule that fixes the repeat failure mode. At rs5413864115 a genome that deleted two A's has core `TCAG`, one edit from `C_alt` and two from `C_ref`, so a lowest-distance rule names it `alt`; `AAAA` (`TCAAAAG`) is one edit from `C_ref` and would be named `ref`. Both are third alleles and both are `other` here. The same holds for the TGG repeat rs5413863413, where two deleted units score 3 edits from `C_alt` and 6 from `C_ref` [V].
5. **Merge anchors into copies.** Anchors whose **aligned variant coordinate** is the same (within the zone length) are one copy: `start`/`end` is the envelope and `anchors` the count; `identity`, `gap_compressed_identity`, `aligned_length` and `flank_edits` are those of the anchor alignment with the most columns (ties: the first anchor in product order). The minimum would report the shortest window, which the variant itself dominates (pi180348: 98.18 % over S1's 165 columns against 98.44 % over S2's 192). Merging by aligned coordinate — not by exact amplicon bounds or by 5′ site — is what keeps a locus from being counted twice when a forward and a reverse set of the same variant are checked together, or when two sets share the allele-specific primer but differ in the common primer.
6. **Ortholog filter.** A copy is orthologous when `ortholog === true` (gene-mode annotation) or when **gap-compressed identity** ≥ `genotype_ortholog_min_identity` (95) — each gap run counting as one edit — and every anchor's size is within ±20 % of its set's reference product. Raw identity would turn a true ortholog carrying one 30–40 bp indel into a paralog; both numbers are reported (`identity`, `gap_compressed_identity`).
7. **Genome allele.** No orthologous copy → step 8. All orthologous copies agreeing → that allele, `source: "amplicon"`, `observed` = the shared core. Disagreeing → `ambiguous`, `observed: null`.
8. **Megablast fallback.** Query `ref[zone.start − 200 − K .. zone.end + 200 + K]`; `blastn -task megablast -db <asm.blastdb.dna> -query - -dust no -soft_masking false -evalue 1e-20 -max_target_seqs 50 -max_hsps 10 -num_threads 1 -outfmt "6 sseqid sstart send sstrand pident length bitscore qstart qend qseq sseq"`, run through `this.spawnLines` in the worker's private cwd, **20 s timeout, no retry, at most `genotype_max_megablast` (30) fallbacks per job**. HSPs are kept when bitscore ≥ 0.9 × best, pident ≥ 95, query cover ≥ 0.8 and the query columns cover the core; `sseq` is already on the query strand [R], and the same exact-core rule applies. Beyond the budget a genome is `missing` with `reason: "fallback_budget"` and results warning `GENOTYPE_FALLBACK_BUDGET`; a spawn failure or timeout is `missing` with `fallback_failed` and `GENOTYPE_FALLBACK_FAILED`.
9. **Unavailable genomes** (`db_unavailable` or `error` in `pangenomeGenome`, `run.js:616-653`) → `allele: "unavailable"` with the status reason, and all predictions `unknown`. An exception inside the caller → `unavailable` with `reason: "call_failed"`, logged, plus results warning `GENOTYPE_FAILED`.

Verified behaviour [V, R]: for the designed set S1 of rs871475760 over all 119 other assemblies the caller gives 73 `ref`, 44 `alt` and 2 `missing` (is36143, pi525695), with 1 copy in 79 genomes, 2 in 37 and 4 in 1 (pi154987) — matching the research megablast table exactly; the inverted duplicate copy of pi180348 (1:15028–15132 − and 1:38342–38446 +) yields one `alt` genome with two copies; pi329250 is called `ref` although Ensembl genotypes that accession `A|A`; the tmp_1_11502 paralog at 1:72.9 Mb carrying GT (92 % identity) is filtered out; gap placement (megablast at 11051 vs dc-megablast at 11056) does not affect the call; and the 42 tied HSPs at tmp_1_11050_G_GA are `other`.

### 5.7 Amplification prediction (resolves finding 2)

For every set and genome, three primers are read, not one.

**Residual mismatches.** For an allele-specific primer, `residual_mm_pos` is `mm_pos` minus the positions recorded in `deliberate_mismatch_positions`; for the common primer it always equals `mm_pos`. Predictions read the residual, so a declared −2 mismatch is not mistaken for a genome difference.

**On-locus products.** A product of either pair is on-locus when it overlaps an orthologous copy (§5.6 step 6); every other product is off-locus.

**Allele-specific primers.** `ref_primer` is read from the REF pair's on-locus products and `alt_primer` from the ALT pair's, using that pair's product likelihood (`classify.classifyAmplicon`) and the allele-specific site's alignment:

| Condition (allele-specific primer, per product) | `status` |
|---|---|
| no product | `no_product` |
| product `unlikely` | `blocked` |
| `mm_pos` null (approximate alignment) | `unknown` |
| residual contains 1, and the primer's `in_shift_tract` is true (§4.12) | `uncertain` |
| residual contains 1 | `terminal_mismatch` |
| residual contains 2 or 3 (but not 1) | `weak` |
| otherwise | `match` |

**The common primer.** `common_primer` is read from **the common site's own alignment, never from a pair's product likelihood**. For every on-locus product of **either** pair, the shared primer's alignment in that product is classified alone: `classify.isIgnored`, `classify.isOverAmplifyingCap` or `classify.isBlocked` on that site → `blocked`; `mm_pos` null → `unknown`; `mm_pos` contains 1 → `terminal_mismatch`; 2 or 3 → `weak`; otherwise `match`. No deliberate-mismatch position is ever subtracted from the common primer. A deliberate mismatch on an allele-specific primer therefore can never be read as a common-primer mismatch: on a REF genome the ALT pair's product of an AS-PCR set is `unlikely` (`[2,1]`), yet the common site inside it is a perfect match, and the REF pair's product says the same [V].

For each primer the best status over its products and copies wins, in the order `match` > `weak` > `terminal_mismatch` > `uncertain` > `unknown` > `blocked` > `no_product`. Ranking `unknown` above `blocked`/`no_product` prevents one copy with no product from hiding an uncertain copy.

**Which product the reported fields come from.** `likelihood`, `mm_pos` and `residual_mm_pos` come from the one product that supplied the winning status. For `ref_primer` and `alt_primer` that is a product of the primer's own pair, with its `classify.classifyAmplicon` likelihood. For `common_primer` it is the on-locus product, of either pair, whose common site gave the winning status, with that product's likelihood. When several products tie on status, the one with the better likelihood (`likely` > `likely_weak` > `unlikely`) is reported, then the first in product order. So on pi180348 (§2.13) the KASP common primer reports `likely` from the ALT pair's product, and an AS-PCR set on a REF genome reports `common_primer {status: "match", likelihood: "likely"}` from the REF pair's product. With `no_product` the three fields are `null`.

**Set prediction** — the first matching row applies:

| # | Situation | `predicted` | `strength` | `reasons` |
|---|---|---|---|---|
| 1 | the genome is `unavailable` | `unknown` | — | |
| 2 | neither pair has an on-locus product (the locus is `missing`) | `none` | — | |
| 3 | common primer `terminal_mismatch` or `blocked` | `no_call` | — | `common_primer_3p_mismatch` |
| 4 | any of the three primers `unknown` or `uncertain` | `unknown` | — | `approx_alignment`, `shift_tract_uncertain` |
| 5 | ref `match`/`weak`, alt not | `ref` | `weak` if the ref or the common primer is `weak`, else `normal` | `common_primer_weak` when the common primer is `weak` |
| 6 | alt `match`/`weak`, ref not | `alt` | likewise | likewise |
| 7 | both `match`/`weak` | `both` | likewise | likewise, plus `third_allele` when the genome's allele is `other` |
| 8 | neither | `none` | — | |

**Off-locus products.** An off-locus product counts when `classify.countsAsAmplicon` accepts it, its allele-specific residual lacks 1, and its common site alone classifies as `match` or `weak`. Each counted product adds 1 to `off_locus_products` and the reason `ref_signal_off_locus` (REF pair) or `alt_signal_off_locus` (ALT pair). It changes the prediction **only by adding the dye that is missing**: a counted ALT-pair product turns `ref` into `both` and a counted REF-pair product turns `alt` into `both` — in the lab that genome lights both dyes; on a genome predicted `none` or `no_call` the off-locus dyes alone decide (`ref`, `alt` or `both`). A same-allele paralog (a REF-pair product on a genome predicted `ref`, or ALT on `alt`) leaves the prediction unchanged and only raises `off_locus_products` and its reason, because the lab still sees one dye. `unknown` is never changed.

**Verified** with `scratch-fix/predict.js`, the worktree's own `realign.js` and `classify.js` [V]:
- **Fake world.** A world built from the real reference was run for KASP S1, KASP S2 and the AS-PCR set of §2.10(d); each set gives:
  - REF genome → `ref`; ALT genome → `alt`;
  - third allele → `none` (allele `other`); missing locus → `none`;
  - two ALT copies → `alt`;
  - a blocked paralog → `ref` with `off_locus_products 0`;
  - an amplifying ALT paralog on a REF genome → `both` (`agrees: false`);
  - a REF paralog on a REF genome → `ref` with `off_locus_products 1`;
  - a common-primer 3′ SNP → `no_call`; a common-primer −2 SNP → `ref`, `strength: "weak"`.
- **Real assemblies.** On the recorded megablast alignments of the 11 research assemblies, S1, S2 and an AS-PCR set on S1's common primer agree 11/11. Reading the common primer from the ALT pair's product instead would give `no_call` on every REF assembly of the AS-PCR set.
- **Shift tracts** (`scratch-fix2/predict_tract.js`, same modules). Each genome has one locus: the reference ±350 bp with the repeat replaced. The sets are rs5413864115 S1 (forward) and S2 (reverse) of §2.10(b), and for rs5413863413 pair 0 of the §4.7 level-0 design runs, neighbours not screened:
  - forward: `CCATTGGTGGTGGTGGTG` (1:13732–13749) with common `CCCCGGTTATAACTGAGGATGG` (1:13813–13834), 103 bp;
  - reverse: `GTACTTCTACCACCACCACCAC` (1:13738–13759) with common `TCCAAGACAAGGTTCATGGAGT` (1:13678–13699), 82 bp.

  In all four sets `as_ref.in_shift_tract` is true and `as_alt.in_shift_tract` false. The results:
  - REF genome → `ref`, `agrees: true`, in all four sets (`alt_primer` `terminal_mismatch`). Under the anchor-inclusive wording both reverse sets were `unknown` here.
  - ALT genome → `unknown` with `shift_tract_uncertain`, `agrees: null` (`ref_primer` `uncertain`).
  - Third repeat lengths are all called `other`:
    - rs5413864115 two-A deletion (`TCAG`) → `alt`, weak, in both orientations;
    - `AAAA` (`TCAAAAG`) → `ref`;
    - rs5413863413 two TGG units deleted → `none`;
    - one unit inserted → `ref`.

    Only `none` agrees, so the others surface under "disagreements only".

**`agrees`** is `true` for (`ref`,`ref`), (`alt`,`alt`), (`other`,`none`) and (`missing`,`none`); `null` when the allele is `ambiguous`/`unavailable`, when `predicted` is `unknown` or `no_call`, or when any primer is `uncertain`; otherwise `false`.

The research's own cases map exactly [R]: REF `[]` → `match`; ALT `[1]` → `terminal_mismatch`; REF with a −3 mismatch `[3]` → after subtraction, residual `[]` → `match`; ALT with a −3 mismatch `[3,1]` → residual `[1]` → `terminal_mismatch` on a product that `classify` already calls `unlikely`, i.e. `blocked`.

### 5.8 Reference control, summaries and partial results

- **Reference control.** The reference genome runs through the same caller with `is_reference: true`. `control.status` is `pass` when its allele is `ref` and the set predicts `ref` on it; `warn` when the prediction is `weak` or an off-locus product appears; `fail` when the allele is not `ref`, the set does not predict `ref`, or the allele-specific pairs have off-target products that contradict the prediction. Any `fail` adds results warning `REFERENCE_CONTROL_FAILED {allele, sets}`.
- **Order.** `genotyping.genomes` is the reference first, then pan-genome genomes in request order, mirroring `pangenome.pairs[].genomes`.
- **Summaries.** `genotyping.summary` counts pan-genome genomes only and its parts add up to `genomes_total`, which equals `pangenome.pairs[i].summary.genomes_total`. Per set, `predicted_ref + predicted_alt + both + none + no_call + unknown = genomes_total` and `agree + disagree + not_comparable = genomes_total`.
- **Specificity-only jobs** have the reference entry only and zero summaries.
- **Partial results.** The first flush after `referenceStage` carries the reference entry; each pan-genome flush appends that genome and recomputes the summaries before `this.flush` (`run.js:699-705`), so `results.genotyping.genomes` grows like `pangenome.genomes`.
- **Size.** About 0.4 kB per copy and 0.2 kB per set × genome: 120 genomes × (2 copies + 5 sets) ≈ 220 kB uncompressed, well inside `max_result_bytes` [I].

### 5.9 Integration points

| File:line (worktree `primer-design`) | Change |
|---|---|
| `check/normalize.js:32` | `BODY_KEYS` gains `'genotyping'` |
| `check/normalize.js:83-143` | `validateShape` calls `genotype.validateShape(body.genotyping)` |
| `check/normalize.js:300`, `:310-323` | mode rule; `genotype.validateLinks(body.genotyping, pairsById)` |
| `check/normalize.js:404-411` | after the reference resolve: `await genotype.validateSets(...)` → `{variant, sets, orientations, mismatchPositions}` |
| `check/normalize.js:455-461` | `request.genotyping` (request-shaped) and `cost.estimate({..., genotyping: true})` |
| `check/cost.js:34-44`, `:104-130` | `DEFAULT_GENOTYPE_CPU_S_PER_GENOME = 0.2`, the extra term, conditional `breakdown.genotyping` |
| `check/index.js:12` | add `algorithmVersionFor(request)`; `ALGORITHM_VERSION` stays `'2'` |
| `jobs/index.js:164` | guarded `algorithmVersionFor` call (§5.3) |
| `check/run.js:211-261` | constructor: `this.geno = request.genotyping ? genotype.prepare(...) : null`, `this.genoEntries = new Map()` |
| `check/run.js:723-774` | `execute`: create `results.genotyping` and prepare before `referenceStage` |
| `check/run.js:556-585` | `referenceStage`: call the caller for the reference before `this.flush` (:584) |
| `check/run.js:616-653` | `pangenomeGenome`: all four return sites (624 unavailable, 634 error, 641 genome, 648 cDNA) record a genotype entry; `markOrthologs` is hoisted so copies can be tested |
| `check/run.js:655-721` | `pangenomeStage`: rebuild `results.genotyping` before `this.flush` (:705) |
| `check/blast.js:8-9`, `:32-53` | megablast args and the genotype outfmt plus its line parser |
| `jobs/worker.js` | no change; `ctx.spawnLines` and `fetchWindow` already cover the new calls |
---

## 6. File-by-file change list (gramene-swagger, worktree `gramene-swagger-genotyping`)

Every line number below was checked against the `primer-design` worktree while writing this spec [V]. Four corrections to note, because they change the work: **there is no startup validation of `primer3_core` anywhere**, so none is added for `ntthal` either (failures surface as `503` on first spawn); **there is no `PrimerCheckPair` definition** — the check's `pairs` schema is inline in `PrimerCheckRequest` (`swagger.yaml:2179-2198`) and stays inline, because `contract.test.js` reads it by path; **`buildRegionTemplate` and `requireRegionLength` are not exported** from `template.js` and must be added to its export list (`template.js:527-539`); and **`ALGORITHM_VERSION` exists twice** (`check/index.js:12` and the fallback `DEFAULT_ALGORITHM_VERSION` at `run.js:51`), so both stay `'2'`.

### 6.1 New modules (`api/helpers/primers/`)

| Path | Exports | Notes |
|---|---|---|
| `variation/client.js` | `createVariationClient({cfg, fetch, now, log})` → `{overlapChunk(species, region, chunkStart, {signal}), variation(species, id, {signal})}`; `VariationError` | limiter, single flight, breaker, LRU + negative cache, byte caps, `redirect: 'error'`, URL assertion (§3.3) |
| `variation/normalize.js` | `recordsToEntries`, `leftNormalize`, `shiftOf`, `tractOf`, `discriminating`, `zoneOf`, `altToRef`, `haplotypes`, `kindOf`, `minimalOf`, `labelOf`, `parseManual`, `submissionSequence`, `IUPAC` | pure; sequence access injected (§3.4–§3.5) |
| `variation/index.js` | `hasVariation`, `variationInfo`, `listVariants`, `lookupVariant`, `resolveDesignVariant`, `neighboursFor` | orchestration over genomes, assemblies, sequence, client, normalize |
| `thermo.js` | `createThermo({cfg, spawn, signal})` → `{hairpin, selfAny, selfEnd, cross, duplex, calls}` | `ntthal` only, pool of 4, memoized per `(salts, metric, s1, s2)`, clamps negatives, validates `^[ACGT]{1,60}$` (§4.10) |
| `genotyping/request.js` | `normalize(body, cfg)`, `ASSAY_DEFAULTS`, `PARAM_KEYS` | §2.7 rules 1–6; reuses `design.parseParams` / `checkParams` |
| `genotyping/presets.js` | `PRESETS`, `LADDER`, `FLOORS`, `levelParams(type, level, userParams)` | frozen tables (§4.6) |
| `genotyping/template.js` | `buildVariantTemplate(variant, resolved, deps)` | REF/ALT templates, `features`, zone in template coordinates, room check, exemption spans (§4.3) |
| `genotyping/guard.js` | `targetTag(zone, orientation, template, gap)`, `commonPrimerOk(set, zone, haplotypes)` | the `SEQUENCE_TARGET` guard and its post-filter (§4.4) |
| `genotyping/mismatch.js` | `TABLE_9_8_1`, `CLASSES`, `terminalPair`, `chooseMismatch`, `applyMismatch` | the literal 32-cell lookup (§4.11) |
| `genotyping/scoring.js` | `scorePair(haplotype, left, right, params, deps)` | one `check_primers` run, widened pair constraints, `problems` extraction, floors (§4.10, §4.6) |
| `genotyping/sets.js` | `candidatesFromPairs`, `deriveAlt`, `discrimination`, `issuesFor`, `scoreSet`, `rankSets`, `setKey`, `checkFragment`, `proposedCheckRequest` | pure |
| `genotyping/order.js` | `oligoName`, `orderRows`, `kaspMix` | pure (§4.17) |
| `genotyping/design.js` | `designGenotyping(body, deps)`, `runOrientation` | orchestration; reuses `design.createDeadline`, `design.buildRecord`, `primer3.run`, `boulder`, `repeat_mask`, `semaphore` |
| `check/genotype.js` | `GENOTYPING_VERSION`, `validateShape`, `validateLinks`, `validateSets`, `prepare`, `emptyResults`, `callGenome`, `semiGlobal`, `mergeCopies`, `gapCompressedIdentity`, `megablastCopies`, `predict`, `writeResults`, `unavailableEntry` | shared by the API (normalize) and the worker (run) |

Splitting the design into `guard.js`, `mismatch.js`, `scoring.js`, `sets.js` and `order.js` keeps each module independently testable; the largest is `sets.js`.

### 6.2 Changed files

| File | Lines | Change | Milestone (§10) |
|---|---|---|---|
| `api/controllers/primers.js` | 56-61, 102-120, 200-209 | add `variation` and `genotyping` to the modules map; three handlers — `listPrimerVariants`, `getPrimerVariant`, `designGenotypingPrimers` — the last copying the abort-on-close pattern of `designPrimers` (`AbortController` at 106, `onClose` 107-111, cleanup at 118); extend the exported object | M4 (the two variant handlers), M6 (`designGenotypingPrimers`) |
| `api/swagger/swagger.yaml` | after 1041; 2154-2198; 2208-2212; 2310-2341; 2395-2450; definitions section from 1063 | three new paths before `/{collection}` (1043); `PrimerCheckRequest.properties.genotyping`; `PrimerWarning.details`; genomes `variation` + `has_variation`; `PrimerCheckResults.properties.genotyping`; all new definitions (§2.5–§2.12) | M4: the two variant paths, the §2.6 definitions, `PrimerWarning.details`, the genomes fields; M6: the design path and the §2.7–§2.8 definitions; M7: `PrimerCheckRequest.properties.genotyping` and the §2.11 definitions; M8: `PrimerCheckResults.properties.genotyping` and the §2.12 definitions |
| `api/helpers/primers/config.js` | 19-99, 103-115, 163-185 | `DEFAULTS` gains `ntthal`, `variation`, `genotyping` and the `check.genotype_*` keys; `ENV_OVERRIDES` gains `NTTHAL`, `PRIMERS_VARIATION_URL` (new `url` type), `PRIMERS_VARIATION_ENABLED` | M4 `variation.*` and its two env keys; M5 `genotyping.*`; M6 `ntthal` / `NTTHAL`; M8 `check.genotype_*` |
| `config/default.yaml` | the `primers:` block | the same keys (§3.1) and the env list comment | as `config.js` |
| `api/helpers/primers/genomes.js` | 287-324, 332-355 | `genomeEntry` adds `has_variation`; `genomesResponse` adds `variation` | M4 |
| `api/helpers/primers/design.js` | 231 | export `parseParams` for `genotyping/request.js` (`checkParams` is already exported, `design.js:685`); no behaviour change | M5 |
| `api/helpers/primers/template.js` | 527-539 | export `buildRegionTemplate` and `requireRegionLength` (currently module-private) | M5 |
| `api/helpers/primers/primer3.js` | 18, 55-62, 73-98 | export `makeWorkDir`, `removeWorkDir` and `SPAWN_ENV` for `thermo.js`; no behaviour change | M6 |
| `api/helpers/primers/check/normalize.js` | 32, 83-143, 300, 310-323, 404-411, 455-461 | §5.9 | M7 |
| `api/helpers/primers/check/cost.js` | 34-44, 104-130 | the genotyping term and its conditional breakdown key | M7 |
| `api/helpers/primers/check/index.js` | 12 | add `algorithmVersionFor(request)`; `ALGORITHM_VERSION` stays `'2'` | M7 |
| `api/helpers/primers/jobs/index.js` | 164 | guarded `algorithmVersionFor` call (§5.3) | M7 |
| `api/helpers/primers/check/run.js` | 211-261, 556-585, 616-653, 655-721, 723-774 | §5.9 | M8 |
| `api/helpers/primers/check/blast.js` | 8-9, 32-53 | `GENOTYPE_OUTFMT`, `buildMegablastArgs({db, threads})`, `parseMegablastLine` | M8 |
| `api/helpers/primers/jobs/worker_main.js` | 60-72 | **bug fix, unrelated to this feature but required by the test plan:** install the SIGTERM/SIGINT handlers in `idle()` *before* the "idling" log line, so the log line is a true readiness signal (§7.1) | M1 |
| `api/helpers/primers/errors.js` | none | `normalizeStatus` (`errors.js:34-37`) already accepts 400–599 | — |
| `app.js` | none | the `/primers` middleware already covers the new paths | — |
| `docs/primer_design_api.md` | new section after the check section; the config table | "Genotyping primers (KASP / AS-PCR)": the three endpoints, the check extension, statuses, warnings and examples; config table gains `NTTHAL`, `PRIMERS_VARIATION_URL`, `PRIMERS_VARIATION_ENABLED` | M9 |
| `docs/genotyping_design_spec.md` | new | this document | M9 |
| `package.json` | none | Node 24 has global `fetch`; no new dependencies | — |

### 6.3 Test and fixture files

| Path | New or changed | Milestone (§10) |
|---|---|---|
| `test/primers/unit/compat_ids.test.js`, `compat_design_normalize.test.js` | new (compatibility, written first — §7.2) | M2 |
| `test/primers/unit/variation_normalize.test.js` | new | M3 |
| `test/primers/unit/variation_client.test.js` | new | M4 |
| `test/primers/unit/genotyping_request.test.js`, `genotyping_guard.test.js` | new | M5 |
| `test/primers/unit/genotyping_design.test.js` | new | M5 (records, explain strings, `pairs_returned`, `not_scored`); completed in M6 (sets, scores, budget, the §2.9 response) |
| `test/primers/unit/thermo.test.js`, `genotyping_mismatch.test.js`, `genotyping_sets.test.js` | new | M6 |
| `test/primers/unit/check_normalize_genotyping.test.js` | new | M7 |
| `test/primers/unit/genotype_caller.test.js` | new | M8 |
| `test/primers/unit/run_genotyping.test.js` | new | M2 (first half, C5), M8 |
| `test/primers/unit/design_genotyping_realdata.test.js` | new (`PRIMERS_REALDATA=1`) | M11 |
| `test/primers/unit/genotype_realdata.test.js` | new (`PRIMERS_REALDATA=1`) | M8 (the 11-assembly table), M11 (full panel, cross-validation) |
| `test/primers/unit/contract.test.js` | changed (§7.7) | M4, M6, M7, M9 |
| `test/primers/unit/controller_wrap.test.js` | changed (§7.3): exports and lazy loading (`:103-111`), `FEATURE_DISABLED` (`:239-261`), client disconnect (`:263-295`) | M4 (the two variant handlers), M6 (`designGenotypingPrimers`) |
| `test/primers/unit/job_id.test.js` | changed (§7.3): the `algorithmVersionFor` case beside `:181-191` | M7 |
| `test/primers/unit/cost.test.js` | changed: C4 (§7.2), then the genotyping term | M2, M7 |
| `test/primers/unit/genomes_response.test.js` | changed | M4 |
| `test/primers/unit/config.test.js` | changed: the new env keys | M4 (`PRIMERS_VARIATION_URL`, `PRIMERS_VARIATION_ENABLED`), M6 (`NTTHAL`) |
| `test/primers/unit/worker_memory.test.js` | changed (`:890`, §7.1) | M1 |
| `test/primers/integration/genotyping_endpoints.test.js`, `check_genotyping.test.js`, `variation_down.test.js` | new | M10 |
| `test/primers/fixtures/design/sorghum_bicolor_1_10500-12100.plus.txt` | new (1,601 nt, covers every example template) | M3 |
| `test/primers/fixtures/variation/overlap_1_10709-11509.json`, `overlap_1_10793-11593.json`, `overlap_1_10883-11685.json`, `overlap_1_11102-11903.json`, `variation_rs871475760.json`, `variation_tmp_1_11502_C_CGT.json`, `variation_rs5413864115.json`, `variation_not_found.json` | new (recorded live responses) | M4 |
| `test/primers/fixtures/primer3/genotyping/*.input.txt`, `*.output.txt`, `index.json` | new (recorded, §7.4) | M5 (design runs), M6 (scoring runs) |
| `test/primers/fixtures/thermo/genotyping.json` | new (recorded `ntthal` values keyed by argv) | M6 |
| `test/primers/fixtures/check_core/genotype/*` | new (megablast rows and assembly windows) | M8 |
| `test/primers/fixtures/contract/requests/genotyping-design-*.json`, `variants-*.json`, `check-genotyping-*.json` | new: hand-written provisionally by the API session in the milestone that adds each endpoint, then replaced by the output of the UI session's gramene-primers `npm run fixtures` (§7.7) | M4 (`variants-*`), M6 (`genotyping-design-*`), M7 (`check-genotyping-*`); replaced in M9 |
| `test/primers/integration/fake_ensembl.js` | new: the loopback fake Ensembl server of §7.6 — a data port serving the recorded bodies and a control port that switches modes (commands in §6.4) | M10 |
| `test/primers/fixtures/primer3/record_genotyping.js` | new: the recorder of §7.4 | M5, extended to scoring runs and thermo in M6 |
| `test/primers/fixtures/contract/compat_job_ids.json` | new (the 11 existing check job ids recorded on `primer-design` HEAD) | M2 |

### 6.4 Operations

**Worktree** (no npm install; `node_modules` is copied, and any later install uses `--ignore-scripts`):

```bash
cd /usr/local/gramene/subsites/sorghum/v11/gramene-swagger-primers
git worktree add --no-track -b genotyping ../gramene-swagger-genotyping primer-design
cd ../gramene-swagger-genotyping
cp -a ../gramene-swagger-primers/node_modules .
cp ../gramene-swagger-primers/package-lock.json .
source ~/.nvm/nvm.sh && nvm use 24.14.1
```

**Dev instance** (loopback only, its own site key and worker, never ports 50011 or 50111):

```bash
source ~/.nvm/nvm.sh && nvm use 24.14.1
W=/usr/local/gramene/subsites/sorghum/v11/gramene-swagger-genotyping; cd $W; mkdir -p .dev   # .dev/ is untracked; never commit it
HOST=127.0.0.1 PORT=50112 SWAGGER_HOST=localhost:50112 SWAGGER_SCHEMES=http \
  PRIMERS_SITE_KEY=sorghum_v11_geno PRIMERS_GLOBAL_MAX_JOBS=1 node app.js > .dev/api.log 2>&1 & echo $! > .dev/api.pid
PRIMERS_SITE_KEY=sorghum_v11_geno PRIMERS_GLOBAL_MAX_JOBS=1 node api/helpers/primers/jobs/worker_main.js > .dev/worker.log 2>&1 & echo $! > .dev/worker.pid
grep -H 'primers site_key=sorghum_v11_geno' .dev/api.log .dev/worker.log   # both must log it
```

**With the loopback fake Ensembl (M10, §7.6).** Start the fake server, then restart the API pointed at it. The worker never calls Ensembl and keeps running:

```bash
node test/primers/integration/fake_ensembl.js --port 50199 --control-port 50198 --mode normal > .dev/fake.log 2>&1 & echo $! > .dev/fake.pid
kill "$(cat .dev/api.pid)"; while kill -0 "$(cat .dev/api.pid)" 2>/dev/null; do sleep 0.2; done
HOST=127.0.0.1 PORT=50112 SWAGGER_HOST=localhost:50112 SWAGGER_SCHEMES=http \
  PRIMERS_SITE_KEY=sorghum_v11_geno PRIMERS_GLOBAL_MAX_JOBS=1 PRIMERS_VARIATION_URL=http://127.0.0.1:50199 \
  node app.js > .dev/api.log 2>&1 & echo $! > .dev/api.pid
curl -s -X POST http://127.0.0.1:50198/mode/slow   # modes: normal | slow | 500 | html404 | oversize | down
```

`down` closes the data port (connection refused), and any other mode reopens it; switching modes never restarts the fake server. Restart the API with the same three commands after the `down` group, so the next group starts with an empty breaker and negative cache. To return to live Ensembl, restart the API without `PRIMERS_VARIATION_URL`. Stop everything with `kill $(cat .dev/*.pid)`.

A separate site key is essential: the existing pm2 worker for port 50111 runs older code and would otherwise claim these jobs. Redis keys are separated by site key (`redis_store.js:40-43`), and `PRIMERS_GLOBAL_MAX_JOBS=1` is only more conservative than the host cap.

**Repointing the pm2 dev apps** — only with the user's explicit go-ahead, after §7.8 passes:

1. Record the current state first: `pm2 describe sorghum_primers_dev_api` and `sorghum_primers_dev_worker` (pm_id 66 and 65), and their `/proc/<pid>/environ`. Today they run with `PRIMERS_SITE_KEY=sorghum_v11_dev`, `PRIMERS_GLOBAL_MAX_JOBS=1`, `PORT=50111`, `SWAGGER_HOST=localhost:50111` [V].
2. Stop and recreate the **worker first**, then the API, from `gramene-swagger-genotyping` with identical env and the Node 24.14.1 interpreter. Worker-before-API matters: a new API with an old worker would produce genotyping jobs whose results lack the block and then serve them from cache for `ttl_done_s` (86,400 s).
3. Confirm both log `primers site_key=sorghum_v11_dev` and that `https://data.sorghumbase.org/sorghum_v11a/primers/genomes` answers.
4. Keep the rollback command ready (the same two `pm2 start` invocations with `cwd` pointing back at `gramene-swagger-primers`).
5. **Do not run `pm2 save`.**

**Rollback / feature flags.** `PRIMERS_VARIATION_ENABLED=0` disables the variants endpoints (`503 FEATURE_DISABLED`) while manual genotyping design keeps working; `PRIMERS_ENABLED=0` disables everything, as today.

**Cleanup after testing:** `redis-cli -p 6380 -n 1 --scan --pattern 'primers:sorghum_v11_geno:*' | xargs -r redis-cli -p 6380 -n 1 DEL`.
---

## 7. Test plan

### 7.1 Baseline, and the flaky test (resolves finding 12)

**Today's baseline, measured twice while writing this spec** [V]: `node --test "test/primers/unit/**/*.test.js"` → **503 tests, 45 suites, 487 pass, 0 fail, 16 skipped** (all 16 need `PRIMERS_REALDATA=1`), ≈ 8.2 s. Both runs were green.

Note that the unit suite is **not** fully offline: `redis_store.test.js` writes to `redis://localhost:6380/1` (keys `primers:devtest_jobs*`) whenever that Redis answers, and it did in both runs [V]. New tests must not add such a dependency.

**The known-flaky test is a real race in the product code, not test noise.** `worker_memory.test.js:874` ("`node worker_main.js` from cwd `/` … exits 0 on SIGTERM") waits for the `idling` line on stderr (`:885`), then sends SIGTERM (`:888`) and asserts `ex.code === 0` (`:890`). In `worker_main.js` `idle()`, the log line is written at `:64` but the `process.once('SIGTERM', …)` handler is only installed at `:71`; the measured gap between them is 148–210 µs (median 157 µs) [V]. A SIGTERM that lands in that gap kills the child by default action, and the exit event then reports `code null, signal 'SIGTERM'` — exactly the observed failure. Injecting a 30 ms stall in the gap reproduces it 10/10; signalling on the first byte reproduces it 10/10; with the handlers installed first, it passes 10/10 [V]. A quiet host hides it; a loaded full-suite run can expose it.

**Fix (milestone M1):** reorder `idle()` in `api/helpers/primers/jobs/worker_main.js` so both signal handlers are installed *before* the log line, making the line a true readiness signal. No other test depends on the ordering (`redis_store.test.js:648,692` wait for a running job or the worker lock, which only exist after `worker.start()` at `:137`, i.e. after the handlers at `:134`). The assertion at `:890` is tightened to `should(ex).match({code: 0, signal: null})` so a future failure shows the signal.

**Gate.** The verification checklist compares against the recorded baseline — "no failures, and the pass/skip counts equal the baseline plus the new tests" — rather than "503 pass". If a failure appears, rerun the single file alone before calling it a regression.

### 7.2 Compatibility tests, written first

These are written and green **before** any feature code, with the reference values recorded from `primer-design` HEAD:

| # | Test | Pins |
|---|---|---|
| C1 | `compat_ids.test.js` | normalize each of the 11 existing `check-*.json` contract fixtures with the stub catalog and assemblies, compute `jobs.jobId(norm.request, norm.dbs, algo)`, and compare with `fixtures/contract/compat_job_ids.json`, recorded on `primer-design` HEAD. Every id must be unchanged after the feature lands. |
| C2 | `compat_design_normalize.test.js` | `design.normalize` on each of the 11 existing `design-*.json` fixtures deep-equals a snapshot recorded on `primer-design` HEAD |
| C3 | `boulder.test.js` (extended) | `extractPairs` on the existing recorded Primer3 outputs deep-equals today's result, with no new keys |
| C4 | `cost.test.js` | every existing estimate is unchanged and `breakdown` still has exactly its four keys when `genotyping` is not requested (`cost.test.js:52`) |
| C5 | `run_genotyping.test.js` (first half) | a non-genotyping `CheckRun` on the fake world produces exactly today's `Object.keys(results)` — no `genotyping` key |
| C6 | full suite | the §7.1 baseline holds |

### 7.3 New unit tests, with concrete expected values

All values below are verified in this document and are pinned as assertions.

| File | Case | Expected |
|---|---|---|
| `variation_normalize.test.js` | rs871475760 record + FASTA fixture | `key 1:11109:C:A`, `kind snv`, `shift 0`, `tractOf` null, zone `[11109,11109]`, both discriminating `{11109,C,A,11109}`, `label "1:11109 C/A"`, `ref_verified true` |
| | tmp_1_11193_C_T | `key 1:11193:C:T`, `ems true` (source `EMS_PMID38100514_Jiao`) |
| | rs5413864115 (`A/-` at 11283) | `key 1:11282:CA:C`, `kind deletion`, `shift 2`, `tractOf` `[11283,11285]`, zone `[11282,11286]`, forward `{11285,A,G,11286}`, reverse `{11283,A,C,11282}`, `minimal {11283,11283,A,-}` |
| | tmp_1_11502_C_CGT (start 11503, end 11502) **and** rs5413863549 (start = end = 11503) | both → `1:11502:C:CGT`, one entry, ids ordered with the requested id first; forward `{11503,A,G,null}`, reverse `{11502,C,T,null}`; zone `[11502,11503]`; `tractOf` null |
| | rs5413863413 `TGG/-` | `key 1:13735:TTGG:T`, `shift 13`, `tractOf` `[13736,13751]`, forward `{13749,G,A,13752}`, reverse `{13738,G,T,13735}`, zone `[13735,13752]` |
| | rs5413863494 `C/T/G`; rs5413863901 `C/T/*` | two entries with `other_alts`; one entry plus issue `STAR_ALLELE` |
| | manual `{11283,"A","-"}` and `{11282,"CA","C"}` | the same entry |
| | manual ref `A` at 1:11109 | `REF_MISMATCH {given:"A", genome:"C"}` |
| | `submissionSequence` for the four example variants | the exact strings of §4.17 |
| | synonyms `[".", "tmp_1_11109_C_A"]` | `["tmp_1_11109_C_A"]` |
| | window 1:11180–11290 | the four entries of §2.3, deep-equal |
| `variation_client.test.js` (fake `fetch`) | 200 overlap array | validated records; a malformed record is dropped and counted |
| | an id containing `,` and `*`, and a 255-character id | accepted, `encodeURIComponent`-ed, URL starts with the configured base |
| | a record id with `,`/`*`/200 chars | **kept** (lenient record validation); `VARIATION_RECORDS_SKIPPED` is not raised |
| | 400 `{"error":"rs0000000001 not found for sorghum_bicolor"}` | `404 UNKNOWN_VARIANT`, cached 5 min, breaker untouched |
| | 400 other text, 404 HTML, non-JSON, body over the cap | `503 {reason:"invalid_response"}`, not cached, breaker +1 |
| | 503 / 429 / timeout / `ECONNREFUSED` | `503` with the matching `reason`, negatively cached 30 s |
| | 3 transport failures in 60 s, then a 4th call | `503 {reason:"breaker_open"}` with **no** fetch |
| | 5 concurrent calls, `max_concurrent 4`, slow fetch | the 5th waits, then `503 {reason:"queue_full", retry_after_s:5}` |
| | two concurrent calls for the same key | a single fetch (single flight) |
| | a 25 kb window | exactly 3 chunk fetches; a second overlapping window reuses them |
| `thermo.test.js` | argv builder | contains `-mv 50 -dv 1.5 -n 0.6 -d 50 -t 37 -r` |
| | parsing | `"No secondary structure…"` → 0; `-11.46` → 0; `42.616970` → 42.62 |
| | memoization; invalid input `ACGN…` | one spawn per distinct key; throws before spawning |
| | gated on the binaries | FAM+REF hairpin 42.62, HEX+ALT 51.50, tails alone 0, FAM-REF × HEX-ALT `ANY` 6.89, duplex `ATCTTTGACTAGCGAGAAATTCGG` × plus footprint 54.69 |
| `genotyping_mismatch.test.js` | all 32 cells of Table 9.8.1 | the table of §4.11 |
| | rs871475760 −2, both orientations, both alleles | `…TTCGG`, `…TTCGT`, `…CTCTATCC`, `…CTCTATCA` |
| | tmp_1_11193_C_T −2 | reverse REF `…ACTCTCG` (A→C) and ALT `…ACTCTTA` (A→T); forward both T→G |
| | tmp_1_11193_C_T −3 | reverse both T→G; forward both C→A; `source` mentions the −3 extrapolation |
| | indel where the −k base differs between haplotypes | no mismatch, warning `MISMATCH_NOT_APPLICABLE` |
| `genotyping_guard.test.js` | forward/reverse tag construction for the five example variants | `SEQUENCE_TARGET` = `402,10`, `391,10`, `405,10`, `390,10`, and for rs5413863413 forward `{t(13753)},10` |
| | post-filter | a common primer at 1:11108–11131 is rejected (`overlap`) for rs871475760; a common primer at 13751–13780 is rejected for rs5413863413, because it covers the ALT 3′ base at 13752 |
| `genotyping_sets.test.js` | `deriveAlt` for the five rows of §4.9 | exact ALT sequences, blocks, `inserted_bases`, product sizes 65/65, 92/92, 88/87, 129/128, 72/74 (with the §2.9–§2.10 common primers) |
| | `setKey` | `f9df650ad116`, `1accc54c262d`, `f5a5c6f2ebcf`, `c791a4956fd3`, `a5277232d8ab`, `cb3ef66afd37`, `53942cb55348`, `7f9af6b1c938` |
| | discrimination via `realign` + `classify` | own `[]`/`likely`, other `[1]`/`likely_weak`; with a −2 mismatch `[2]`/`likely` and `[2,1]`/`unlikely` |
| | floors | an ALT primer at 51.9 °C is dropped (`below_floor`); at 54.89 °C with 19.35 % GC it is kept with `ALT_PRIMER_SUBOPTIMAL` |
| | issues, quality and score (§4.16) | rs871475760 kasp S1 `f9df650ad116` 9.18 `usable`, S2 `1accc54c262d` 19.63 `poor`; tmp_1_11193_C_T S1 22.51 before S2 24.84; rs5413864115 S1 15.48, S2 23.95; tmp_1_11502_C_CGT S1 9.56 (`AS_TM_IMBALANCE` not counted as an issue); as_pcr S1 7.92; `WEAK_TERMINAL_CLASS` counts 0; the eight scored forward rs871475760 candidates rank 19.63 first and Primer3's pair 0 (22.30) last |
| | exact arithmetic and rounding (§2.1) | tmp_1_11193_C_T S1 `common_minus_as` from Primer3's `57.102` and `56.427` is exactly `0.675` → `0.68`, never `0.67`; scores keep their exact sums (`9.179713` → `9.18`, `23.953151` → `23.95`); the rounding helper gives `0.125` → `0.13`, `−0.125` → `−0.13` and `2.675` → `2.68` (binary floating point gives 2.67); GC of the 31-nt S2 ALT primer is 6/31 → `19.35`; `primer3_penalty` `7.179713` → `7.1797` |
| | score rules S1 and S6 (§4.16) | tmp_1_11193_C_T S1 has exactly one allele-specific and one common `NEIGHBOUR_IN_PRIMER` (22.51, not 23.51); rs871475760 AS-PCR S1 has two `ALT_PRIMER_SUBOPTIMAL`, both from the mismatch runs, and none from the matched ALT primer's run (7.92, not 8.92, which would promote `f2adeabf4265`); rs871475760 S2's penalty term is its own pair 2's 14.634369, not pair 0's 14.376148 |
| | shift tract (§3.5, §4.12) | rs5413864115 S1 and S2: `as_ref.in_shift_tract` true, `as_alt.in_shift_tract` false, one `SHIFT_TRACT_DISCRIMINATION` each; S2 scores 23.95 (an anchor-inclusive tract would give 24.95) |
| | `proposedCheckRequest` | stops at 5 sets, 10 pairs, 13 primers; deep-equals §2.9 `check.request` |
| | `orderRows` | names `rs871475760_S1_REF_FAM`, `…_ALT_HEX`, `…_COM`; label sanitising `1:11109:C:A` → `1_11109_C_A` |
| `genotyping_request.test.js` | `{id, position}`; `ref "-" alt "-"`; unknown keys | `INVALID_VARIANT id_or_manual`; `INVALID_VARIANT alleles`; `INVALID_REQUEST` |
| | defaults | kasp → tails `ref_fam_alt_hex`, mismatch `none`; as_pcr → `none` / `auto` |
| | pinned params | `params.min_tm: 58` never relaxed; `settings.pinned` contains `min_tm` |
| | product minimum | `max_size 30` → `[[61,120]]`; user `[[40,90]]` with `max_size 25` → `[[51,90]]` |
| `genotyping_design.test.js` (stubbed sequence, variation client, recorded Primer3 and thermo) | rs871475760 kasp, 2 sets | deep-equals the §2.9 response except warning wording |
| | tmp_1_11193_C_T | S1 forward `f5a5c6f2ebcf` L1, S2 reverse `c791a4956fd3` L1, `EMS_TARGET`, two `RELAXED_CONSTRAINTS`, the four explain strings of §2.10(a) |
| | rs5413864115 manual Ensembl style | S1 `a5277232d8ab` forward L1, S2 `cb3ef66afd37` reverse L1, `SHIFTABLE_INDEL`, `DENSE_NEIGHBOURS 4` |
| | tmp_1_11502_C_CGT | forward `blocked` by rs5413863234 at distance 4; S1 `53942cb55348`, `AS_TM_IMBALANCE 1.16`, `NEIGHBOUR_IN_PRIMER` for the common primer (rs5413863452 at 6), `DUPLICATE_VARIANT_IDS`, `DENSE_NEIGHBOURS 5` |
| | as_pcr rs871475760 | S1 `7f9af6b1c938` with the mismatch block of §2.10(d) |
| | scoring cap and budget (§4.8) | **Default caps** (rs871475760 kasp): forward attempt `pairs_returned 20` = `not_scored 12` + `sets 8`; reverse `common_neighbour_3p 6` + `not_scored 6` + `sets 8`; 18 Primer3 runs, 123 `ntthal` calls, no `DESIGN_BUDGET_EXHAUSTED`. **Same request with `genotyping.max_primer3_runs: 12`:** `200` and no exception; forward unchanged; reverse `sets_found 2`, attempt `not_scored 18` + `sets 2`; warning details deep-equal `DESIGN_BUDGET_EXHAUSTED {orientation: "reverse", primer3_runs: 12, thermo_calls: 87, max_primer3_runs: 12, max_thermo_calls: 272, not_scored: 18, sets_returned: 2}`; the returned sets are still `f9df650ad116` and `1accc54c262d` (reproduced by `scratch-fix/score_sets.js` with `CAP_PRIMER3_RUNS=12` and by `scratch-fix2/score2.js` with `CAP_P3=12`) |
| | `template_only` | no Primer3 and no thermo spawn; `sets: []`, `check: null`; **the semaphore is not acquired** |
| | Ensembl 503 during a manual design | 200 with `NEIGHBOURS_UNAVAILABLE`, `ids: []` |
| | Ensembl 503 during a design by id | 503, and **no** semaphore acquisition (asserted through a stub semaphore) |
| | deadline | a Primer3 stub that never resolves → `504 DEADLINE_EXCEEDED` |
| | mask stub covering the AS window | the Primer3 template has the AS window unmasked; `VARIANT_IN_REPEAT` |
| `genotype_caller.test.js` | cores and haplotypes (`prepare`) | the five rows of the §5.6 table: rs871475760 `TCT`/`TAT` with 33-nt haplotypes 1:11093–11125; tmp_1_11193_C_T `TCT`/`TTT`; rs5413864115 `TCAAAG`/`TCAAG`, haplotypes 1:11266–11301 (36/35 nt); tmp_1_11502_C_CGT `CCA`/`CCGTA`, K 15; rs5413863413 `ATTGGTGGTGGTGGTGGTA`/`ATTGGTGGTGGTGGTA`, K 22 |
| | core end on a gap (§5.6 step 3) | rs871475760 genome lacking 1:11110 → `other`, observed `TCG`; lacking 1:11108 → `ref`, observed `TCT`; rs5413864115 lacking 1:11281 → `ref`; rs5413863413 lacking 1:13752 → `other`; no aligned base outward → `missing` |
| | third lengths of an insertion | tmp_1_11502_C_CGT with `GTGT` inserted → `other` (`CCGTGTA`); `G` only → `other` (`CCGA`) |
| | truncated anchor | an anchor that covers the core but only 43 % of `R_s` is excluded by the step-6 identity test, never read as an allele |
| | exact-core rule at rs5413864115 | `TCAAAG` → ref; `TCAAG` → alt; `TCAG` (two A's deleted) → **other**; `TCAAAAG` → **other**; each with `observed` reported |
| | exact-core rule at rs5413863413 | one unit deleted → alt; two units deleted → other; one unit inserted → other |
| | shift-tract predictions (§4.12, §5.7): the four sets and synthetic genomes of the §5.7 shift-tract bullet | `in_shift_tract` is `as_ref` true, `as_alt` false for rs5413864115 S1 and S2 and for both rs5413863413 sets. REF genome → `ref`, `agrees: true`, in all four. ALT genome → `unknown` with `shift_tract_uncertain`, `agrees: null`. Two-A deletion → allele `other`, predicted `alt` (weak); `AAAA` → `other`, `ref`; two TGG units deleted → `other`, `none`; one unit inserted → `other`, `ref` |
| | primer-call fields (§5.7) | on the §2.13 genomes every `ref_primer`, `alt_primer` and `common_primer` object equals the example, including pi180348's `common_primer.likelihood` `likely` from the ALT pair's product |
| | SNP third allele | genome `T` at a C/A site → `other` |
| | `semiGlobal` | query fully consumed, free ends, core mapped through a 1 bp indel in the flank |
| | copy merging | S1 anchor 1:31523–31587 and S2 anchor 1:31536–31627 merge into one copy 31523–31627 with `anchors 2`, `aligned_length 192` and `identity 100` from the S2 alignment, the longest; the two pi180348 copies report 192 columns / 98.44 % and 193 columns / 98.45 %; a forward and a reverse set of the same locus also merge |
| | ortholog test | raw identity 94.9 % with one 35 bp indel but gap-compressed 99 % → ortholog; size +21 % → paralog; `ortholog: true` wins |
| | megablast parsing | the two pi180348 rows of the research `megablast.tsv` → `alt`, `alt`; budget 30 exhausted → `missing`/`fallback_budget` |
| `check_normalize_genotyping.test.js` | the §2.11 body | normalized block equals the body (uppercased, left-aligned) with **no** added keys; `orientation` is not stored |
| | `mode: "transcript"` | `INVALID_REQUEST {field:"genotyping", reason:"mode"}` |
| | swapped alleles; a primer shifted by 1 nt; 2 edits; a common primer inside the zone; missing/differing `expected` | `GENOTYPING_SET_INVALID` with reasons `alleles_swapped`, `not_at_variant`, `too_many_edits`, `common_in_zone`, `expected_required`, `expected_differs` |
| | the as_pcr mismatch primers of §2.10(d) | accepted; `deliberate_mismatch_positions [2,2]` |
| | job ids | the §2.11 body hashes with algo `'2+g1'`; the same body without `genotyping` keeps its pre-change id; the stubbed check module without `algorithmVersionFor` still works (`job_id.test.js:183`) |
| `job_id.test.js` (M7) | `jobs.submit` beside `:181-191` | with a check module that has `algorithmVersionFor`, a genotyping body hashes with algo `'2+g1'` and gets a different id from the same body without `genotyping`; the existing test's stub (only `ALGORITHM_VERSION` and `normalize`) still gives today's id |
| `controller_wrap.test.js` (M4, M6) | exports `:103-111`, `FEATURE_DISABLED` `:239-261`, client disconnect `:263-295` | the controller exports `listPrimerVariants` and `getPrimerVariant` (M4) and `designGenotypingPrimers` (M6) without loading `variation/index.js` or `genotyping/design.js`; all three answer `503 FEATURE_DISABLED` when primers are disabled; a disconnect aborts the genotyping design signal as it does for `designPrimers` (M6) |
| `run_genotyping.test.js` (the `verdicts` fake world, real bgzip re-alignment) | synthetic species run with **three sets** — a KASP reverse set, a KASP forward set and the AS-PCR set of §2.10(d) (mismatch −2 on both AS primers). Genomes: REF, ALT, third allele, no locus, duplicated ALT locus, an 85 % paralog whose common site is blocked, an amplifying ALT paralog next to a REF locus, an amplifying REF paralog next to a REF locus, a private SNP at the common primer's 3′ base, and one at its −2 base | calls `ref`, `alt`, `other`, `missing`, `alt` (2 copies), `ref` (paralog excluded), `ref`, `ref`, `ref`, `ref`. Predictions, identical for each of the three sets: `ref`, `alt`, `none`, `none`, `alt`, `ref` (`off_locus_products 0`), `both` (`alt_signal_off_locus`, `agrees: false`), `ref` (`off_locus_products 1`, `ref_signal_off_locus`), **`no_call`** (`common_primer_3p_mismatch`, `agrees: null`), `ref` with `strength: "weak"` (`common_primer_weak`). For the AS-PCR set on the REF genome `alt_primer` is `blocked` while `common_primer` is `match` (§5.7). Summaries add up; the partial flush after the reference contains only the reference |

### 7.4 Recorded fixtures

- **Primer3.** `test/primers/fixtures/primer3/record_genotyping.js`, run manually on squam, writes one `.input.txt` / `.output.txt` pair per run of the §7.3 design cases plus the `check_primers` scoring runs, and an `index.json` keyed by the **sha256 of the exact serialized input**. The stub `primer3.run` serves output only for a byte-identical input and **fails the test on any unrecorded input**, so a fixture can never silently drift from the code.
- **Thermo.** The same recorder writes `fixtures/thermo/genotyping.json`, keyed by the full argv. The stub behaves the same way.
- **Equivalence tests** (gated on the binaries existing): replay 50 recorded thermo entries against the real `ntthal`; assert that `check_primers` on the REF haplotype reproduces the design run for every recorded set (Tm ±0.01, `*_TH` ±0.01, end stability ±0.01).
- **Variation.** Live Ensembl responses for the four example windows and the three example ids, recorded once, plus a `not_found` body.

### 7.5 Real-data tests (`PRIMERS_REALDATA=1`, serial)

| Test | Expectation |
|---|---|
| `design_genotyping_realdata.test.js` | real FASTA, `primer3_core`, `ntthal`, variation client stubbed from the recorded fixtures: the eight set keys of §7.3 reproduce; wall time per design < 3 s |
| `genotype_realdata.test.js`, 11 research assemblies | S1 copies and calls: bicolorv5 1:31523–31587 ref; austrcf317961 ref; pi565121 ref; pi655972 ref; pi180348 1:15028–15132 (−) and 1:38342–38446 (+) alt; pi276837 alt ×2; pi510757 alt; pi656027 alt ×2; pi329250 ref; rio alt ×2; s3691 ref. All 11 predictions agree for S1, S2 and an AS-PCR set built on S1's common primer (`scratch-fix/predict.js` part B) |
| `genotype_realdata.test.js`, full panel (`PRIMERS_REALDATA_SLOW=1`) | S1 over 119 assemblies: 73 ref / 44 alt / 2 missing (is36143, pi525695); copies 1:79, 2:37, 4:1 (pi154987); 119/119 agree. pi536008, whose amplicon-window identity is 78.2 %, is expected to reach the megablast fallback and be called `alt` **[U]** |
| tmp_1_11502_C_CGT caller | every orthologous copy `ref`; the pi180348 paralog carrying GT at 1:72.9 Mb (92 %) is excluded, `paralog_copies: 1` |
| Ensembl cross-validation (test-only, never product behaviour) | for the 52 sample-matched assemblies, calls agree with 29 ref + 15 alt Ensembl genotypes; the documented disagreements are exactly pi329250, pi656031, pi533991 (Ensembl `A\|A`, assembly C) and grif16309, pi154987, pi656029 (Ensembl `C\|C`, assembly A) |

### 7.6 Integration tests on the 50112 dev instance

`PRIMERS_IT_BASE=http://127.0.0.1:50112/sorghum_v11`, serial, one job at a time.

A **loopback fake Ensembl server** (`test/primers/integration/fake_ensembl.js`, started and switched, with the matching API restart, by the commands of §6.4) is selected through `PRIMERS_VARIATION_URL=http://127.0.0.1:<port>` (http is accepted only for 127.0.0.1) with modes: `normal` (serves the recorded bodies), `slow` (9 s), `500`, `html404` (an Apache-style page), `oversize` (a body above the cap) and `down` (connection refused). Live-Ensembl smoke tests are gated behind `PRIMERS_IT_ENSEMBL=1`.

| Test | Expectation |
|---|---|
| variants list / lookup | the four keys of §2.3; both aliases plus `DUPLICATE_VARIANT_IDS` for tmp_1_11502_C_CGT |
| id syntax | `GET /primers/variants/tmp_1_13549_TTA_T%2C%2A?system_name=…` resolves or answers 404 `UNKNOWN_VARIANT`, but **never** 400 `PATTERN`; a 254-character id likewise |
| errors | unknown id → 404; `sorghum_rio` → 422 `NO_VARIATION_DATA`; a 60,000 bp window → 400 `VARIANT_WINDOW_TOO_LONG`; `/primers/variants/a%2Fb` → 400 validator |
| design | the §2.9 request → `sets[0].key f9df650ad116`, `sets[1].key 1accc54c262d`, `check.request` equal to §2.9 |
| design errors | §2.10(e) → 400 `REF_MISMATCH`; `rs5413863494` without `alt` → 400 `ALT_REQUIRED` |
| no-store / CORS | `cache-control: no-store` and `access-control-allow-origin: *` on all three new endpoints |
| **slot starvation** | fake Ensembl in `slow` mode: 4 concurrent genotyping designs **and** 4 concurrent ordinary `/primers/design` requests; **no** PCR design may receive `503 BUSY`, and each genotyping design either succeeds or fails with a variation error — this is the regression test for finding 6 |
| Ensembl down | `down` mode: the list answers 503 within 10 s, the next call answers immediately from the breaker, and a manual design still returns 200 with `NEIGHBOURS_UNAVAILABLE` |
| genotyping check, specificity only | 202 → done; `results.genotyping.genomes` holds the reference only, allele `ref`, both sets predict `ref`, `control.status: "pass"` |
| genotyping check, 3 genomes | `summary {ref:2, alt:1}`; pi180348 has two copies; 3/3 agree for both sets; `estimate.cpu_s 95` |
| job id stability | re-POST an existing `check-gene-specificity.json` body → the id pinned in C1 |
| cost | a 14-primer genotyping check over all genomes → 422 `JOB_TOO_LARGE` |

### 7.7 Contract tests and fixtures

`contract.test.js` changes, all against verified anchors:

| Anchor | Change | Milestone (§10) |
|---|---|---|
| `:92-116` | the operations list gains `listPrimerVariants` and `getPrimerVariant` (M4) and `designGenotypingPrimers` (M6); the `consumes` assertion at `:107` passes because both new GETs declare `consumes: [application/json]`, as the existing GETs do | M4, M6 |
| `:118-132` | the strict-definition list gains `PrimerGenotypingRequest`, `PrimerVariantInput`, `PrimerGenotypingAssay`, `PrimerGenotypingParams` (M6) and `PrimerCheckGenotyping`, `PrimerCheckGenotypingVariant`, `PrimerCheckGenotypingSet` (M7) | M6, M7 |
| `:136-187` | documented examples: the §2.9 and §2.10 design request bodies, plus invalid design bodies expected to fail, validated against `POST /primers/genotyping/design` (M6); the §2.11 check body, plus invalid `genotyping` blocks, validated against `PrimerCheckRequest` (M7) | M6, M7 |
| `:317-364` | sync tests, design half: `PrimerGenotypingParams` keys ⊆ `design.PARAM_SPECS` with equal bounds (except `product_size_ranges`); every `PrimerGenotypingRequest` property is known to `genotyping/request.normalize`; the assay enums equal the module's tables | M6 |
| `:366-387` | sync test, check half: `check/normalize` `BODY_KEYS` contains `genotyping`, and the §2.11 body passes `checkNormalize.validateShape` beside the existing full-body assertion. It needs the M7 `normalize.js` change, so it is never part of M6 | M7 |
| `:401-419` | `fixtureCase` treats `variants-*.json` as `{method, path, query}` wrappers (M4), routes `genotyping-design-*.json` to `POST /primers/genotyping/design` (M6) and keeps `check-genotyping-*` as check bodies (M7) | M4, M6, M7 |
| `:447-479` | the handler-layer pass routes genotyping design bodies through `genotyping/request.normalize` (M6) and genotyping check bodies through `check/normalize` up to the catalog stub (M7) | M6, M7 |
| new test | **blocking** example validation: every JSON example in `docs/genotyping_design_spec.md` is validated against the new definitions with an **`x-nullable`-aware walker**, because sway 1.0.0 ignores `x-nullable` when validating responses and would reject every legitimate null (including today's `PrimerCheckJob.error`) | M9 |
| new test | no `UNUSED_DEFINITION` warning for any `Primer*` definition (sway already reports these at `:84-90`) | M4, then re-run in M6, M7 and M8 |

**Ownership.** gramene-primers belongs to the UI session, and the API session never edits it:
- Until the generator changes below land, each API milestone that adds an endpoint (M4, M6, M7) hand-writes that endpoint's request fixtures in `test/primers/fixtures/contract/requests/` from the §2.3–§2.11 examples, under their final file names.
- In M9 the API session hands the generator changes to the UI session.
- The UI session's `npm run fixtures` output then replaces the provisional files, and the diff is reviewed.

The generator changes: `npm run fixtures` gains `genotyping-design-{id,manual-deletion,aspcr,template-only,all-fields}.json`, `variants-{list,lookup}.json` and `check-genotyping-{two-sets,gene-mode}.json`; the stale-file regex becomes `^(design|check|genotyping|variants)-`; the manifest learns the new definition and the two GET operationIds; the generator picks the preset from `mode`/`assay.type` instead of always validating against `pcr`; and the "four design modes" assertion (`fixtures.gen.test.ts:177-178`) is scoped to `design-*` files. That assertion uses `expect.arrayContaining`, so it would not have failed anyway [V] — it is scoped for clarity, not necessity.

### 7.8 Verification checklist

```bash
W=/usr/local/gramene/subsites/sorghum/v11/gramene-swagger-genotyping; cd $W
source ~/.nvm/nvm.sh && nvm use 24.14.1
git -C ../gramene-swagger status --porcelain        # unchanged (empty) before and after
pm2 describe sorghum_swagger11 | grep -E 'pid|uptime'   # unchanged pid
node --test "test/primers/unit/**/*.test.js"        # 0 failures; counts = baseline (503/487/16) + new tests
node --test test/primers/unit/worker_memory.test.js # 3 consecutive green runs after the idle() fix
PRIMERS_REALDATA=1 node --test --test-concurrency=1 "test/primers/unit/**/*.test.js"
# start the 50112 API and worker (§6.4), then:
PRIMERS_IT_BASE=http://127.0.0.1:50112/sorghum_v11 node --test --test-concurrency=1 "test/primers/integration/**/*.test.js"
BASE=http://127.0.0.1:50112/sorghum_v11
curl -s "$BASE/primers/genomes?system_name=sorghum_bicolor" | jq '.variation, (.genomes[]|select(.system_name=="sorghum_rio")|.has_variation)'
curl -s "$BASE/primers/variants?system_name=sorghum_bicolor&region=1&start=11180&end=11290" | jq '[.variants[].key]'
  # ["1:11182:A:G","1:11193:C:T","1:11203:C:T","1:11282:CA:C"]
curl -s -H 'Content-Type: application/json' -d '{"system_name":"sorghum_bicolor","variant":{"id":"rs871475760","alt":"A"},"assay":{"type":"kasp","num_sets":2}}' \
  "$BASE/primers/genotyping/design" | jq '[.sets[].key]'        # ["f9df650ad116","1accc54c262d"]
curl -s -H 'Content-Type: application/json' -d '{"system_name":"sorghum_bicolor","variant":{"id":"tmp_1_11502_C_CGT"}}' \
  "$BASE/primers/genotyping/design" | jq '.orientations.forward.status'   # "blocked"
redis-cli -p 6380 -n 1 --scan --pattern 'primers:sorghum_v11_geno:*' | head   # only the dev site key
ss -ltnp | grep -E '50011|50111'    # same pids as before
```

**Pass criteria:** all suites green against the recorded baseline; the keys and statuses as commented; the live checkout, `sorghum_swagger11` and the two pm2 dev apps untouched; no job under another site key; and the slot-starvation test showing zero `503 BUSY` for ordinary designs.
---

## 8. Front-end spec (for the gramene-primers and gramene-search session)

§2 is authoritative. The UI never derives roles, allele bases, dyes, tails, oligo names, order rows or allele calls — it displays them, and matches check results by target sequences.

### 8.1 Client methods (`src/client.ts`, `PrimersClient` in `src/types.ts:638-650`)

```ts
listVariants(q: VariantListQuery, o?: RequestOptions): Promise<VariantListResponse>;                     // GET /primers/variants
getVariant(id: string, q: { system_name: string }, o?: RequestOptions): Promise<VariantLookupResponse>;  // GET /primers/variants/{encodeURIComponent(id)}
designGenotyping(req: GenotypingDesignRequest, o?: RequestOptions): Promise<GenotypingDesignResponse>;   // POST /primers/genotyping/design, "design" timeout (60 s)
```

All three are **optional members** of `PrimersClient`, so existing custom clients and test doubles still compile; the Genotyping tab is disabled when `designGenotyping` is missing. `listVariants` sends `types` as CSV, sends `include_ems` only when `false`, and memoizes per `system_name|region|start|end|types|include_ems|limit` for 10 minutes, evicting on error like `listGenomes`; a cancelled caller does not cancel the shared request. Response guards mirror `client.ts:174-181`: `designGenotyping` requires `template` to be an object and defaults `sets`/`warnings` to `[]`; `listVariants` requires `variants[]`.

### 8.2 TypeScript types (`src/types.ts`, additive)

```ts
// Written out from the §2.6–§2.12 definitions; a type-parity test (§8.8) compares these keys with the swagger definitions.
// Reused unchanged from types.ts: DesignMode, PrimerWarning, CheckRequest, CheckName, CheckParams, RepeatMaskMode, RequestOptions.
export type DesignerMode = DesignMode | 'genotyping';          // DesignMode itself is NOT widened
export type VariantKind = 'snv' | 'mnv' | 'insertion' | 'deletion' | 'complex';
export type AssayType = 'kasp' | 'as_pcr';
export type TailScheme = 'none' | 'ref_fam_alt_hex' | 'ref_hex_alt_fam';
export type OligoRole = 'as_ref' | 'as_alt' | 'common';
export type SetOrientation = 'forward' | 'reverse';
export type SetQuality = 'good' | 'usable' | 'poor';
export type MismatchClass = 'max' | 'strong' | 'medium' | 'weak';
export type Likelihood = 'likely' | 'likely_weak' | 'unlikely';
export type IssueSeverity = 'info' | 'warn' | 'high';
export type GenotypeAllele = 'ref' | 'alt' | 'other' | 'ambiguous' | 'missing' | 'unavailable';
export type PrimerCallStatus = 'match' | 'weak' | 'terminal_mismatch' | 'uncertain' | 'blocked' | 'no_product' | 'unknown';
export type PredictedGenotype = 'ref' | 'alt' | 'both' | 'none' | 'no_call' | 'unknown';
export interface Span { start: number; end: number }

// ---- variants (§2.3, §2.4, §2.6)
export interface VariantListQuery { system_name: string; region: string; start: number; end: number;
  types?: VariantKind[]; include_ems?: boolean; limit?: number }
export interface VariantSource { name: 'ensembl'; release: string }
export interface VariantRecord { id: string; source: string; ems: boolean }
export interface VariantIssue { code: 'REF_MISMATCH' | 'STAR_ALLELE' | 'ALLELE_TOO_LONG' | 'UNSUPPORTED_ALLELE' | 'REPEAT_TOO_LONG';
  message: string; details?: Record<string, unknown> }
export interface AlleleSite { position: number; ref_base: string; alt_base: string; alt_maps_to: number | null }
export interface PrimerVariant { key: string; ids: string[]; synonyms: string[]; label: string; kind: VariantKind; region: string;
  vcf: { position: number; ref: string; alt: string }; minimal: { start: number; end: number; ref: string; alt: string };
  alleles: string[]; multiallelic: { alleles: string[]; other_alts: string[] } | null; shift: number;
  zone: Span; discriminating: { forward: AlleleSite; reverse: AlleleSite };
  records: VariantRecord[]; ems: boolean; consequence: string | null;
  ref_verified: boolean; designable: boolean; issues: VariantIssue[];
  requested_id?: string | null; submission_sequence?: string }
export interface VariantListResponse { system_name: string; region: string; start: number; end: number; source: VariantSource;
  total: number; returned: number; truncated: boolean; variants: PrimerVariant[]; warnings: PrimerWarning[] }
export interface VariantLookupResponse { requested_id: string; system_name: string; source: VariantSource;
  variants: PrimerVariant[]; warnings: PrimerWarning[] }

// ---- design request (§2.7)
export type VariantInput = { id: string; alt?: string } | { region: string; position: number; ref: string; alt: string };
export interface GenotypingAssay { type: AssayType; orientation: 'both' | 'forward' | 'reverse'; tails: TailScheme;
  deliberate_mismatch: 'none' | 'auto'; mismatch_position: 2 | 3; num_sets: number; max_relaxation: 0 | 1 | 2;
  neighbour_policy: 'avoid_3p' | 'ignore' }
export interface GenotypingParams { opt_size?: number; min_size?: number; max_size?: number;
  opt_tm?: number; min_tm?: number; max_tm?: number; opt_gc?: number; min_gc?: number; max_gc?: number;
  max_tm_diff?: number; max_poly_x?: number; gc_clamp?: number; max_end_stability?: number;
  salt_monovalent?: number; salt_divalent?: number; dntp_conc?: number; dna_conc?: number;
  product_size_ranges?: [number, number][] }
export interface GenotypingDesignRequest { system_name: string; variant: VariantInput; assay?: Partial<GenotypingAssay>;
  avoid_repeats?: boolean; repeat_mask_mode?: RepeatMaskMode; template_only?: boolean; label?: string; params?: GenotypingParams }

// ---- design response (§2.8)
export interface PrimerGenomic { region: string; start: number; end: number; strand: 1 | -1; blocks: Span[] }   // PrimerGenomicLocation
export interface NeighbourHit { key: string; ids: string[]; label: string; start: number; end: number; alleles: string;
  ems: boolean; distance_from_3p: number | null }
export interface LikelihoodAt { mm_pos: number[]; likelihood: Likelihood }
export interface DeliberateMismatch { position: 2 | 3; original_base: string; new_base: string; template_base: string;
  terminal_pair: string; terminal_mismatch_class: MismatchClass; added_mismatch_class: MismatchClass; source: string }
export interface TemplateSpan { start: number; end: number; sequence: 'ref' | 'alt' }
export interface TailedStructures { hairpin_th: number; self_any_th: number; self_end_th: number }
export interface GenotypingOligo { role: OligoRole; allele: string | null; three_prime_base: string; haplotype: 'ref' | 'alt' | 'both';
  target_seq: string; matched_seq: string; dye: 'FAM' | 'HEX' | null; tail_seq: string | null; order_seq: string;
  len: number; order_len: number; tm: number; tm_method: 'primer3' | 'ntthal_duplex'; matched_tm: number | null; gc: number;
  hairpin_th: number; self_any_th: number; self_end_th: number; end_stability: number; primer3_problems: string | null;
  template: TemplateSpan; genomic: PrimerGenomic; inserted_bases: number;
  deliberate_mismatch: DeliberateMismatch | null;
  discrimination: { own_allele: LikelihoodAt; other_allele: LikelihoodAt; terminal_mismatch_class: MismatchClass; in_shift_tract: boolean } | null;
  tailed: TailedStructures | null; neighbours: NeighbourHit[] }
export interface GenotypingProduct { size: number; template: TemplateSpan; genomic: PrimerGenomic; inserted_bases: number }
export interface PairThermo { compl_any_th: number; compl_end_th: number }
export interface TailedThermo { ref_alt_any_th: number; ref_alt_end_th: number; ref_common_any_th: number; ref_common_end_th: number;
  alt_common_any_th: number; alt_common_end_th: number }
export interface SetIssue { code: string; severity: IssueSeverity; message: string; details: Record<string, unknown> }
export interface OrderRow { name: string; set_id: string; set_key: string; role: OligoRole; allele: string | null; dye: 'FAM' | 'HEX' | null;
  order_seq: string; target_seq: string; tail_seq: string | null; length: number; tm: number; gc: number;
  orientation: SetOrientation; product_size_ref: number; product_size_alt: number; variant_key: string; notes: string }
export interface CheckPairInput { id: string; left: string; right: string; expected: { region: string; start: number; end: number } }
export interface GenotypingSetRef { id: string; ref_pair: string; alt_pair: string }
export interface GenotypingSet { id: string; key: string; rank: number; orientation: SetOrientation; relaxation_level: 0 | 1 | 2;
  quality: SetQuality; score: number;
  primers: { as_ref: GenotypingOligo; as_alt: GenotypingOligo; common: GenotypingOligo };
  products: { ref: GenotypingProduct; alt: GenotypingProduct };
  thermo: { ref_common: PairThermo; alt_common: PairThermo; tailed: TailedThermo | null };
  tm_balance: { as_tm_diff: number; common_minus_as: number };
  neighbour_sites: number; primer3_penalty: number; warnings: PrimerWarning[]; issues: SetIssue[];
  check: { set: GenotypingSetRef; pairs: CheckPairInput[] }; order: OrderRow[] }
export interface ExplainSide { raw: string; [reason: string]: string | number }
export interface GenotypingAttempt { level: number; changes: Partial<GenotypingParams>;
  explain: { left: ExplainSide; right: ExplainSide; pair: ExplainSide }; pairs_returned: number;
  rejected: { force: number; overlap: number; common_neighbour_3p: number; alt_scoring_failed: number; below_floor: number; duplicate: number };
  not_scored: number; sets: number }
export interface OrientationResult { status: 'ok' | 'no_sets' | 'blocked' | 'skipped';
  reason: 'neighbour_at_3p' | 'too_close_to_end' | 'n_in_primer_window' | 'not_requested' | 'below_floor' | 'budget_exhausted' | null;
  discriminating_position: number; relaxation_level: number | null; sets_found: number; blockers: NeighbourHit[]; attempts: GenotypingAttempt[] }
export interface GenotypingTemplate { system_name: string; region: string; start: number; end: number; strand: 1;
  length: number; alt_length: number; seq: string; alt_seq: string;
  masked: boolean; mask_source: 'softmask' | 'blast_depth' | null; mask: [number, number][]; masked_fraction: number;
  features: { variant: Span; zone: Span; discriminating: { forward: number; reverse: number }; alt_offset: number; exempt: [number, number][] } }
export interface KaspMix { stock_uM: number; as_ref_uL: number; as_alt_uL: number; common_uL: number; water_uL: number; total_uL: number; source: string }
export interface NeighbourSummary { data: 'ensembl' | 'none' | 'unavailable'; window: Span; variants: number; non_ems: number; ems: number; dense_non_ems: number }
export interface GenotypingSettings { preset: AssayType; params: GenotypingParams; pinned: string[];
  ladder: { level: number; changes: Partial<GenotypingParams> }[]; floors: { as_min_tm: number; as_min_gc: number } }
export interface GenotypingEngine { primer3: string | null; thermo: string | null; genotyping_design: string; variation_source: string | null }
export interface GenotypingDesignResponse { variant: PrimerVariant; template: GenotypingTemplate;
  assay: GenotypingAssay & { ems_target: boolean; kasp_mix: KaspMix | null };
  neighbours: NeighbourSummary; orientations: { forward: OrientationResult; reverse: OrientationResult } | null; sets: GenotypingSet[];
  check: { request: CheckRequest; set_ids: string[]; unique_primers: number; omitted_set_ids: string[] } | null;
  settings: GenotypingSettings; engine: GenotypingEngine; warnings: PrimerWarning[] }

// ---- check request and results (§2.11, §2.12)
export interface CheckGenotypingInput { variant: { region: string; position: number; ref: string; alt: string }; sets: GenotypingSetRef[] }
export interface GenotypeVariantInfo { key: string; region: string; position: number; ref: string; alt: string; shift: number; zone: Span;
  flank: number; core: { ref: string; alt: string }; haplotypes: { ref: string; alt: string } }
export interface GenotypeCopy { region: string; start: number; end: number; strand: 1 | -1; variant_position: number | null;
  identity: number; gap_compressed_identity: number; aligned_length: number; observed: string | null; flank_edits: number | null;
  call: 'ref' | 'alt' | 'other' | 'missing'; anchors: number; ortholog: boolean | null; source: 'amplicon' | 'megablast' }
export interface GenotypeGenome { system_name: string; display_name: string; is_reference: boolean; allele: GenotypeAllele;
  observed: string | null; source: 'amplicon' | 'megablast' | null; copies: GenotypeCopy[]; orthologous_copies: number; paralog_copies: number;
  reason: 'variant_not_covered' | 'no_orthologous_copy' | 'fallback_failed' | 'fallback_budget' | 'db_unavailable' | 'blast_error' | 'call_failed' | null }
export interface GenotypePrimerCall { status: PrimerCallStatus; likelihood: Likelihood | null; mm_pos: number[] | null; residual_mm_pos: number[] | null }
export type GenotypeReason = 'common_primer_3p_mismatch' | 'common_primer_weak' | 'ref_signal_off_locus' | 'alt_signal_off_locus'
  | 'shift_tract_uncertain' | 'approx_alignment' | 'no_orthologous_copy' | 'third_allele';
export interface GenotypeSetGenome { system_name: string; ref_primer: GenotypePrimerCall; alt_primer: GenotypePrimerCall; common_primer: GenotypePrimerCall;
  predicted: PredictedGenotype; strength: 'normal' | 'weak' | null; agrees: boolean | null; reasons: GenotypeReason[]; off_locus_products: number }
export interface GenotypeSetSummary { genomes_total: number; predicted_ref: number; predicted_alt: number; both: number; none: number;
  no_call: number; unknown: number; weak: number; agree: number; disagree: number; not_comparable: number }
export interface GenotypeSetResult { id: string; ref_pair: string; alt_pair: string; orientation: SetOrientation; deliberate_mismatch_positions: number[];
  specificity: { ref_pair: { verdict: string; consistent_with_allele: boolean }; alt_pair: { verdict: string; consistent_with_allele: boolean };
    off_target_count: number };
  control: { status: 'pass' | 'warn' | 'fail'; allele: GenotypeAllele; reasons: string[] };
  reference: GenotypeSetGenome; summary: GenotypeSetSummary; genomes: GenotypeSetGenome[] }
export interface GenotypingResults { algorithm_version: string; variant: GenotypeVariantInfo;
  summary: { genomes_total: number; ref: number; alt: number; other: number; ambiguous: number; missing: number; unavailable: number };
  genomes: GenotypeGenome[]; sets: GenotypeSetResult[] }

// ---- UI state (§8.3): what was submitted, so results match by sequence (§8.6)
export interface SubmittedGenotypingSet { id: string; set_key: string; ref_pair: string; alt_pair: string;
  primers: { as_ref: string; as_alt: string; common: string } }   // uppercase target sequences
// CheckRequest (types.ts:297-306) gains:  genotyping?: CheckGenotypingInput;
// CheckResults  (types.ts:546-558) gains: genotyping?: GenotypingResults;
// GenomeEntry   (types.ts:241-253) gains: has_variation?: boolean;
// GenomesResponse (types.ts:255-260) gains: variation?: { available: boolean; source: 'ensembl' | null; release: string | null };
```

**`DesignMode` is deliberately not widened.** Widening it would leak `'genotyping'` into `CheckRequest.mode`, `PrimerTemplate.mode` and `estimateCheckCpu`, and `PrimerDesigner` would then be able to send `mode: 'genotyping'` to `/primers/check`, which swagger rejects with `ENUM_MISMATCH`. Instead `DesignerMode = DesignMode | 'genotyping'` is introduced **and applied consistently** to `PrimerDesignerProps.modes`/`defaultMode` (`types.ts:712-713`), `DesignerContext.modes`/`defaultMode`, `ALL_MODES` (`state.ts:19`), `availableModes` (`state.ts:42-52`) and `initialDesignerState`; otherwise a TypeScript host cannot pass the new mode. The widened exported types are recorded in the 0.2.0 CHANGELOG.

Reader tolerance stays as documented (`types.ts:1-5`): unknown fields are ignored, a missing `genotyping` block means "not a genotyping job", a missing `has_variation` means `false`.

### 8.3 State additions (still `v: 1`)

```ts
export interface GenotypingState {
  variantId?: string; alt?: string; variantKey?: string;
  manual?: { region: string; position: number; ref: string; alt: string };
  window?: { region: string; start: number; end: number };            // <= 50,000 bp
  filters?: { types?: VariantKind[]; includeEms?: boolean; query?: string };
  assay?: Partial<GenotypingAssay>; params?: Partial<GenotypingParams>;
  label?: string; avoidRepeats?: boolean; repeatMaskMode?: RepeatMaskMode;
  designed?: boolean; selectedSetKey?: string; checkedSetKeys?: string[];   // <= 5, /^[0-9a-f]{12}$/
  check?: { checks: CheckName[]; genomes?: string[]; params?: Partial<CheckParams>; jobId?: string; submitted?: SubmittedGenotypingSet[] };
  view?: { tab: 'sets' | 'alleles' | 'specificity' | 'pangenome' | 'order' };
}
// PrimerDesignerState: mode: DesignerMode;  genotyping?: GenotypingState;
```

`normalizeDesignerState` (`state.ts:145-211`) gains `cleanGenotyping(raw)`: ids are checked against `^[A-Za-z0-9][A-Za-z0-9_.:,*-]{0,254}$` (the same pattern as the server — a narrower one would drop about 2 % of real Ensembl ids), alleles against `^([ACGT]{1,50}|-)$`, keys against `^[^:]+:\d+:[ACGT]+:[ACGT]+$`, `jobId` against `^[0-9a-f]{32}$`; integers must be positive and the window ≤ 50,000 bp; invalid fields are dropped. An older build reading `mode: 'genotyping'` falls back to its default mode and ignores the slice [R].

### 8.4 Mode wiring

| Place | Change |
|---|---|
| `components/util.ts:46-51` `MODE_LABELS` | becomes `Record<DesignerMode, string>` with `genotyping: 'Genotyping'` |
| `ModeTabs.tsx:16-36`, `PrimerDesigner.tsx:338` | the tab is disabled with the reason "choose a genome" when no genome is known, and when the client lacks `designGenotyping` |
| `PrimerDesigner.tsx:414-462` | a new branch renders `<GenotypingPanel>`; `RepeatOptions` is reused; `ParamsPanel` takes `presets={['kasp','as_pcr']}` and the genotyping key list (it hard-codes `['pcr','qpcr']` today at `ParamsPanel.tsx:188`) |
| `src/components/designChecks.ts:55-75`, `util.ts:71-88` | a genotyping case: inputs are ready when a designable entry is chosen or the manual fields validate; estimated template length = 2 × 400 + allele lengths |
| `presets.ts` | `GENOTYPING_PRESETS` mirror the §4.6 level-0 tables; the ladder is display-only; `changedGenotypingParams` diffs against the assay type's level 0 |
| `validate.ts` | `validateVariantInput` (the two styles of §2.1, ≤ 50 nt, ref ≠ alt) and `validateGenotypingParams`, including the hint that the effective product minimum is `2 × max_size + 1` |
| `request.ts` | `buildGenotypingRequest(state, ctx)` sends only changed assay keys and params; `buildGenotypingCheckRequest` (§8.6) |

### 8.5 Components

**`GenotypingPanel`** — genome line with the `variation` badge ("Ensembl 115 variants" / "manual entry only"); `VariantPicker` or `ManualVariantInputs`; the selected-variant card (label, key, ids with "+n", EMS badge, consequence, multi-allelic alt chooser, `issues`, and `shift` → "indel can slide 2 bp"); `AssayOptions`; buttons **Check reference** (`template_only`) and **Design sets**.

**`VariantPicker`** — window defaults to the gene span ± 2 kb clamped to 50 kb, editable with the server's own validation; `listVariants` debounced 300 ms and aborted on change; a virtualized table with position, label, kind, ids, source badges (EMS hollow), consequence and status (designable ✓ or the issue code); keyboard up/down/enter; filters for kind, "hide EMS" and a client-side text search. **Selecting a row stores `variantKey` and `vcf`, and designs from `{region, position, ref, alt}`** — the id travels only as a label, so ids containing `,` or `*` can never break a design. `truncated: true` shows "showing 2,000 of N; narrow the window". Errors map as in §1.6, with a countdown from `retry_after_s`.

**`ManualVariantInputs`** — region, position, ref, alt, with the hint "use `-` for an inserted or deleted allele; for an insertion, position is the base after the insertion point"; **Check reference** calls `template_only` and shows either "genome has C at 1:11109 ✓" plus the normalized label/key/ids, or marks the ref field from `details.genome`.

**`AssayOptions`** — assay type (KASP | Gel AS-PCR), orientation, tails (hidden for as_pcr behind "advanced"), deliberate mismatch and its position (−2 / −3, with −3 labelled "extrapolated"), number of sets (1–10), max relaxation, neighbour policy, avoid repeats, and a collapsible `ParamsPanel` marking pinned params "not relaxed".

**`SetsTable`** (replaces `PairsTable` in this mode) — one row per set keyed by `set.key`, with: checkbox (disabled when it would exceed 5 sets, 10 pairs or 20 distinct primers, or when the mirrored `estimateCheckCpu` for the currently chosen genomes is `over_limit`, with the live cost line — over three genomes all 20 primers fit, over the full panel the limit bites at 14); set id, orientation arrow and level badge; REF primer (`target_seq`, 3′ base highlighted, dye chip, copy buttons for target and order sequences); ALT primer; common primer; product `ref/alt` bp; Tm REF/ALT/common with ΔTm; a neighbours badge from `neighbour_sites`; quality chip (`good`/`usable`/`poor`) and warning chips by code (text, never colour alone); and, after a check, the agreement chip (e.g. "119/119 agree") plus the two specificity verdicts. A deliberate-mismatch primer underlines the substituted base with the tooltip "−2 A→G (weak extra; terminal mismatch maximum) — Little 1995 Table 9.8.1".

**`OrientationExplain`** — one card per orientation with `status`, `blockers` ("blocked by rs5413863234 (G/A) 4 nt from the 3′ end") and the level-by-level `attempts` table reusing the existing explain renderer.

**`SetDetail`** — per oligo: role, allele detected, 3′ base, genomic footprint (blocks and inserted bases), Tm with method, GC, structures, tailed structures, discrimination ("own allele: likely; other allele: 3′-terminal mismatch → weak; terminal class maximum"), neighbours with distances; plus set dimers and the three order rows with the KASP mix note.

**`TemplateMap`** (reused) — a variant marker at `template.features.variant` with separate F and R discriminating ticks when they differ, a shaded band for `template.features.zone`, neighbour ticks from `listVariants` (filled = non-EMS, hollow = EMS), and one lane per set drawing the REF and ALT arrows stacked on the same footprint with their 3′ tips coloured by allele, plus the common-primer arrow facing the other way. This replaces the one-left/one-right assumption at `TemplateMap.tsx:34-49` and `:374-378` for this mode; a blocked orientation is hatched.

**`AlleleMatrix`** (results tab "Alleles", reusing the grid keyboard handler of `PangenomeMatrix.tsx:307`; `:49-79` is `buildMatrixRows`, the row model) — rows: the reference control first, then genomes in results order; columns: the allele, then one per submitted set. The allele cell shows a glyph, colour **and** text (REF ● "C", ALT ▲ "A", other ◆ with the observed core, ambiguous ?, missing –, unavailable ×) with the copy count as a superscript. The set cell shows the predicted signal (REF ●, ALT ▲, both ◐, none ○, **no_call ⊘**, unknown ·), a ✓/✗ for `agrees`, a "weak" marker when `strength: "weak"`, and ⚠ when `off_locus_products > 0`; the tooltip lists `reasons`. The detail panel shows the copies table (region:start–end, strand, identity and gap-compressed identity, the observed core aligned under `variant.core.ref`/`.alt`, call and source), `paralog_copies`, and per set the three primer statuses with `mm_pos` and `residual_mm_pos`. A "disagreements only" filter selects `agrees === false`.

**Exports** (`ExportMenu` items built at `PrimerDesigner.tsx:320-334`): order sheet TSV and CSV (the `order` rows of the checked sets plus `#` note lines for the KASP mix and the submission string), genotype calls TSV (`genomes × sets`), the KASP service submission line, and a primers FASTA named from `order[].name`, which is unique by construction and fixes the duplicate `G_P1_F` problem [R]. Every cell passes through `tsvCell` (`exporters.ts:24-32`).

### 8.6 Check integration

```ts
export function buildGenotypingCheckRequest(i: {
  design: GenotypingDesignResponse; setKeys: string[]; checks?: CheckName[];
  mode?: 'region' | 'gene'; geneId?: string; genomes?: string[] | null;
  allGenomes?: GenomesResponse | null; params?: Partial<CheckParams> | null;
}): CheckRequest;
```

`pairs` is the concatenation of the chosen sets' `check.pairs` in rank order; `genotyping.sets` their `check.set`; `genotyping.variant` is `{region, ...design.variant.vcf}`. Tails are never sent — only `target_seq`, which already contains any deliberate mismatch. Limits `{maxSets: 5, maxPairs: 10, maxUniquePrimers: 20}` throw `CheckRequestError` with `TOO_MANY_SETS` / `TOO_MANY_PAIRS` / `TOO_MANY_PRIMERS`. A request whose mirrored `estimateCheckCpu({unique_primers, genomes, genotyping: true})` is `over_limit` throws `OVER_CPU_LIMIT`, so the pan-genome gate depends on the chosen genomes, never on a fixed primer count. The number 13 survives only as the server's cap for the full-panel `design.check.request` (§4.18). The simplest path is to post `design.check.request` unchanged.

**Matching results to sets never uses ids or ranks.** `matchGenotypingResults(sets, job, submitted)` builds the triple `pairKey(ref.left, ref.right) + '#' + pairKey(alt.left, alt.right)` from uppercase target sequences in `job.request.pairs` (grouped through `job.request.genotyping.sets`), and looks up each design set's identical triple; `results.genotyping.sets` is then attached by the submitted set id, and pair-level results by pair id. Unmatched design sets are `notChecked`; submitted triples with no design set are `orphanIds` ("results for primers not in the current design"). This preserves the rank-independence of today's `matchCheckResults` (`results.ts:38-81`) and survives a re-design that renumbers sets.

**Defensive state:** a `done` job whose `request.genotyping` is set but whose `results.genotyping` is missing means an older worker served it — show "Allele check unavailable (server update in progress)" with a Re-run button rather than an empty matrix.

### 8.7 Cost mirroring (`src/cost.ts`)

```ts
export const GENOTYPE_CPU_S_PER_GENOME = 0.2;   // mirrors check.genotype_cpu_s_per_genome
// CheckCpuInput gains  genotyping?: boolean;  CheckCpuEstimate gains genotyping_cpu_s: number;
// genotyping = input.genotyping ? (1 + pan.length) * GENOTYPE_CPU_S_PER_GENOME : 0;
```

Parity tests use the **real** genome sizes (reference 708,735,318 and the 119 others summing 83,774,241,795) and the exact server outputs, not rounded inputs: 6 primers full panel → 2,690; 12 → 5,356; 13 → 5,800; 14 → 6,245 (`over_limit`); 6 primers over the three §2.11 genomes (719,894,357 / 690,603,914 / 773,056,139) → 95 [V]. The UI text is "4 sets · 12 distinct primers · ≈ 5,356 CPU-s".

### 8.8 Fixtures, playground and gramene-search

- **Offline fixtures** (`test/fixtures/genotyping/`): the §2.3, §2.4, §2.9, §2.10 and §2.13 documents copied verbatim, plus a partial job (reference + one genome). The design fixtures are real data; the illustrative blocks of §2.13 are replaced with a real 50112 job before hand-off (§7.6).
- **`fakeClient.ts`** gains `listVariants`, `getVariant` and `designGenotyping` serving those fixtures, plus error modes `VARIATION_SOURCE_UNAVAILABLE`, `NO_VARIATION_DATA` and `REF_MISMATCH`.
- **Component tests:** a genotyping `describe` in `designer.test.tsx` (pick → design → check → alleles → export); order-sheet TSV regexes in `check.test.tsx`; `AlleleMatrix`, `SetsTable` and the genotyping lanes in `parts.test.tsx`; the new components in `a11y.test.tsx`; the mode set and grid in `playground.test.tsx`.
- **Playground:** `pages.ts` gains `genotyping-rs871475760`, `genotyping-insertion-blocked`, `genotyping-no-variation` and `genotyping-ensembl-down`; `examples/playground/mockClient.ts:225-249` gains the three methods and a job simulation that emits the reference and then three genomes at 1 s intervals; `examples/playground/App.tsx:120` includes the new mode.
- **Type parity:** a `types.parity.test.ts` loads the gramene-swagger definitions of §2.6–§2.12 from the fixture manifest and asserts that every property name of each definition appears in the matching §8.2 interface and vice versa.
- **gramene-search:** `src/components/results/details/Primers.js` adds `'genotyping'` to `MODES` (`:22`, `:56-57`); no new props, because the component reads `variation` from `listGenomes`. `uiViewState.js` and `viewSnapshot.js` need no change. `demo.js:146` still points `ensemblRest` at release 108, which is irrelevant here because variants come through swagger. A VEP-tab deep link is deferred (§9.3).

### 8.9 Acceptance checklist for the UI session

1. On SORBI_3001G000200, Genotyping mode lists the four variants of §2.3 in 1:11180–11290; hiding EMS leaves two.
2. Designing rs871475760 as KASP with 2 sets shows S1 reverse (65 bp, `usable`, score 9.18) and S2 forward (92 bp, `poor`, score 19.63, a high tail-hairpin chip on the ALT primer plus `ALT_PRIMER_SUBOPTIMAL` and `NEIGHBOUR_IN_PRIMER`); the map shows two lanes, the variant marker and the zone band.
3. tmp_1_11502_C_CGT shows the forward orientation blocked with the neighbour text, and S1 with the `AS_TM_IMBALANCE` chip.
4. Gel AS-PCR shows the mismatch description and `tm_method` "template duplex".
5. Checking S1 and S2 posts exactly the §2.11 body; the cost line shows 95 CPU-s for three genomes and 2,690 for the full panel.
6. With the fixture job, the Alleles tab shows the reference control, pi180348 with two copies, and 3/3 agree for both sets; "disagreements only" empties the grid.
7. A set whose common primer has a 3′ mismatch in one genome shows ⊘ `no_call` there, with the reason in the tooltip.
8. The order-sheet TSV for S1 equals the §1.7 rows plus the note lines; FASTA names are unique.
9. A saved-view round trip keeps the variant, assay, checked set keys and job id; an old build ignores the slice.
10. `npm run typecheck && npm test && npm run build && npm run lint:pkg` pass and the new components are axe-clean.
---

## 9. Risks, open questions and deferred items

### 9.1 Risks

| Risk | Impact | Mitigation |
|---|---|---|
| KASP presets and the ladder are validated in silico only; LGC's own design rules (Kraken) are unpublished [R] | unknown wet-lab failure rate | the response reports the relaxation level, tail structures and every issue; the per-assembly check exposes neighbour and paralog problems before ordering; revisit the presets after the first plates (§11 asks the user to confirm the assay defaults) |
| Tailed-oligo hairpins are routinely above 47 °C — the HEX tail alone gives 51.50 °C on an ordinary primer [V] | warning fatigue | two thresholds: warn at 47 °C, high at 55 °C (the final KASP annealing temperature); a warn issue adds only 1 to the score and a high issue 3, so routine warn-level hairpins barely move the ranking |
| Ensembl REST availability or format drift | pickers unusable | limiter, single flight, breaker, negative caching, strict response validation (a bad body becomes a source error, never wrong data), and manual entry always works; the fake-Ensembl integration modes pin the behaviour |
| Neighbour data exists for sorghum_bicolor only, and the EMS exemption is an assumption about germplasm | undetected neighbours elsewhere; an EMS neighbour could matter in a mutant × mutant cross | `neighbours.data` says which case applies; EMS neighbours are still reported per primer; for an EMS **target** natural neighbours only stop blocking, they are still reported and penalised |
| Shiftable indels and homopolymers: TP-ARMS succeeded in under half of such assays [R] | poor discrimination despite an in-silico `likely_weak` | `SHIFTABLE_INDEL` and `SHIFT_TRACT_DISCRIMINATION` warnings; the allele-specific primer that ends inside the shift tract (§3.5, §4.12) is reported `uncertain` on the other allele's genomes, which then predict `unknown` and are excluded from concordance |
| Haploid assemblies cannot show heterozygosity, and some assemblies are mosaics (pi154987 and pi656029 carry two haplotypes [R]) | a `ref`/`alt` call for a line that is heterozygous in the field | documented in the contract; copies that disagree give `ambiguous`; Ensembl germplasm genotypes are deliberately not shown |
| Ortholog thresholds (gap-compressed identity 95 %, size ±20 %) may still send a diverged ortholog to the fallback (pi536008 at 78.2 % raw identity) | extra CPU, possibly `missing` | the fallback exists for exactly this; the real-data expectation is marked [U] and the thresholds are config keys |
| The genotyping cost coefficient (0.2 CPU-s per genome) is not measured end to end | the estimate is off by a few CPU-s | small next to BLAST; measured in integration and adjustable through `check.genotype_cpu_s_per_genome` |
| The submit-time FASTA read adds a failure mode to the API process | +≤ 5 ms per genotyping submit | reuses the sequence handle LRU; errors map to existing codes |
| `SEQUENCE_TARGET` shrinks the common-primer search space near a template edge | fewer candidates for variants close to a contig end | the gap shrinks before it blocks, and the orientation is reported `skipped {too_close_to_end}` rather than silently failing |
| A genotyping check consumes a large share of the 6,000 CPU-s budget: only 4 sets fit the full panel | users may want more | `check.omitted_set_ids` names what was left out, and the client mirrors the cost so Submit is disabled before the server refuses |
| The `results.genotyping` document grows with genomes × sets | large job documents | ≈ 220 kB for 120 genomes and 5 sets, well inside `max_result_bytes` [I] |
| Only 8 candidates per orientation are scored, so a better-scoring pair further down Primer3's 20 can be missed | a slightly worse best set | Primer3 returns pairs in penalty order and the penalty is the largest score term; `not_scored` shows the truncation per attempt; `max_scored_per_orientation` is a config key, and each extra candidate costs at most 3 `check_primers` runs and 17 `ntthal` calls, so the budget caps must be raised with it (§4.8) |

### 9.2 Open questions for the user

The six defaults of §11 are the decisions that need a yes/no. Beyond those:

1. **`expected` for hand-pasted sets.** v1 requires `expected` on every pair a set uses. Should the server instead derive it from the reference alignment, so primers designed elsewhere can be checked?
2. **Gene mode as the UI default.** In the gene tab, should the check run in `gene` mode (ortholog-annotated primary products) rather than `region`?
3. **Wet-lab confirmation.** Ordering rs871475760 S1 and genotyping BTx623 (expect FAM only), PI510757 (HEX only) and Rio (HEX only) would validate the whole chain; worth doing before the presets are frozen.

### 9.3 Deferred (explicitly not in v1)

- Tetra-primer ARMS and CAPS/dCAPS (excluded by the user's assay decision).
- Listing Ensembl germplasm genotypes (excluded by the user; used only to cross-validate in tests).
- **Tm rebalancing by moving a 5′ end** by up to ±3 nt. It measurably helps (1.16 → 0.13 °C on tmp_1_11502 in a prototype [R]) but gives the two pairs of one set different `expected` windows, which the submit-time validation forbids; revisit together with that rule.
- **A design-time paralog megablast.** The check's specificity and pan-genome stages already find paralogs, and the design-time scan would add a second BLAST per design for information the user gets a minute later.
- IUPAC/degenerate bases in primers for dense neighbourhoods: Primer3 never places a primer over an IUPAC base [R], so this needs post-design substitution and rescoring.
- Batch design (many variants, one order sheet), multi-allelic 4-dye assays, and heterozygosity prediction.
- The gramene-search VEP-tab deep link, and service-specific submission formats beyond the `[X/Y]` string.
- Caching design responses, and per-user quotas.
---

## 10. Implementation order

Each milestone is independently testable and small enough to review on its own: a milestone's tests use only the code of that milestone and earlier ones. Nothing is committed, pushed or deployed without the user's explicit go-ahead.

| # | Milestone | Contents | Done when |
|---|---|---|---|
| **M0** | Worktree and baseline | create `gramene-swagger-genotyping` from `primer-design` (§6.4), copy `node_modules`, record the unit baseline (503 / 487 pass / 16 skipped [V]) | the suite runs green in the new worktree and the baseline is written into the checklist |
| **M1** | Flaky-test fix | reorder `idle()` in `jobs/worker_main.js` so the SIGTERM/SIGINT handlers precede the "idling" log line; tighten `worker_memory.test.js:890` to match `{code: 0, signal: null}` (§7.1) | that file passes 3 consecutive runs and the full suite stays green; the race is documented in the commit message |
| **M2** | Compatibility tests | C1–C6 of §7.2, with the 11 job ids and the design-normalize snapshots recorded from `primer-design` HEAD | all six pass **before** any feature code exists |
| **M3** | Variant resolution, pure half | `variation/normalize.js` (left normalization, shift, discriminating positions, zone, kinds, labels, submission strings) with the §7.3 table as its test | every row of §3.4 and §7.3 passes offline, with no network and no FASTA beyond the fixture window |
| **M4** | Variant resolution, I/O half, and its HTTP surface | `variation/client.js` (limiter, single flight, breaker, LRU, byte caps, status mapping, lenient record ids) and `variation/index.js`; config keys `variation.*` with `PRIMERS_VARIATION_URL` / `PRIMERS_VARIATION_ENABLED` (`config.js`, `config/default.yaml`); `genomes.js` `has_variation` and `variation`; controller handlers `listPrimerVariants` and `getPrimerVariant`; **swagger:** the two variant paths, the §2.6 definitions, `PrimerWarning.details` and the genomes fields; **`contract.test.js`:** the two operationIds, `variants-*.json` routing and hand-written provisional `variants-*.json` fixtures (§7.7); `genomes_response.test.js`, `config.test.js`; **`controller_wrap.test.js`:** the two variant handlers (§7.3) | the fake-fetch tests of §7.3 pass; the two endpoints answer on 50112 with the §2.3 and §2.4 examples; sway reports no new warning (no `UNUSED_DEFINITION`); `contract.test.js` and `controller_wrap.test.js` pass |
| **M5** | Design core | `genotyping/{request,presets,template,guard}.js` and the orientation runner: template window, zone, `SEQUENCE_TARGET` guard, mask exemption, ladder with hard floors, the §4.8 filters with the scoring cap and budget reservation, explain reporting — Primer3 real, thermo stubbed, ALT derivation not yet; config keys `genotyping.*`; `template.js` exports `buildRegionTemplate` / `requireRegionLength`; `design.js` exports `parseParams`; the recorder `test/primers/fixtures/primer3/record_genotyping.js` (§7.4) | the rs871475760, tmp_1_11193_C_T, rs5413864115 and tmp_1_11502_C_CGT records reproduce the recorded inputs byte for byte, and the pairs, explain strings, `pairs_returned` and `not_scored` of §2.9–§2.10 at `PRIMER_NUM_RETURN` 20 |
| **M6** | Sets, scoring and the design endpoint | `sets.js`, `mismatch.js`, `scoring.js`, `thermo.js`, `order.js`: ALT derivation, the literal Table 9.8.1, `check_primers` scoring, floors, discrimination, tails, issues/quality/score/ranking, order rows, `check.request`; config key `ntthal` / `NTTHAL` (with its `config.test.js` case); `primer3.js` exports for `thermo.js`; the recorder extended to scoring runs and thermo; controller handler `designGenotypingPrimers`; **swagger:** the `/primers/genotyping/design` path and the §2.7–§2.8 definitions; **`contract.test.js`:** the operationId, the strict definitions `PrimerGenotypingRequest`, `PrimerVariantInput`, `PrimerGenotypingAssay`, `PrimerGenotypingParams`, the §2.9–§2.10 design request bodies (`:136-187`), the design half of the sync tests (`:317-364`; the `BODY_KEYS` assertion belongs to M7), `genotyping-design-*.json` routing (`:401-419`, `:447-479`) and provisional fixtures; **`controller_wrap.test.js`:** `designGenotypingPrimers` (§7.3) | `POST /primers/genotyping/design` on 50112 returns the §2.9 response; the eight keys and scores of §4.16 reproduce; the budget, rounding, score-rule and shift-tract tests of §7.3 pass; sway reports no new warning; `contract.test.js` and `controller_wrap.test.js` pass with no M7 code present |
| **M7** | Check extension, submit side | `check/genotype.js` validation half, `check/normalize.js` wiring, `algorithmVersionFor` and the guarded `jobs/index.js:164` call, the cost term; **swagger:** `PrimerCheckRequest.properties.genotyping` and the three §2.11 definitions; **`contract.test.js`:** those definitions in the strict list, the §2.11 check body (`:136-187`), the check half of the sync tests (`BODY_KEYS` contains `genotyping`, `:366-387`) and `check-genotyping-*.json` routing with provisional fixtures; **`job_id.test.js`:** the `algorithmVersionFor` case (§7.3) | C1 still passes (existing job ids unchanged); the §2.11 body is accepted over HTTP on 50112 and every §7.3 rejection reason fires; `cost.test.js` and `job_id.test.js` pass with and without genotyping; sway and `contract.test.js` clean |
| **M8** | Check extension, worker side | `check/genotype.js` caller half plus the `run.js` integration: anchors, alignment, exact-core calls, copy merging, ortholog test, megablast fallback with its budget, three-primer prediction, reference control, summaries, partial flushes; `check/blast.js` megablast args and parser; config keys `check.genotype_*`; **swagger:** `PrimerCheckResults.properties.genotyping` and the §2.12 definitions | the fake-world test of §7.3 gives ref / alt / other / missing / 2-copy alt / paralog-excluded / `no_call` for the KASP and AS-PCR sets, and the 11-assembly real-data table reproduces; sway and `contract.test.js` clean |
| **M9** | Docs, spec and example validation | `docs/primer_design_api.md`; this spec copied to `docs/genotyping_design_spec.md`; the blocking example-validation test (the `x-nullable`-aware walker over every JSON example of the spec); the gramene-primers fixture generator, as the hand-off of §7.7 — the API session writes no gramene-primers code, and once the UI session has run `npm run fixtures` it replaces its provisional request fixtures with the generated ones | every documented example validates; the generated request fixtures round-trip through both normalizers; no swagger or contract-test change is left for this milestone |
| **M10** | Integration on 50112 | the dev API and worker, the loopback fake-Ensembl server `test/primers/integration/fake_ensembl.js` and its modes (commands in §6.4), the full §7.6 table including the slot-starvation test | every row passes; no ordinary design receives `503 BUSY` while Ensembl is slow |
| **M11** | Real data | `PRIMERS_REALDATA=1` suites, including the 119-assembly caller run and the Ensembl cross-validation | §7.5 expectations hold; pi536008 [U] is resolved either way and recorded |
| **M12** | Hand-off | freeze §8 as the front-end brief, replace the illustrative blocks of §2.13 with a real 50112 job, publish the fixture set; **then, only with the user's go-ahead**, repoint the pm2 dev apps (worker first) so the UI session can build against `sorghum_v11a` | the UI session has real fixtures and a live endpoint; the live checkout, `sorghum_swagger11` and port 50011 are untouched throughout |

M3/M4 and M5/M6 are the natural review boundaries; M7 and M8 can be reviewed separately because the submit side is pure and the worker side is not.

---

## 11. Defaults to confirm with the user

| # | Default | Recommended value | One-line rationale |
|---|---|---|---|
| 1 | Which allele gets FAM | **REF = FAM, ALT = HEX** (switchable per design) | LGC's manual labels the allele-1 tail FAM and plots FAM on the X axis, and REF is allele 1; PolyMarker's opposite convention stays available as `ref_hex_alt_fam`. |
| 2 | Deliberate-mismatch position | **−2** (penultimate), with −3 selectable and labelled extrapolated | Little's Table 9.8.1 is *defined* for the penultimate base, so −2 is the only position with a sourced base table; Liu 2012's better empirical rate at −3 comes with no usable pair notation. |
| 3 | Check defaults for an EMS target | **specificity + reference control only** (pan-genome off by default, one click to enable) | An EMS allele is private to one BTx623-background line, so the 119-assembly allele matrix is nearly all `ref` and costs thousands of CPU-s; the user may prefer it on when EMS lines are crossed out. |
| 4 | Number of sets returned | **6** (`num_sets`, 1–10) | Enough to choose from in the table, while only 4 sets (12 primers) fit a full-panel check at 5,356 CPU-s against the 6,000 limit. |
| 5 | A non-EMS neighbour inside the common primer's last 5 nt | **reject that pair** (the orientation continues with other candidates) | A 3′ mismatch under the common primer silences *both* alleles, which is a no-call rather than a wrong call; with `num_return 20` there are usually alternatives. Warning-only is the alternative if too many loci lose their best set. |
| 6 | `avoid_repeats` for genotyping designs | **false in v1** (same as `/primers/design`), recommended `true` once reviewed | Measured on the rs871475760 window: enabling it masks 200 of 801 bases (`blast_depth`, 24.97 %) and would remove the verified gel AS-PCR common primer at 1:10880–10903 [V] — good for locus uniqueness, but it changes which sets exist, so the user should choose knowingly. |

---

## 12. Decision log

| Decision | Alternatives considered | Why the alternatives were rejected |
|---|---|---|
| Dedicated endpoints (`/primers/variants`, `/primers/variants/{id}`, `/primers/genotyping/design`) | a `variant` mode on `/primers/design`; variants as a side effect of `template_only` previews | A mode widens `PrimerDesignRequest.mode` **and** leaks into `PrimerCheckRequest.mode` and the client's `DesignMode`; per-mode required fields cannot be expressed in a strict swagger definition; indel ALT pairs in `pairs[]` break `product_size = right.end − left.start + 1` and `len = end − start + 1`, which the existing UI renders. Serving variants through a design preview would also hold a Primer3 slot for a 1 MB Ensembl fetch. |
| Sets as first-class objects with four sequence fields per oligo | returning pairs and letting clients assemble sets | Every consumer would re-derive roles, dyes, tails and order names, and tails would eventually leak into check requests (they fail the check's `^[ACGT]{15,36}$` pattern). |
| Canonical VCF `key` as identity; ids for lookup and display | keying designs and jobs on the Ensembl id | 88 of 3,976 real ids in one 50 kb window contain `,` or `*` or exceed 128 characters [V], and the same event appears under several ids (`tmp_1_11502_C_CGT` = `rs5413863549`). |
| Id pattern `^[A-Za-z0-9][A-Za-z0-9_.:,*-]{0,254}$`, lenient validation for record ids | narrower patterns (a 128- or 200-character limit, or no `,` and `*`) | They reject 53 and 88 of those ids respectively [V]; on the record side that silently removes real neighbours from screening, which is a correctness bug, not a cosmetic one. |
| `SEQUENCE_TARGET` guard just outside the zone, computed from both ALT anchor mappings | product-size minimum plus a REF-only post-filter; `SEQUENCE_PRIMER_PAIR_OK_REGION_LIST`; a second Primer3 step with the AS primer fixed | The product minimum misses multi-base shiftable indels (rs5413863413's zone reaches 13752, 3 bp past the REF span) and MNVs. `OK_REGION_LIST` works too but needs per-orientation index arithmetic and an interval cap; a second step doubles the Primer3 runs. The TARGET guard was verified to fix the real overlap case in one run [V]. |
| No `PICK_ANYWAY` in design runs; ladder floors Tm ≥ 52 °C, GC ≥ 15 % for derived primers | a `PICK_ANYWAY` rung with a ranking penalty | KASP anneals at 55 °C; a 45 °C primer produces wrong homozygous calls, not visible failures. Level 2 already yields a real 52.87 °C primer at the hardest example locus [V]. |
| Score ALT and mismatch primers with Primer3 `check_primers` on their own haplotype | `oligotm`/`ntthal` plus a JS re-implementation of `END_STABILITY` and the pair terms | The JS route re-implements published thermodynamics and diverged from Primer3 once on `PAIR_COMPL_END_TH` (13.18 vs 0). `check_primers` reproduces the design run exactly [V] and costs 9 ms. `ntthal` is kept only where Primer3 cannot help: tailed structures and mismatch duplex Tm. |
| Literal 32-cell Table 9.8.1 lookup | the class-rule shortcut | The shortcut differs from the table in 2 of 32 cells and changes a real answer: tmp_1_11193_C_T reverse at −2 needs REF A→C and ALT A→T, not C for both [V]. |
| Exact-core allele calling; anything else is `other` | lowest haplotype edit distance, with ties → `other` and a divergence cap | Distance names a third repeat length after the nearer allele: a two-A deletion at rs5413864115 is 1 edit from ALT and would be called `alt`; `AAAA` would be called `ref` [V]. |
| Prediction reads all three primers, plus off-locus products | reading only the allele-specific primer's `mm_pos` | A private SNP under the common primer's 3′ end gives a lab no-call; reporting that genome as a concordant REF/ALT call is the most damaging kind of wrong answer. Paralogs carrying the other allele light both dyes. |
| Ensembl work strictly before `semaphore.acquire`, under a limiter, breaker and single flight | resolving inside `runDesign` (where the semaphore is already held) | `design()` takes the slot at `design.js:671` before `runDesign`; four genotyping designs waiting 8 s each inside slots hand `503 BUSY` to every ordinary PCR user after their 10 s wait. |
| `request.genotyping` stored request-shaped; version via `jobId`'s `algo` (`'2+g1'`), guarded call site | adding `algorithm`/`orientation` to the stored request; bumping `ALGORITHM_VERSION` to `'3'` | Server-added keys make `job.request` fail `PrimerCheckRequest` validation; bumping the global version invalidates every existing check job for a change that cannot affect their results. The `typeof` guard is needed because `job_id.test.js:183` injects a stub check module [V]. |
| One `GENOTYPING_SET_INVALID` code with `reason` values | separate codes per failure mode | Matches the existing `GENOME_NOT_CHECKABLE {details.reasons}` precedent and keeps the UI's error mapping small. |
| 503 with `details.reason` for every upstream failure | a new 502 `VARIATION_SOURCE_INVALID` | 502 is outside the documented status set and outside the UI's retry handling, which treats 5xx as retryable but expects `retry_after_s`. |
| Warn-only Tm balance in v1 | moving a 5′ end by ±3 nt and rescoring | It would give the two pairs of a set different `expected` windows, which the submit-time validation forbids; deferred together with that rule. |
| `avoid_repeats` default false | defaulting to true with the AS zone exempt | Measured: on the example window masking removes a quarter of the template and one verified example set [V]. Recommended as a confirm-item rather than changed silently. |
| Rank by the §4.16 formula; `AS_TM_IMBALANCE`, `COMMON_TM_OUT_OF_RANGE` and `info` issues are not counted as issues (verification) | rank by Primer3 penalty with issues as a tie-break, which would have kept the former S2 `f970d5fa9d74` | Penalty-first puts a 30-mer with two ≥ 55 °C tailed hairpins (22.30) ahead of a 29-mer with one (19.63), exactly the failure the tail thresholds exist to catch; counting the two Tm issues on top of their own terms punished one imbalance twice. Every example is recomputed by `scratch-fix/score_sets.js`. |
| Score at most 8 candidates per orientation, reserve each candidate's worst case, caps 54 runs / 272 calls (verification) | scoring every returned pair (up to 40 kasp / 80 as_pcr runs and about 350 `ntthal` calls); `PRIMER_NUM_RETURN` 10 | Scoring everything overruns the 1.5 s target in the worst case and would make `DESIGN_BUDGET_EXHAUSTED` routine; NUM_RETURN 10 starves orientations whose best pairs fall to the zone and neighbour filters. The capped worst case measured 1.15 s. |
| Common-primer status from the common site's own alignment, best over both pairs' on-locus products (verification) | reading it from the REF pair's product, the ALT pair's product, or a product likelihood | Any product-based reading predicts `no_call` for AS-PCR sets, because one of the two products is `unlikely` by design; the ALT-pair reading gave `no_call` on 6 of the 11 real assemblies. |
| Core = left-aligned VCF span with its anchor base, extended by `shift`, ±1 bp; haplotypes = core ± K (verification) | zone ±1 bp; haplotypes = variant ± K | The zone ±1 runs one base past the first post-tract base (`TCAAAGT` at rs5413864115); both call alleles identically, so the definition that matches the pinned values was kept. The core is not the shift tract (next rows). |
| Copy `identity` and `aligned_length` from the longest anchor alignment (verification) | the minimum over anchors | The minimum reports the shortest window, dominated by the variant itself (pi180348 98.18 % over 165 columns against 98.44 % over 192). |
| Shift tract = the reference bases an indel slides through, anchor base excluded, used only by `in_shift_tract` (re-verification) | the anchor-inclusive span that also bounds the core | The anchor is not a repeat base. With it, the reverse ALT primer of every shiftable deletion (235 of 235 in 1:11080–61079) is `uncertain`, so every REF genome predicts `unknown` for reverse sets, and rs5413864115 S2 gains an issue (24.95). |
| Neighbour issues once per primer site: the allele-specific pair and the common primer (re-verification) | once per neighbour; once per oligo; once per set | Both allele-specific primers cover one genomic site, so per oligo counts it twice (tmp_1_11193_C_T S1 23.51). Per set, a second site under the common primer would cost nothing. Per neighbour, one dense site would dominate the score, and `DENSE_NEIGHBOURS` already reports density. |
| `ALT_PRIMER_SUBOPTIMAL` from each ordered oligo's own scoring run (re-verification) | also from the matched ALT primer's run when a deliberate mismatch is applied | The issue describes what is synthesized. The matched ALT primer is never ordered and its Tm is already floored through `matched_tm`; counting it would add a warn to every AS-PCR set whose matched ALT primer runs cool, and would reorder rs871475760 AS-PCR. |
| Exact decimal arithmetic and one half-away-from-zero rounding for every reported number (re-verification) | rounding binary floating-point results | Floats turn 57.102 − 56.427 = 0.675 into 0.67; over the 65 example candidates the shortcut changes 7 values. |

---

## 13. Critic findings resolution

Every blocker and major finding from the three critiques, plus the twelve items named in the brief. "Major" and "minor" are the critics' own severities; minors are listed where the resolution is not obvious from the majors.

| # | Finding (source) | Resolution | Section |
|---|---|---|---|
| 1 | **Brief 1 / assay major (all specs):** the allele caller names a third repeat length after the nearer allele | exact-core rule: `ref`/`alt` only on an exact match of the core — the anchored VCF span extended by the shift, ±1 bp (inserted bases included; §5.6); everything else is `other`, and `observed` is always reported; flank edits are confidence only. Verified on the two-A deletion, `AAAA`, the TGG repeat and a third SNP allele | §5.6 step 4, §7.3 |
| 2 | **Brief 2 / assay major (all specs):** amplification prediction ignores the common primer | the common primer is classified from its own site's alignment, best over the on-locus products of both pairs (row V4); distance 1 or a blocked site → `no_call` with `common_primer_3p_mismatch`, distance 2–3 → `weak`; `agrees` becomes `null`, never `true` | §5.7, §2.12 |
| 3 | **Brief 3 / assay major (reuse):** the relaxation ladder can emit sub-annealing primers | `PICK_ANYWAY` is never used in a design run; hard floors Tm ≥ 52 °C (all AS primers) and GC ≥ 15 % (derived primers); otherwise the orientation reports `no_sets` with `explain` | §4.6 |
| 4 | **Brief 4 / assay major (consumer):** the common-primer guard uses only the REF footprint and can overlap the AS primer | the zone includes both discriminating positions and both ALT anchor mappings; a `SEQUENCE_TARGET` keeps the common primer ≥ 10 nt away; a post-filter re-checks both haplotypes and the ALT `check_primers` run proves the common primer exists on the ALT haplotype | §3.5, §4.4 |
| 5 | **Brief 5 / assay minor (consumer, biology):** the deliberate-mismatch base comes from a class shortcut, and the position mixes sources | literal 32-cell Table 9.8.1, per-primer rows; verified values for both example variants at −2 and −3; default −2 with a sourced rationale, −3 labelled extrapolated, both raised as a confirm-item | §4.11, §11 item 2 |
| 6 | **Brief 6 / feasibility majors (all three specs):** Ensembl calls hold a design slot and lack limits | all Ensembl work runs before `semaphore.acquire` under a 4-way limiter, single flight, a 3-failure/60 s breaker, bounded LRU + negative caches, 8 s timeouts, byte caps and 10 kb chunking; outages degrade to `NEIGHBOURS_UNAVAILABLE` and never touch `/primers/design`; a slot-starvation integration test enforces it | §3.3, §4.2, §7.6 |
| 7 | **Brief 7 / contract + feasibility majors:** id patterns reject real ids; designs key on ids | pattern `^[A-Za-z0-9][A-Za-z0-9_.:,*-]{0,254}$` in swagger, handlers and UI state; record ids validated leniently and never dropped; designs key on the VCF `key`, with the picker designing from `vcf` | §2.4, §3.3, §3.4, §8.5 |
| 8 | **Brief 8 / feasibility major (biology) and minor (reuse):** per-variant FASTA/overlap fetches | one FASTA read per operation with pure normalization over it; one shared 10 kb chunk cache for listing, resolve, design neighbours and the UI track | §3.3, §3.6 |
| 9 | **Brief 9 / contract major (biology):** description-only response objects, `x-nullable` blindness, unused definitions | every response object has a `properties` block; `x-nullable` wherever null occurs; a **blocking** example test with an `x-nullable`-aware walker; no unreferenced `Primer*` definition | §2.6, §2.8, §2.12, §7.7 |
| 10 | **Brief 10 / contract + feasibility minors:** stored check requests contain server-added keys | `request.genotyping` is stored exactly as sent (values normalized only); orientation is re-derived; the version rides in `jobId`'s `algo` as `'2+g1'` through a `typeof`-guarded `algorithmVersionFor`; existing ids and `ALGORITHM_VERSION '2'` unchanged and pinned by test C1 | §5.3, §7.2 |
| 11 | **Brief 11 / feasibility minor (reuse):** indel ALT products break pair invariants | sets are their own response shape; `products.ref`/`products.alt` carry sizes, template spans, genomic blocks and `inserted_bases`; `pairs[]` is untouched | §2.8, §4.9 |
| 12 | **Brief 12 / feasibility minor (all):** the baseline suite is flaky | root-caused to a real race in `worker_main.js` `idle()` (handlers installed after the readiness log line; 148–210 µs window, reproduced 10/10 under an injected stall); fixed in M1; the gate compares against a recorded baseline rather than "503 pass" | §7.1, §10 M1 |
| 13 | **Contract major (reuse):** variant mode overloads the design request/response | dedicated endpoints; `PrimerDesignRequest`, `design.MODES`, the enum-sync test and the 22 existing fixtures untouched | §2, D1 |
| 14 | **Contract major and assay major (reuse):** set identity exists only by id convention; no submit validation | explicit `genotyping.sets[]`, orientation and allele side derived from the shared primer, and full submit-time validation with one error code and named reasons | §2.11, §5.2 |
| 15 | **Contract major (reuse):** browsing variants through the design endpoint can return `503 BUSY` | `GET /primers/variants` is a separate endpoint outside the semaphore with its own timeout and cache | §2.3, §4.2 |
| 16 | **Assay major (reuse):** under `auto`, an EMS target ignores every neighbour | natural neighbours are always detected and reported; for an EMS target they warn (high severity) and cost ranking points instead of blocking; only EMS neighbours are informational | §4.15 |
| 17 | **Assay major (reuse):** off-locus products ignored in the prediction | an amplifying off-locus product of the **other** allele's pair adds that dye (`both`) with `*_signal_off_locus`; a same-allele paralog only raises `off_locus_products` (corrected in row V14) | §5.7 |
| 18 | **Assay minor (all):** wrong-allele primers inside a shift tract are reported with false confidence | `in_shift_tract` at design time, from the anchor-excluded shift tract (row RV1), `uncertain` at check time, excluded from concordance counts | §3.5, §4.12, §5.7 |
| 19 | **Assay minor (consumer, biology, reuse):** copies double-counted or dropped at window edges | copies merge by **aligned variant coordinate**, not by amplicon bounds or 5′ site; the core boundary moves to the nearest aligned base rather than forcing a fallback | §5.6 steps 3, 5 |
| 20 | **Assay minor (all):** ortholog identity counts every gap column | gap-compressed identity for the ortholog test, with both numbers reported | §5.6 step 6, §2.12 |
| 21 | **Assay minor (consumer):** discrimination reported only as `mm_pos` | `terminal_mismatch_class` per allele-specific primer, with an info issue when both classes are weak | §4.12 |
| 22 | **Assay minor (reuse):** `distance_3p` for common primers measured from the 5′ end | defined and computed along the primer from its own 3′ end, with the forward/reverse mapping spelled out | §4.15 |
| 23 | **Assay minor (reuse):** tail-structure warnings fire on every metric at 47 °C | hairpins warn at 47 °C and are high at 55 °C; dimers warn at 47 °C; a warn issue adds 1 to the score and a high issue 3 (row V1) | §4.14 |
| 24 | **Assay minor (reuse):** `common_minus_as` computed but never used | warned outside −1…+3 °C and scored, and compared against nominal Tm values on both sides | §4.13, §4.16 |
| 25 | **Assay minor (consumer):** repeat masking default | kept `false` with the measured consequence, raised as a confirm-item; the AS-zone exemption is implemented either way | §4.5, §11 item 6 |
| 26 | **Contract minors:** unused definition, YAML fragments, `PrimerWarning.details`, 422→400 statuses, no 502, `consumes` on GETs, fixture routing, cost parity, generator presets, TS mode leak, matching by id | all adopted: every fence is at definitions level, `details` declared, request-vs-genome conflicts are 400, upstream failures are 503 with `reason`, both GETs declare `consumes`, fixtures routed by prefix and wrapper, parity pinned to real genome sizes, the generator picks the preset from the assay, `DesignMode` stays narrow while `DesignerMode` is applied consistently, and results match by sequence triple | §2.5–§2.15, §5.4, §7.7, §8.2, §8.6, §8.7 |
| 27 | **Feasibility minors:** code-reference offsets, pm2 recipe, Ensembl in tests, megablast budget, conditional `breakdown.genotyping`, fixture drift, scope | anchors re-verified in this document (with four corrections), a reversible worker-first pm2 recipe, a loopback fake Ensembl with five modes, a 30-fallback budget at 20 s with no retry, a conditional breakdown key, stubs that fail on unrecorded input, and the design-time paralog scan and Tm rebalancing deferred | §6, §6.4, §7.4, §7.6, §5.4, §5.6, §9.3 |
| V1 | **Verification (major):** the §4.16 formula did not reproduce the example scores, sets and keys (failed claims: `f970d5fa9d74` not the best forward set; scores 8.18, 22.30, 21.47, 24.81) | one formula stated exactly: which issues count (never the two Tm issues or `info`), unrounded terms, the product term at every kasp level, 2-decimal rounding. Every score, rank, key, pair id, order row, check fragment and test pin regenerated by `scratch-fix/score_sets.js`. New sets: rs871475760 S2 `1accc54c262d` (19.63), tmp_1_11193_C_T S2 `c791a4956fd3`, tmp_1_11502_C_CGT S1 `53942cb55348`; rs871475760 S1 now carries `ALT_PRIMER_SUBOPTIMAL` (9.18). §4.14 and §9.1 no longer claim that only high issues move a set | §4.16, §2.9–§2.11, §2.13, §4.14, §7.3, §7.6, §7.8, §8.9, §9.1 |
| V2 | **Verification (major):** explain strings, `pairs_returned` and `sets_found` were recorded at `PRIMER_NUM_RETURN` 10 while §4.7 prescribes 20 (failed claims: §2.9 and §2.10(a) attempts) | 20 kept; every attempt of §2.9–§2.10 re-recorded with exactly the §4.7 records (`scratch-fix/p3/`), e.g. rs871475760 forward `considered 1463, unacceptable product size 1443, ok 20`; attempts gain `not_scored` | §4.7, §2.8, §2.9, §2.10 |
| V3 | **Verification (major):** the process budget was smaller than the prescribed work, with no rule for which candidates are scored (failed claim: §4.1 step 11 "≤ 34 runs") | at most 8 candidates scored per orientation, in Primer3 order after the filters. Each candidate reserves 1 or 3 runs and 15 + 2 calls. Caps are 54 runs / 272 calls, the exact worst case over kasp and as_pcr, both orientations, three levels, tails and mismatch. Exhaustion gives `not_scored` plus `DESIGN_BUDGET_EXHAUSTED` in a `200`, never a `500`. The worst case measured 1.14–1.15 s with the real binaries | §3.1, §4.1, §4.8, §4.16, §2.8, §2.15, §7.3, §9.1 |
| V4 | **Verification (major):** the common primer's status had no defined source product, so AS-PCR genomes could all predict `no_call` | common status from the common site's own alignment (`isIgnored` / `isOverAmplifyingCap` / `isBlocked`, then `mm_pos`), best over the on-locus products of both pairs; deliberate mismatches never reach it; an ordered prediction table. KASP S1, S2 and the AS-PCR set give the intended results on the fake world and agree 11/11 on the real assemblies | §5.7, D10, §2.12, §7.3 |
| V5 | **Verification (major):** the core formula contradicted its pinned values, and the §2.13 haplotypes were not core ± K (failed claims: §5.6 core; §2.13 31-nt haplotype) | core = `[vcf.position − 1, vcf.position + len(ref) − 1 + shift + 1]`, the ALT core the same span on the ALT haplotype, haplotypes = core ± K; five variants recomputed (§5.6 table); third repeat lengths still `other`; §2.13 haplotypes now 33 nt | §5.6, §2.12, §2.13, §7.3 |
| V6 | **Verification (major):** swagger changes sat in M9 although M4, M6 and M7 need them over HTTP; several files had no milestone | each path, its definitions, controller handler, config keys and contract-test edits move to the milestone that first exercises them: M4 variants and genomes, M6 design, M7 check request, M8 results. The recorder goes to M5/M6, `blast.js` to M8, the fake Ensembl to M10 and §6.3. M9 keeps docs, the spec copy, example validation and the fixture-generator hand-off; §6.2 gains a milestone column | §10, §6.2, §6.3 |
| V7 | **Verification (minor):** D6 placed the overlapping common primer at 1:11086–11109 (failed claim) | 1:11108–11131 | §0 D6 |
| V8 | **Verification (minor):** the example `exempt` spans did not follow §4.5 (failed claim) | `[[365,37],[401,37]]` (§2.9) and `[[367,38],[400,38]]` (§2.10(b)) | §2.9, §2.10(b) |
| V9 | **Verification (minor):** examples omitted reportable neighbours (failed claims: §2.9 S1 common; §2.10(c) S1 warnings) | §2.9 S1 common lists EMS `tmp_1_11069_G_A` at 22; the new §2.10(c) S1 common primer carries `rs5413863452` at 6 with its `NEIGHBOUR_IN_PRIMER`; neighbour lists are sorted by distance | §2.9, §2.10(c), §2.10(d), §2.8 |
| V10 | **Verification (minor):** `Y` at 1:11318 contradicted the multi-allelic rule (failed claim) | a `*` allele is ignored when deciding whether a neighbour is biallelic; the string stays and no warning is added | §4.17 |
| V11 | **Verification (minor):** 6 primers over 3 genomes cost 95, not 94 (failed claim) | 95 | §5.4 |
| V12 | **Verification (minor):** internal contradictions (failed claims: D8 GC 20 %; `oligotm` in §2.1/D7; the 19.35 % ALT primer is S2's; section references) | D8 says GC ≥ 15 % for derived primers; §2.1 and D7 say Primer3 `check_primers` plus `ntthal`; §2.10(a) prose names S2; references to §4.15/§4.16 corrected; GC and Tm rounding defined in §2.1 | §0, §2.1, §2.8, §2.10(a), §2.10(d), §3.7, §4.8 |
| V13 | **Verification (minor):** a core end aligned to a gap was undefined | the end moves outward to the nearest aligned genome base, and with no such base the anchor is `missing`; tests added (1:11110 deleted → `other`, 1:11108 deleted → `ref`) | §5.6 step 3, §7.3 |
| V14 | **Verification (minor):** any amplifying off-locus product turned the prediction into `both`, even a same-allele paralog | only a product of the other allele's pair adds its dye; a same-allele paralog raises `off_locus_products` and its reason only | §5.7, D10, §7.3 |
| V15 | **Verification (minor):** a scoring-run `PRIMER_ERROR` would fail the whole request through the design-run mapping | `400 PRIMER3_INPUT_ERROR` applies to design runs only; in scoring runs a `PRIMER_ERROR` increments `rejected.alt_scoring_failed` | §4.7, §2.8 |
| V16 | **Verification (minor):** the TS types had placeholders and undefined names | every interface written out from the §2.6–§2.12 definitions, plus a type-parity test | §8.2, §8.8 |
| V17 | **Verification (minor):** the client capped checks at 13 primers whatever the genome count | gate on the mirrored `over_limit` plus 20 primers and 10 pairs; 13 stays only as the full-panel `check.request` cap | §8.5, §8.6 |
| V18 | **Verification (minor):** wrong anchors (failed claims: `checkParams` already exported; `PangenomeMatrix.tsx:49-79`) | `design.js` exports only `parseParams`; keyboard handler at `:307`; `src/components/designChecks.ts`; `examples/playground/mockClient.ts` and `App.tsx` | §6.2, §8.4, §8.5, §8.8 |
| V19 | **Verification (minor):** M9 edited gramene-primers, which the UI session owns | the API session never edits gramene-primers; it hand-writes provisional request fixtures in each endpoint milestone, and the UI session regenerates them | §7.7, §6.3, §10 |
| V20 | **Verification (minor):** the dev block lacked Node 24.14.1 and the fake-Ensembl start, mode and restart commands | `nvm use 24.14.1`, the fake server with data and control ports, mode switching and the API restart are spelled out | §6.4, §7.6 |
| V21 | **Verification (minor):** the shared fetch carried the first caller's deadline, contradicting single flight | the shared fetch uses only `AbortSignal.timeout`; each caller races it against its own deadline | §3.3 |
| V22 | **Verification (minor):** §2.13 `aligned_length` 205 and identity 98.2 matched no defined alignment | copy values come from the longest anchor alignment: 192 / 100 % (reference, v5, pi329250), 192 / 98.44 % and 193 / 98.45 % (pi180348) | §5.6 step 5, §2.12, §2.13 |
| RV1 | **Re-verification (major):** §5.6 defined the shift tract to include the VCF anchor base. That gave rs5413864115 S2 a second `SHIFT_TRACT_DISCRIMINATION` (24.95, not the pinned 23.95) and made §5.7 predict `unknown` for reverse sets on every REF genome | two terms. The shift tract (§3.5) excludes the anchor base and serves only `in_shift_tract`, now defined per primer (§4.12) and used by the §5.7 status table. The core (§5.6) is the anchored VCF span extended by the shift, ±1 bp. S2 stays 23.95; REF genomes predict `ref` in both orientations at rs5413864115 and rs5413863413; the two-A deletion and `AAAA` stay `other` | §3.5, §4.12, §5.6, §5.7, §2.8, §2.12, §2.15, §7.3, §9.1, §12, D9 |
| RV2 | **Re-verification (minor, row V1):** three scoring rules were implicit: neighbour-issue grouping, whether the matched ALT primer's run raises `ALT_PRIMER_SUBOPTIMAL` under a deliberate mismatch, and the penalty of the candidate's own Primer3 pair rather than pair 0 | §4.16 rules S1–S7 with a per-code table of how many issues a set can carry; §4.15 "one issue per primer site"; §2.8 `primer3_penalty` and §2.15 wording. All eight example scores and ranks recomputed from the rules, unchanged | §4.15, §4.16, §2.8, §2.15, §7.3 |
| RV3 | **Re-verification (minor):** tmp_1_11193_C_T S1 `common_minus_as` regressed from 0.68 to 0.67, a floating-point artefact | one rounding rule for every reported number: exact decimal arithmetic from unrounded inputs, then one half-away-from-zero rounding; 0.68 restored and pinned | §2.1, §2.8, §2.10(a), §4.13, §4.16 S7, §7.3 |
| RV4 | **Re-verification (minor, row V6):** M6's `contract.test.js` needed M7's `BODY_KEYS` change, and the `:136-179` bodies, `controller_wrap.test.js` and `job_id.test.js` had no milestone | sync tests split (design half `:317-364` in M6, check half `:366-387` in M7); documented bodies split between M6 and M7; `controller_wrap.test.js` in M4 and M6, `job_id.test.js` in M7; §6.3 and §7.7 gain milestone columns | §10, §6.3, §7.3, §7.7 |
| RV5 | **Re-verification (minor, row V12):** stale references in F12, §2.6 and the §2.9 introduction | §4.17, §7.7 and §7.4/§6.3 | §0.1, §2.6, §2.9 |
| RV6 | **Re-verification (minor):** `common_primer.likelihood` had no defined source product | every primer call reports `likelihood`, `mm_pos` and `residual_mm_pos` from the product that supplied its status, ties broken by likelihood, then product order; §2.13 reproduces unchanged | §5.7, §2.12, §7.3 |
| RV7 | **Re-verification (minor):** the §7.3 budget pin omitted `thermo_calls` and `max_thermo_calls` | deep-equal details with `thermo_calls: 87, max_thermo_calls: 272`; §2.15 says the counters are those at the end of the request | §7.3, §2.15 |
| RV8 | **Re-verification (minor):** flow A4 and the §2.7 description gave `avoid_repeats` a default of true | both say false, a default to confirm (§11 item 6) | §1.1, §2.7 |
| RV9 | **Re-verification (minor):** the §4.16 column headed "level" held `5 × level` | columns named after the rules (`S2 5 × level`, …), plus an exact-sum column | §4.16 |

---

## Appendix A. Verification log

Scratch directory: `genotyping-design/scratch-synth/`. Everything was read-only: no repository file was modified, no pm2 process or port was touched, and no check job was submitted. The dev API on port 50111 was used only for `template_only` reads and `GET /primers/genomes`.

| Script / artefact | What it did | Key outputs used |
|---|---|---|
| `tpl_10650_11960.json`, `tpl_11150_11750.json`, `tpl_13600_13820.json` | `POST /primers/design` `template_only` on the 50111 dev API | reference bases for every example; 11109 = C, 11282–11286 = CAAAG, 11500–11503 = GCCA, 13730–13755 = TGCCATTGGTGGTGGTGGTGGTAGAA |
| `tpl_masked.json` | the same window with `avoid_repeats: true` | mask `[[1,200]]`, `blast_depth`, `masked_fraction 0.2497`, warning `BLAST_DEPTH_MASK` (§4.5, §11 item 6) |
| `verify_seqs.py` | checked every example primer against the genome, derived the ALT primers, recomputed the set keys and the submission string | all 16 primer/coordinate assertions pass; keys `f9df650ad116`, `f970d5fa9d74`, `7f9af6b1c938` reproduce |
| `little.py` | implemented the literal Table 9.8.1 and computed the mismatch bases | the §4.11 table; identical to an independent derivation of the same eight substitutions |
| `p3run.py` (runs A–I) | `primer3_core` 2.6.1 with the prescribed records | F2–F4: the overlap failure and its fix, the guard in both orientations and for the deletion, the level-2 floor at 52.87 °C |
| `p3check.py`, the 11193 scoring runs | `check_primers` on REF/ALT/mismatch haplotypes | ALT and mismatch values, products 65 / 87 / 117 / 116 / 126, `PROBLEMS` strings, and equality with the design run |
| the `ntthal`/`oligotm` batch | Tm, hairpins, self and cross dimers, tailed oligos, mismatch duplexes | every thermodynamic value in §2.9, §2.10 and §4.14 |
| `caller.py` | the exact-core rule against distance-based calling on real repeat sequence | F5: `other` for a two-A deletion, `AAAA`, two deleted TGG units and one inserted unit |
| `zone.py` | discriminating positions, ALT anchor mappings and zones for five variants | F6, including rs5413863413 `[13735, 13752]` |
| the id-pattern check over `ov_50k.json` | 3,978 real Ensembl ids | F10: 88 with `,`, 46 with `*`, 11 over 128 characters, longest 255; the adopted pattern accepts all |
| live Ensembl overlap calls (four windows plus two spot checks) | neighbours, EMS sources, both insertion conventions, normalized key counts | F11, the `neighbours` counts 57/45/12, 54/45/9, 53/45/8, 59/49/10 and the dense counts 0, 1, 4, 5 |
| the submission-string regeneration | genome plus live neighbours, IUPAC coding | F12: all four strings reproduce character for character |
| the `cost.js` run with real genome sizes | the server's own estimator | F9: 2,690 / 5,356 / 5,800 / 6,245 / 95 |
| the worktree audit (read-only) | ~90 file:line anchors, the flaky test and the baseline | §6 anchors and their four corrections; the `idle()` race, reproduced 10/10 under an injected stall; baseline 503 / 487 / 16 twice |

Values in this log that the verification fixes replaced are kept as a record of what was run. They are the former S2 key `f970d5fa9d74`, the former §2.10(a) S2 and §2.10(c) S1 sets with their products 116 and 117, and every explain string recorded at `PRIMER_NUM_RETURN` 10. Appendix B lists the current values and the scripts that reproduce them.

**Not verified, and marked [U] or [I] where used:** the megablast fallback on pi536008; the wet-lab performance of the presets and of the deliberate mismatch at −3; the genotyping CPU coefficient end to end; the `specificity`, `pangenome`, `primers` and `timings_ms` blocks of the §2.13 example (illustrative, to be replaced by a real 50112 job before the fixtures are handed over); and the behaviour of the swagger router for path ids containing `,` and `*`, which §7.6 tests explicitly.

---

## Appendix B. Verification fixes

Scratch directory: `genotyping-design/scratch-fix/`, with the exact re-check commands in `scratch-fix/RECHECK.md`. The rules of Appendix A held:
- no repository file was edited, no pm2 process touched, and no check job submitted;
- the 50111 API was used only for `template_only` reads of 1:10500–12100 and 1:13500–13900, byte-identical to the verifier's copies;
- the unedited spec is kept as `scratch-fix/spec.before-fix.md`.

Row numbers match §13.

| # | Gap or failed claim (`verify.json`) | What changed | Proof: script → result |
|---|---|---|---|
| V1 | §4.16 formula vs examples; failed claims on `f970d5fa9d74`, 8.18, 22.30, 21.47 / 24.81 | exact formula. New sets: rs871475760 S2 `1accc54c262d`, tmp_1_11193_C_T S2 `c791a4956fd3`, tmp_1_11502_C_CGT S1 `53942cb55348`. Scores, keys, pair ids, order rows, check fragments and pins regenerated | `score_sets.js`, the reference implementation of §4.3–§4.16 → the eight scores of the §4.16 table. `emit_examples.js` → the §2.9, §2.10(a), §2.10(c) and §2.11 fragments, spliced verbatim. `check_examples.js` → every example value, key and pin |
| V2 | `PRIMER_NUM_RETURN` 20 vs the recorded explain strings; failed claims §2.9 and §2.10(a) | all attempts re-recorded at 20; `not_scored` added | `score_sets.js` writes the exact records and outputs to `p3/`; `check_examples.js` compares every attempt |
| V3 | process budget; failed claim §4.1 step 11 | cap of 8 per orientation, reservation rule, caps 54 / 272, a warning instead of an error | `timing.js` → 54 runs + 272 calls in 1,141 and 1,152 ms. `CAP_PRIMER3_RUNS=12 node score_sets.js rs871475760_kasp` → exit 0, `DESIGN_BUDGET_EXHAUSTED`, reverse `not_scored 18` + `sets 2` |
| V4 | common-primer status | own-site status, ordered prediction table, AS-PCR tests | `predict.js`: part A → all fake-world expectations for KASP S1, KASP S2 and AS-PCR; part B → 11/11 agree for all three sets, where the old ALT-pair reading gives `no_call` in 6/11 for AS-PCR; 0 failures |
| V5 | core definition; failed claims §5.6 core and §2.13 haplotypes | the anchored VCF span extended by the shift, ±1 bp (called "shift tract ±1" at the time; Appendix C, C1 separates the two terms); haplotypes core ± K | `cores.py` → the five cores and haplotypes of §5.6; `TCAG`, `TCAAAG`→`TCAAAAG`, 3 and 6 TGG units, `GTGT` and `G` all `other`; 0 failures |
| V6 | milestones vs swagger | swagger, handlers, config and contract edits in M4/M6/M7/M8; M9 keeps docs, the spec copy, example validation and the generator hand-off | review of §10 against the milestone column of §6.2 (no script) |
| V7 | D6 coordinates | 1:11108–11131 | reverse complement of 1:11108–11131 in `tpl_10500_12100.json` = `TCTTTGACTAGCGAGAAATTCAGA` (checked in the final re-verification below) |
| V8 | `exempt` spans | computed from §4.5 | `check_examples.js` recomputes `exempt` for §2.9 and §2.10(b) |
| V9 | omitted neighbours | EMS neighbour added to §2.9 S1; the new §2.10(c) S1 carries its common-primer warning | `score_sets.js` over the verifier's live-Ensembl overlap `ov_10600_11960.json`; `check_examples.js` |
| V10 | `Y` at 1:11318 | `*` ignored when testing biallelic | rule decision; the verifier's `neighbours.py` found 11318 to be the only differing base |
| V11 | cost 94 → 95 | 95 | `scratch-verify/cost_verify.js` → 95 |
| V12 | contradictions (D8, D7, §2.1, §2.10(a), references) | text corrected; GC and Tm rounding defined | `check_examples.js` applies the rounding to every Tm and GC; text review |
| V13 | core end on a gap | outward-move rule | `cores.py` → 1:11110 deleted `other` (`TCG`); 1:11108 and 1:11281 deleted `ref`; no outward base `missing` |
| V14 | off-locus rule | other-allele dye only | `predict.js` part A → ALT paralog `both`; REF paralog `ref` with `off_locus_products 1` |
| V15 | `PRIMER_ERROR` mapping | design runs only | rule text; the verifier reproduced the error with `scratch-synth/chk_NEG_alt_on_refhap.in.txt` |
| V16 | TS types | written out from the §2.6–§2.12 YAML | no script; §8.8 specifies the type-parity test that automates it |
| V17 | client primer limit | `over_limit` gate | `scratch-verify/cost_verify.js` → 6 primers over 3 genomes 95 CPU-s; 14 primers over the full panel 6,245 |
| V18 | code anchors | corrected | the verifier's read-only anchor checks of gramene-swagger-primers and gramene-primers |
| V19 | gramene-primers ownership | UI session owns it; provisional fixtures | text |
| V20 | dev instance commands | `nvm use`, fake Ensembl, restart | text |
| V21 | single flight vs deadline | shared timeout only | text |
| V22 | §2.13 `aligned_length` / identity | longest anchor alignment | `predict.js` part C → 192 / 100 %, 192 / 98.44 %, 193 / 98.45 %; gap-compressed identical (one 1-column gap run) |

**Final re-verification** after all edits:
- `check_examples.js`, `cores.py` and `predict.js` as above;
- `scratch-verify/swagger_verify.js` → re-parses every JSON and YAML block, merges the fragments into a copy of `swagger.yaml`, runs sway, and validates the example requests and responses;
- `scratch-synth/validate3.js` → the `x-nullable`-aware walker over the §2.9 and §2.13 examples.

The results are recorded in `scratch-fix/final_checks.log`.

---

## Appendix C. Re-verification fixes

Scratch directory: `genotyping-design/scratch-fix2/`. The exact re-check commands and their expected outputs are in `scratch-fix2/RECHECK2.md`. The rules of Appendices A and B held:
- no repository file was edited, no pm2 process touched and no check job submitted;
- the 50111 API was used only for `template_only` reads of 1:10500–12100, 1:13300–14200, 1:13400–14000, 1:10000–59999 and 1:59000–62100 (the first and third are byte-identical to the re-verifier's copies, and all five agree where they overlap);
- the unedited spec is kept as `scratch-fix2/spec.before-fix2.md`.

Rows C1–C9 are §13 rows RV1–RV9.

| # | Issue (`reverify.json`) | What changed | Proof: script → result |
|---|---|---|---|
| C1 | the shift tract included the anchor base: S2 scored 24.95 and reverse sets predicted `unknown` on REF genomes | §3.5 `tract(v)` without the anchor; §4.12 per-primer rule; §5.6 core named separately; §5.7 status row uses `in_shift_tract`; D9, §2.8, §2.12, §2.15, §9.1, §12 | `tract56.py` → the five §5.6 rows with tracts and flags; zone = tract ± 1 and the deletion/insertion geometry for all 432 shiftable indels of 1:11080–61079 (the anchor-inclusive wording would flag the reverse ALT primer of all 235 deletions); the two-A deletion, `AAAA`, two TGG units deleted and one inserted are `other`; 0 failures. `predict_tract.js` → REF `ref` and ALT `unknown` for rs5413864115 S1/S2 and both rs5413863413 sets, with both reverse sets `unknown` on REF under the old wording; the rs5413863413 primers match the genome; 0 failures. `score2.js` → S2 `cb3ef66afd37` 23.95 (exact 23.953151) with one `SHIFT_TRACT_DISCRIMINATION`; the re-verifier's `rescore.js` with `TRACT=noanchor` → 23.95 |
| C2 | scoring rules not written down | §4.16 S1–S7 and the per-code issue table; §4.15 grouping; §2.8 `primer3_penalty`; §2.15 | `score2.js`, written from S1–S7 → all 8 scores, keys, orientations, levels and counted issues; the eight forward rs871475760 candidates; budgets 18/123, 20/150, 20/141, 9/66, 29/4; 29 checks, 0 failures. Ranks unchanged, so no set id, pair id, order row, submission string, check fragment or pin moved |
| C3 | `common_minus_as` 0.67 | §2.1 rounding rule; 0.68 in §2.10(a); §4.13, §4.16 S7, §7.3 | `score2.js` → exact 0.675000 → 0.68, and the list of the 7 values that float rounding would change over 65 candidates. `fx/score_sets.js`, `fx/check_examples.js` and `fx/emit_examples.js` use the same decimal-safe rounding |
| C4 | M6 needed M7; three test edits without a milestone | §10 M4/M6/M7; milestone columns in §6.3 and §7.7; §7.3 rows for `job_id.test.js` and `controller_wrap.test.js` | anchors read from the worktree (`contract.test.js:317,342,366`, `job_id.test.js:181`, `controller_wrap.test.js:103,239,263`); `rv/staged_swagger.js` → M4, M6, M7 and M8 add no error or warning |
| C5 | stale cross-references | F12 → §4.17; §2.6 → §7.7; §2.9 → §7.4/§6.3 | text |
| C6 | `common_primer.likelihood` undefined | §5.7 product rule; §2.12 descriptions; §7.3 row | `likelihood213.js` → every `ref_primer`, `alt_primer` and `common_primer` object of §2.13 (reference and three genomes, S1 and S2) reproduces under the rule; 0 failures |
| C7 | budget pin incomplete | `thermo_calls: 87`, `max_thermo_calls: 272` | `CAP_P3=12 node score2.js` and `CAP_PRIMER3_RUNS=12 node fx/score_sets.js` → identical details |
| C8 | `avoid_repeats` default | false in A4 and §2.7 | text |
| C9 | §4.16 column name | renamed, exact sums added | `rv/compare_spec.js` reads the renamed column |

**Re-run of the re-verifier's own scripts.** Copies are in `scratch-fix2/rv/` with their inputs linked, so the originals and their logs stay untouched; `diff` against `scratch-reverify/` shows every change. Expected values changed only where the spec changed on purpose:
- `rescore.js` runs with `TRACT=noanchor`, the reading the spec now states;
- `compare_spec.js` runs with `MODE=group_noanchor_ordered_pairi`, and its §4.16 column-5 expectation is `5 × relaxation_level`, following the renamed column;
- `cores56.py`, `predict57.js`, `spot7.py`, `budget.py`, `staged_swagger.js` and `swagger_verify.js` run unchanged.

The fix agent's checks are re-run from copies in `scratch-fix2/fx/`. `score_sets.js`, `check_examples.js` and `emit_examples.js` use the §2.1 decimal-safe rounding, so the regenerated §2.10(a) fragment differs from `scratch-fix/examples/` only in `common_minus_as` 0.68. `predict.js` reads the megablast table by absolute path, because the copy runs one directory deeper.

---

## Implementation deviations (M3–M7)

This section was added when the spec was copied into `docs/` (milestone M9a). It lists where the implementation on the
`genotyping` branch departs from the spec above, milestone by milestone, as recorded when each milestone was built.
`docs/primer_design_api.md` documents the implemented behaviour.

### M3 — variant normalization (`variation/normalize.js`)

| Topic | The spec said | What was built |
| --- | --- | --- |
| Exports | §6.1 lists the exports of `variation/normalize.js`. | Also exports `tractOf`, `sequenceWindow`, `SequenceWindowError`, `inWindow` and `DEFAULTS`. |
| `recordsToEntries` | §3.4 builds entries with the requested id ordered first. | Returns `{entries, skipped: {count, reasons}}` and never adds `requested_id` or `submission_sequence`: the callers do. `variation/index.js` puts `requested_id` on the design variant, and the lookup reports it at the top level only; the design adds `submission_sequence`. |
| Unshiftable repeats | §2.6 declares `zone` and `discriminating` as plain objects. | An entry whose shift exceeds `max_shift`, or whose repeat reaches the region end, gets issue `REPEAT_TOO_LONG` with `zone` and `discriminating` set to `null`. Swagger marks both `x-nullable`. |
| Id order | §3.4: the requested id first, then `rs` ids, then the rest alphabetically. | `rs` ids are ordered by number. |
| Manual alleles | §2.6: `alleles` are the site's alleles as reported. | A manual entry reports its alleles in minimal form. |

### M4 — Ensembl client, variant resolution, variants endpoints, genomes flags

| Topic | The spec said | What was built |
| --- | --- | --- |
| Chunk end | §3.3: `chunkEnd = chunkStart + chunk_bp − 1`. | `overlapChunk` also takes `regionLength` and clamps the last chunk to the region end. |
| Record ids | §3.3: ids are validated leniently, so no record is dropped because of its id. | An id that fails even the lenient rule keeps its record, with `id: null`. |
| Limiter queue | §3.3: a waiter that does not start within `queue_wait_ms` fails with `queue_full`. | The queue is also capped at 256 waiters; overflow fails at once with `queue_full`. |
| Redirects | §3.3: `redirect: 'error'`; the status mapping does not list 3xx. | A 3xx answer maps to reason `invalid_response`. |
| Design by id with variation disabled | §6.4: `PRIMERS_VARIATION_ENABLED=0` disables the variants endpoints (`503 FEATURE_DISABLED`). | A design by id on a genome listed in `variation.species` also answers `503 FEATURE_DISABLED`. |
| Lookup entries | §2.4: one entry per designable alt. | A lookup returns every entry of the id, the non-designable ones included, with their issues. |
| Lookup-only ids | §3.7: id lookups take `source` from the overlap records at the same site. | An id the overlap does not list falls back to its lookup mapping, with `source: null`. |
| Nullable fields | §2.6: `zone`, `discriminating` and `PrimerVariantRecord.source` are not `x-nullable`. | Swagger marks all three `x-nullable`. |

### M5 — design core

| Topic | The spec said | What was built |
| --- | --- | --- |
| Blocker distances | §2.8: `distance_from_3p` is `null` in orientation blockers, where no primer exists yet. | It is filled in for blockers, as in the §2.10(c) example (`distance_from_3p: 4`). |
| Param rules | §2.7 rule 5: the cross-field param rules, `400 INVALID_PARAMS {param}`. | The rules are checked at every ladder level; a failure at a relaxed level also carries `details.level`. |

### M6 — sets, scoring and the design endpoint

| Topic | The spec said | What was built |
| --- | --- | --- |
| `engine.thermo` | §2.8: the `ntthal` version string. | `"ntthal "` followed by the `primer3_core` version, because `ntthal` prints no version. |
| Mismatch-set thermo | §2.8: `thermo.ref_common` and `alt_common` of a set. | For deliberate-mismatch sets they come from the mismatch runs. |
| Tailed-pairing issues | §2.15: `TAILED_STRUCTURE` details are `{oligo, metric, value, severity}`. | An issue about a tailed cross-dimer carries `details {pair, metric, value}`; the matching warning adds `severity`. |
| `AS_TM_IMBALANCE` high threshold | §2.15: high above 2.0 °C. | 2.0 is a constant (`sets.js AS_TM_DIFF_HIGH`) with no config key. |
| `ntthal` pool | §4.1 and §6.1: a pool of 4. | The pool of 4 (`genotyping.thermo_concurrency`) is per request, not process-wide. |
| Module layout | §6.1: the neighbour and pair-screen helpers in `genotyping/design.js`. | They live in `genotyping/sets.js`; `design.js` re-exports them. |

### M7 — check, submit side

| Topic | The spec said | What was built |
| --- | --- | --- |
| Rule order | §5.2 part B lists `too_many_edits` before `alleles_swapped`. | `alleles_swapped` is checked before `too_many_edits` (and again after), because the literal order contradicts the §7.3 expectations. |
| `common_in_zone` | §5.2: the common primer's footprint must not intersect `[zone.start − 1, zone.end + 1]`. | It also rejects a common primer that lies on the allele-specific side of the variant. |
| `expected_mismatch` | §5.2: the REF primer's 5′ end must equal `expected.start` (forward) or `expected.end` (reverse). | It also fires when `expected.region` differs from the variant's region. |
| Extra shape rules | §2.11: the `PrimerCheckGenotyping` definition. | Two more `400 INVALID_REQUEST` rules: REF must differ from ALT, and set ids must be unique. |
| Region edge | §5.2 part B: `400 REF_MISMATCH` or `VARIANT_TOO_REPETITIVE` from the positional checks. | Reaching the region edge while left-aligning gives `400 REGION_OUT_OF_BOUNDS`. |
| Beyond the read | §2.14: `VARIANT_TOO_REPETITIVE {region, position, shift, max}`. | Sliding past the ±1,700 bp read gives `VARIANT_TOO_REPETITIVE` with `shift: null`. |
| Where the prepared sets live | Not specified. | `check/genotype.js prepare()`'s output is stored in the job document at `resolved.genotyping`. |

### Found while documenting (M9a)

| Topic | The spec said | What was built |
| --- | --- | --- |
| Manual design with variation disabled | §3.8: `NO_VARIATION_DATA` for a design on a genome without variation data. | With `PRIMERS_VARIATION_ENABLED=0`, a manual design on `sorghum_bicolor` also gets `neighbours.data: "none"` and warning `NO_VARIATION_DATA`, whose message says the genome has no known-variant data. |
| `PrimerVariant.requested_id` | §2.6: "lookup and design responses only". | Only design responses carry it on the variant; the lookup response has it at the top level. Swagger describes the implemented behaviour. |
| `THERMO_UNAVAILABLE` | §2.14: `{binary, retry_after_s}`. | `retry_after_s` is 60 for a missing or non-executable binary and 5 when the process could not be started. |
| Lookup errors in swagger | §2.5: `GET /primers/variants/{variant_id}` 422 lists `NO_VARIATION_DATA`, `AMBIGUOUS_VARIANT_MAPPING`, `VARIANT_NOT_ON_ASSEMBLY`. | The handler can also answer `422 NO_SEQUENCE` and `AMBIGUOUS_ASSEMBLY`, and `503 MONGO_UNAVAILABLE`, which the swagger descriptions of that operation do not list. |

### M8 — check, worker side

<!-- M8 deviations: to be added by M9b once the worker-side allele caller (results.genotyping, §2.12–§2.13, §5.6–§5.9) lands. -->

To be completed when the worker-side caller lands (M8).
