# Genotyping front-end brief

Hand-off from the API session (milestone M12) to the UI session that builds KASP / allele-specific PCR into
**gramene-primers** and switches it on in **gramene-search**. It describes what was *built*, not what the spec asked for.

- **API worktree:** `/usr/local/gramene/subsites/sorghum/v11/gramene-swagger-genotyping`, branch `genotyping`, HEAD `141953f`.
  Unit suite: `node --test "test/primers/unit/**/*.test.js"` → **731 tests, 711 pass, 0 fail, 20 skipped** (Node 24.14.1).
- **Base path:** `/sorghum_v11` (production `https://data.sorghumbase.org/sorghum_v11`). The genotyping endpoints exist only
  where this branch is deployed; the pm2 repoint of the dev apps is a separate step that needs the user's go-ahead, so
  **confirm the live base URL with the user before building against it**. gramene-search enables the tab from
  `PRIMERS_API` (`src/demo.js:167`).
- **UI session owns** `gramene-primers` and `gramene-search` entirely. The API session wrote no UI code.

### Read order and authority

| Source | Use it for | Authority |
| --- | --- | --- |
| [`docs/primer_design_api.md`](primer_design_api.md), section [Genotyping primers](primer_design_api.md#genotyping-primers-kasp--allele-specific-pcr) | the implemented contract | **Wins over everything else** |
| `test/primers/fixtures/docs/*.json` | real request/response pairs the UI must render | Wins over prose |
| `docs/genotyping_design_spec.md` §8 | the original front-end plan: types, component list, state shape | Starting point only; its "Implementation deviations" section lists where it is wrong |
| `/home/olson/.claude/plans/i-want-to-add-radiant-lerdorf.md` | why the feature exists, decisions, risks | Background |

Do not copy the whole spec into code comments; cite doc anchors instead.

---

## 1. What exists server-side

| Endpoint | One-line semantics | New? |
| --- | --- | --- |
| `POST /primers/design` | PCR/qPCR design on one template | **Unchanged** — same request, same response, same warnings. No UI change required |
| `GET /primers/genomes?system_name=` | same-species genomes and capabilities | **Additive:** top-level `variation {available, source, release}` and `genomes[].has_variation`. Everything else byte-identical ([variation fields](primer_design_api.md#get-primersgenomes-variation-fields)) |
| `GET /primers/variants?system_name=&region=&start=&end=` | known Ensembl variants in a window, left-aligned, merged, REF-checked ([docs](primer_design_api.md#get-primersvariants)) | New |
| `GET /primers/variants/{variant_id}?system_name=` | one Ensembl id → one entry per alternative allele, plus `synonyms` ([docs](primer_design_api.md#get-primersvariantsvariant_id)) | New |
| `POST /primers/genotyping/design` | synchronous (≤ 45 s): ranked KASP / AS-PCR **sets** for one variant, with order rows and a ready check request ([docs](primer_design_api.md#post-primersgenotypingdesign)) | New |
| `POST /primers/check` | queue a specificity / pan-genome check | **Additive:** optional `genotyping {variant, sets}` block ([docs](primer_design_api.md#post-primerscheck-the-genotyping-block)) |
| `GET /primers/check/{job_id}` | poll a check | **Additive:** `results.genotyping` for genotyping jobs ([docs](primer_design_api.md#results-resultsgenotyping)) |

Facts that shape the UI:

- **No auth, `Cache-Control: no-store`, CORS `*`.** `Retry-After` is exposed to browsers.
- **The browser never calls Ensembl.** Variants always come through this API. Only `sorghum_bicolor` has variation data
  (1 of 120 genomes in the real `genomes` capture).
- **Identity is the VCF key** `region:position:REF:ALT` (`1:11109:C:A`). Ensembl ids are for lookup and display only:
  design from `variant.vcf`, never from the id, because real ids contain `,` and `*` and run to 255 characters.
- **Existing check job ids are unchanged.** A request carrying `genotyping` hashes with algorithm `"2+g1"`; everything else
  still hashes with `"2"`.
- **Jobs without a `genotyping` request have no `results.genotyping` key at all** (not `null`).
- Design is synchronous but shares the 4-slot Primer3 semaphore with `/primers/design` (`503 BUSY` after a 10 s wait).
  All Ensembl work happens *before* the slot is taken, so a slow Ensembl never starves ordinary designs.

---

## 2. Client methods, types and files to change (gramene-primers)

### 2.1 Client methods (`src/client.ts`, `PrimersClient` in `src/types.ts:638-650`)

```ts
listVariants(q: VariantListQuery, o?: RequestOptions): Promise<VariantListResponse>;                    // GET /primers/variants
getVariant(id: string, q: { system_name: string }, o?: RequestOptions): Promise<VariantLookupResponse>; // GET /primers/variants/{encodeURIComponent(id)}
designGenotyping(req: GenotypingDesignRequest, o?: RequestOptions): Promise<GenotypingDesignResponse>;  // POST /primers/genotyping/design
```

| Rule | Why |
| --- | --- |
| All three are **optional** members of `PrimersClient` | existing hosts and test doubles (including `examples/playground/mockClient.ts`) still compile; disable the Genotyping tab when `designGenotyping` is missing |
| `designGenotyping` uses the **design** timeout (60 s), the others `other` (15 s) | the server deadline is 45 s |
| `listVariants` sends `types` as CSV, sends `include_ems` only when `false`, memoizes per `system_name\|region\|start\|end\|types\|include_ems\|limit`, evicts on error, and one caller's abort must not cancel the shared promise | mirror `listGenomes` (`client.ts:209-224`, `withCallerSignal` at `:238`) |
| Response guards mirror `client.ts:174-181` | `designGenotyping`: `template` must be an object, default `sets`/`warnings` to `[]`. `listVariants`/`getVariant`: `variants` must be an array |
| Errors: reuse `toApiError` unchanged | handler bodies are `{message, code, details}`; `details.retry_after_s` already becomes `retryAfterMs` (`client.ts:64-65`) |

### 2.2 File-by-file

| File | Add / change |
| --- | --- |
| `src/types.ts` | the block in §2.3 below. `CheckRequest` (`:297`) gains `genotyping?`; `CheckResults` (`:546`) gains `genotyping?`; `GenomeEntry` (`:241`) gains `has_variation?`; `GenomesResponse` (`:255`) gains `variation?`; `PrimersClient` (`:638`) gains the three optional methods; `PrimerDesignerState` gains `genotyping?`; `PrimerDesignerProps.modes`/`defaultMode` (`:711-712`) and `PrimerDesignerFeatures` (`:697`, add `genotyping?`) widen to `DesignerMode` |
| `src/client.ts` | the three methods, a `variantsCache` beside `genomesCache` |
| `src/request.ts` | `buildGenotypingRequest(state, ctx)`, `buildGenotypingCheckRequest(...)` (§4), `GENOTYPING_CHECK_LIMITS = {maxSets: 5, maxPairs: 10, maxUniquePrimers: 20}`; `CheckRequestErrorCode` gains `TOO_MANY_SETS` and `OVER_CPU_LIMIT`. `CHECK_LIMITS` (`:22`) stays as it is |
| `src/validate.ts` | `validateVariantInput` (the two input styles, alleles `^([ACGTacgt]{1,50}\|-)$`, `ref !== alt`), `validateGenotypingParams` + `GENOTYPING_PARAM_LIMITS` (the accepted subset only — `num_return`, `max_ns` and the junction params are **rejected** by the server), and the hint that the effective product minimum is `2 × max_size + 1` |
| `src/presets.ts` | `GENOTYPING_PRESETS` (`kasp`, `as_pcr` level 0), `GENOTYPING_LADDER` (display only), `changedGenotypingAssay`, `changedGenotypingParams`. `ParamsPanel.tsx:188` hard-codes `['pcr','qpcr']` — parameterize it |
| `src/state.ts` | `ALL_MODES` (`:19`) and `availableModes` (`:42`) take `DesignerMode`; `normalizeDesignerState` (`:145`) gains `cleanGenotyping(raw)`; `designerIdentity` (`:234`) unchanged unless you want the variant key in it |
| `src/results.ts` | `matchGenotypingResults(sets, job, submitted)` and `submittedGenotypingSets(req)` beside `matchCheckResults` (`:44`) and `submittedPairs` (`:84`). Reuse `pairKey` (`:15`) |
| `src/cost.ts` | `GENOTYPE_CPU_S_PER_GENOME = 0.2`; `CheckCpuInput` gains `genotyping?: boolean`; `CheckCpuEstimate` gains `genotyping_cpu_s` (§4) |
| `src/exporters.ts` | `orderSheetToTSV` + `ORDER_TSV_HEADER`, `genotypeCallsToTSV` + `GENOTYPE_CALLS_TSV_HEADER`, `orderRowsToFasta` (names from `order[].name`, unique by construction). Every cell through `tsvCell` (`:24`) |
| new `src/genotyping.ts` | display metadata the way `src/pangenome.ts` does it: `ALLELE_META`, `PREDICTION_META`, `PRIMER_STATUS_META`, `alleleMatrixRows`, `emptyGenotypeSummary`, `isDisagreement` |
| `src/components/**` | new `GenotypingPanel`, `VariantPicker`, `ManualVariantInputs`, `AssayOptions`, `SetsTable`, `SetDetail`, `OrientationExplain`, `AlleleMatrix`, `AlleleLegend`. Changed: `PrimerDesigner.tsx` (mode branch near `:414-462`, tabs at `:265-274`, exports at `:322-334`), `ModeTabs.tsx:16`, `util.ts:46` `MODE_LABELS`, `designChecks.ts` (a genotyping case in `analyzeDesignInputs`/`restoreBlockers`), `TemplateMap.tsx` (§3.4), `Warnings.tsx` `WARNING_TEXT` (`:5`, add the new codes), `components/index.ts` barrel, `src/index.ts` public exports |
| `examples/playground/**` | `pages.ts` (`MockVariant`, new pages), `mockClient.ts` (the three methods + a genotyping job simulation), `App.tsx:120` mode list |
| `src/styles/primers.css` | matrix/lane/dye styles; keep glyph + text, never colour alone |

### 2.3 Types — corrected against the real captures

Start from spec §8.2 (`docs/genotyping_design_spec.md:3279-3421`); it is close. The table below is what the **captures**
show, and it wins wherever the two differ. Key sets are exact and complete.

| Object | Exact keys (capture-verified) |
| --- | --- |
| design response | `variant, template, assay, neighbours, orientations, sets, check, settings, engine, warnings` |
| `variant` (design) | `key, requested_id, ids, synonyms, label, kind, region, vcf, minimal, alleles, multiallelic, shift, zone, discriminating, records, ems, consequence, ref_verified, designable, issues, submission_sequence` — listing/lookup entries are the same **without** `requested_id` and `submission_sequence` |
| `template` | `system_name, region, start, end, strand, length, alt_length, seq, alt_seq, masked, mask_source, mask, masked_fraction, features`; `features = {variant, zone, discriminating {forward, reverse}, alt_offset, exempt}` |
| `assay` | the 8 request fields + `ems_target, kasp_mix`; `kasp_mix = {stock_uM, as_ref_uL, as_alt_uL, common_uL, water_uL, total_uL, source}` |
| `orientations.{forward,reverse}` | `status, reason, discriminating_position, relaxation_level, sets_found, blockers, attempts`; attempt = `level, changes, explain, pairs_returned, rejected, not_scored, sets` |
| `sets[]` | `id, key, rank, orientation, relaxation_level, quality, score, primers, products, thermo, tm_balance, neighbour_sites, primer3_penalty, warnings, issues, check, order` |
| oligo (`primers.as_ref` / `as_alt` / `common`, identical key sets) | `role, allele, three_prime_base, haplotype, target_seq, matched_seq, dye, tail_seq, order_seq, len, order_len, tm, tm_method, matched_tm, gc, hairpin_th, self_any_th, self_end_th, end_stability, primer3_problems, template, genomic, inserted_bases, deliberate_mismatch, discrimination, tailed, neighbours` |
| `order[]` row | `name, set_id, set_key, role, allele, dye, order_seq, target_seq, tail_seq, length, tm, gc, orientation, product_size_ref, product_size_alt, variant_key, notes` |
| `check` | `request, set_ids, unique_primers, omitted_set_ids`; `request = {system_name, mode, checks, pairs, genotyping}`, `genotyping = {variant, sets}` |
| `results.genotyping` | `algorithm_version, variant, summary, genomes, sets` |
| `genotyping.genomes[]` | `system_name, display_name, is_reference, allele, observed, source, copies, orthologous_copies, paralog_copies, reason`; copy = `region, start, end, strand, variant_position, identity, gap_compressed_identity, aligned_length, observed, flank_edits, call, anchors, ortholog, source` |
| `genotyping.sets[]` | `id, ref_pair, alt_pair, orientation, deliberate_mismatch_positions, specificity, control, reference, summary, genomes`; prediction row = `system_name, ref_primer, alt_primer, common_primer, predicted, strength, agrees, reasons, off_locus_products`; primer call = `status, likelihood, mm_pos, residual_mm_pos` |

**Every field the captures actually show as `null`** — type these as nullable and render them:

| Nullable | When |
| --- | --- |
| `variant.requested_id` | manual variant (`capture-genotyping-design-manual-deletion`) |
| `variant.multiallelic` | biallelic site (all captures) |
| `variant.discriminating.{forward,reverse}.alt_maps_to` | an inserted ALT base has no reference coordinate (insertion capture) |
| `variant.zone`, `variant.discriminating` (whole objects) | issue `REPEAT_TOO_LONG` only — **not in any capture**, but real listings can contain it |
| `template.mask_source` | `avoid_repeats` off (the default) |
| `orientations.*.reason` | `status: "ok"` |
| `orientations.*.relaxation_level` | `status: "blocked"` / `"skipped"` (insertion capture, forward) |
| oligo `matched_tm`, `primer3_problems`, `deliberate_mismatch`, `tailed`, `tail_seq`, `dye`, `allele`, `discrimination` | no mismatch / no tail / common primer (`allele`, `discrimination`, `dye`, `tail_seq`, `tailed` are always null on `common`) |
| `order[].allele`, `order[].dye`, `order[].tail_seq` | the common-primer row |
| `assay.kasp_mix` | `type: "as_pcr"` |
| `orientations`, `check` | `template_only: true` → `orientations: null`, `sets: []`, `check: null`; `check` is also `null` when no set was found |
| `engine.primer3`, `engine.thermo`, `engine.variation_source` | `template_only`, or no variation source |
| `copies[].ortholog` | always `null` in region mode (only gene mode annotates orthologs) |
| `genomes[].reason`, `genomes[].observed`, `genomes[].source` | `reason` is null for a successful call; `observed`/`source` are null for `missing`/`unavailable`/`ambiguous` |
| primer call `likelihood`, `mm_pos`, `residual_mm_pos` | `status: "no_product"`, or an approximate alignment (`status: "unknown"`) |
| `strength`, `agrees` | see §5 — both are routinely `null` |
| job fields | `queue_position` (running/done), `results.transcriptome` (non-transcript), `error` |

Other capture-derived corrections to §8.2: `deliberate_mismatch_positions` is `[]` for KASP; `tm_method` is `"primer3"`
in every capture (`"ntthal_duplex"` only with a deliberate mismatch); `settings.pinned` is `[]` unless the client sends
params; `warnings` and `issues` are always arrays (possibly empty); `neighbours[].distance_from_3p` is filled in even for
orientation `blockers` (spec said `null`).

### 2.4 State stays `v: 1` and additive

```ts
export interface GenotypingState {
  variantId?: string; alt?: string; variantKey?: string;
  manual?: { region: string; position: number; ref: string; alt: string };
  window?: { region: string; start: number; end: number };          // <= 50,000 bp
  filters?: { types?: VariantKind[]; includeEms?: boolean; query?: string };
  assay?: Partial<GenotypingAssay>; params?: Partial<GenotypingParams>;
  label?: string; avoidRepeats?: boolean; repeatMaskMode?: RepeatMaskMode;
  designed?: boolean; selectedSetKey?: string; checkedSetKeys?: string[];   // <= 5, /^[0-9a-f]{12}$/
  check?: { checks: CheckName[]; genomes?: string[]; params?: Partial<CheckParams>;
            jobId?: string; submitted?: SubmittedGenotypingSet[] };
  view?: { tab: 'sets' | 'alleles' | 'specificity' | 'pangenome' | 'order' };
}
// PrimerDesignerState: mode: DesignerMode;  genotyping?: GenotypingState;
```

- **`v: 1` can stay.** Every addition is a new optional key, and `normalizeDesignerState` (`state.ts:145-211`) already
  drops or clamps anything it does not recognise; an older build reading `mode: 'genotyping'` falls back to its first
  available mode (`state.ts:149`) and ignores the slice. Do **not** bump `v` — gramene-search stores this state per gene
  (`src/bundles/uiViewState.js:147`) and in saved views (`src/bundles/viewSnapshot.js:104`), and a bump discards them.
- **Do not widen `DesignMode`.** Introduce `DesignerMode = DesignMode | 'genotyping'` and apply it only to the designer
  surface. Widening `DesignMode` would let `mode: 'genotyping'` reach `POST /primers/check`, which swagger rejects with
  `ENUM_MISMATCH`, and would leak into `PrimerTemplate.mode` and `estimateCheckCpu`.
- Store **no response data**: `checkedSetKeys` hold `sets[].key` (12 hex, stable across re-designs), and `submitted`
  holds uppercase target sequences so results can be matched after a reload.

---

## 3. The UI surface

### 3.1 Variant picker

| Element | Behaviour |
| --- | --- |
| Availability | `GET /primers/genomes` → `variation.available` (and `genomes[].has_variation` for the chosen genome). False ⇒ hide the picker, offer manual entry; badge "Ensembl 115 variants" / "manual entry only" |
| Window | default gene span ± 2 kb, clamped to **50,000 bp** (`400 VARIANT_WINDOW_TOO_LONG {length, max}`); `start ≤ end ≤` region length (`400 REGION_OUT_OF_BOUNDS`) |
| Listing | `listVariants` debounced 300 ms, aborted on change. `total` / `returned` / `truncated` + warning `VARIANTS_TRUNCATED {returned, total, limit}` → "showing 2,000 of N; narrow the window". Sorted by `vcf.position`, then `key`; `limit` default 2000, max 5000 |
| Row | `vcf.position`, `label` (`1:11283-11283 A/-`), `kind`, `ids` (first + "+n"), `records[].source` badges (EMS hollow; `ems: true` means *every* record is EMS), `consequence`, and status: `designable: true` ✓ or the first `issues[].code`. Non-designable rows are listed on purpose — show why |
| Filters | `types` (CSV of `snv,mnv,insertion,deletion,complex`), `include_ems=false` to drop all-EMS entries, plus a client-side text search over ids and labels |
| Selection | store `variantKey` **and** `vcf`; design with `variant: {region, position, ref, alt}` from `vcf`. Keep the id as a label only |
| Deep link / id lookup | `getVariant(id)` → `{requested_id, system_name, source, variants[], warnings}`. It returns **one entry per alternative allele, including non-designable ones**; `synonyms` lists Ensembl's synonyms of the requested id; warning `DUPLICATE_VARIANT_IDS {key, ids}` when several ids describe one event. URL-encode the id (real ids contain `,` and `*`: `tmp_1_13549_TTA_T%2C%2A`) |
| Multi-allelic site | `multiallelic {alleles, other_alts}` non-null ⇒ make the user pick `alt`; designing by id without it is `400 ALT_REQUIRED {id, alts}`, a wrong one `400 ALT_NOT_AT_SITE {id, alt, alleles}` |
| Manual entry | VCF style (neither allele `-`) or Ensembl style (exactly one `-`; for an insertion `position` is the base *after* the insertion point). Offer **Check reference** = `template_only: true`, which verifies REF and returns the normalized `variant`, `template`, `assay`, `neighbours`; on `400 REF_MISMATCH` mark the ref field from `details.genome` |
| `shift > 0` | show "the indel can slide N bp" — the forward and reverse primers then end at different bases |

### 3.2 Assay options

| Control | Values | Default (kasp / as_pcr) | Note |
| --- | --- | --- | --- |
| `type` | `kasp`, `as_pcr` | `kasp` | any combination of the rest is accepted on either type |
| `orientation` | `both`, `forward`, `reverse` | `both` / `both` | the unrequested one reports `status: "skipped", reason: "not_requested"` |
| `tails` | `none`, `ref_fam_alt_hex`, `ref_hex_alt_fam` | `ref_fam_alt_hex` / `none` | FAM `GAAGGTGACCAAGTTCATGCT`, HEX `GAAGGTCGGAGTCAACGGATT`; hide behind "advanced" for as_pcr |
| `deliberate_mismatch` | `none`, `auto` | `none` / `auto` | |
| `mismatch_position` | `2`, `3` | `2` / `2` | label −3 "extrapolated" (Little 1995 Table 9.8.1 is defined for −2) |
| `num_sets` | 1–10 | 6 | only 4 sets fit a full-panel check |
| `max_relaxation` | 0–2 | 2 | a relaxed set carries `relaxation_level` and warning `RELAXED_CONSTRAINTS` |
| `neighbour_policy` | `avoid_3p`, `ignore` | `avoid_3p` | `ignore` turns blocking into issue `NEIGHBOUR_AT_3P` (high) |
| `avoid_repeats`, `repeat_mask_mode` | bool, `n_mask`/`three_prime` | `false`, `n_mask` | reuse `RepeatOptions`; the AS window is always exempt (`template.features.exempt`) |
| `params` | the subset in `GENOTYPING_PARAM_LIMITS` | preset level 0 | `product_size_ranges` is 1–4 ranges of 20–1000; a param the client sends is **pinned** (never relaxed) and echoed in `settings.pinned` |

Show `settings.ladder` and `settings.floors` (`as_min_tm 52`, `as_min_gc 15`) read-only, and `assay.ems_target` when
every record of the target is EMS (warning `EMS_TARGET`: neighbours warn instead of blocking).

### 3.3 Sets table (`SetsTable`, replacing `PairsTable` in this mode)

One row per set, keyed by **`set.key`** (12 hex, `sha256` of orientation + the three target sequences — stable across
re-designs; `set.id` is positional and is *not* an identity).

| Column | Source | Rendering |
| --- | --- | --- |
| select | — | disabled when it would exceed 5 sets / 10 pairs / 20 distinct primers, or when the mirrored cost is `over_limit` |
| set | `id`, `orientation`, `relaxation_level` | "S1 · reverse · L0"; an arrow for orientation |
| REF primer | `primers.as_ref` | `target_seq` with the 3′ base highlighted; `allele` chip ("C"); `dye` chip (`FAM`/`HEX`/none); copy buttons for **`target_seq`** and **`order_seq`** separately |
| ALT primer | `primers.as_alt` | same; `allele` may be multi-base for an indel (`"CGT"` in the insertion capture); `inserted_bases > 0` marks bases with no reference coordinate |
| Common primer | `primers.common` | `allele`, `dye`, `tail_seq`, `tailed`, `discrimination` are all `null` here |
| tail | `tail_seq` / `order_seq` | **`target_seq` is what anneals and the only sequence a check ever receives; `order_seq = tail_seq + target_seq` is what the vendor synthesizes.** Never send `order_seq` to `/primers/check` — a tailed oligo fails the check's `^[ACGTacgt]{15,36}$` pattern |
| product | `products.ref.size` / `products.alt.size` | "65 / 65 bp"; an indel differs by `template.features.alt_offset` |
| Tm | `primers.*.tm`, `tm_balance` | REF/ALT/common with `as_tm_diff` and `common_minus_as`; `tm_method` `"ntthal_duplex"` ⇒ label "template duplex", with `matched_tm` as the perfect-match Tm |
| neighbours | `neighbour_sites` | count badge; detail from each oligo's `neighbours[]` (`distance_from_3p`, 1 = the 3′ base; EMS entries are informational) |
| quality | `quality` | `good` / `usable` / `poor` chip. `poor` = a high-severity issue; `usable` = a warn issue or `relaxation_level ≥ 1`. Quality does **not** change the order |
| issues | `issues[]` (`{code, severity, message, details}`) | one chip per code with text, never colour alone. `warnings[]` repeats the warn/high ones with fuller messages (and `details.severity` for `TAILED_STRUCTURE`/`MISMATCH_STRUCTURE`) |
| score | `score`, `primer3_penalty` | lower is better; the rows are already in rank order (`rank`, 0-based) |
| after a check | §4 | agreement chip ("11/11 agree") + the two `specificity` verdicts |

`OrientationExplain`: one card per orientation with `status`, `reason`, `blockers[]` ("blocked by rs5413863234 (G/A),
4 nt from the 3′ end") and the `attempts[]` table (`level`, `changes`, `explain`, `pairs_returned`, `rejected.*`,
`not_scored`, `sets`) through the existing explain renderer (`src/explain.ts`). `pairs_returned` always equals
`rejected` + `not_scored` + `sets`.

### 3.4 Template map: the variant track

`TemplateMap` currently assumes exactly one left and one right primer per row (`packPairLanes` at `TemplateMap.tsx:34-49`,
the two arrows at `:375-376`). For genotyping it needs a variant track and three-primer lanes:

| Element | Source (template coordinates, 1-based, plus strand) |
| --- | --- |
| variant marker | `template.features.variant {start, end}` |
| zone band | `template.features.zone {start, end}` — only an AS primer may touch it |
| F / R discriminating ticks | `template.features.discriminating {forward, reverse}` — draw both when they differ (they do for every indel: deletion 403 vs 401, insertion 402 vs 401) |
| ALT haplotype | `template.alt_seq`, `alt_length`, `features.alt_offset` (`+2` insertion, `−1` deletion, `0` SNV) |
| mask exemption | `features.exempt` as `[start, length]` pairs (e.g. `[[365,37],[401,37]]`) |
| neighbour ticks | the `neighbours[]` of the three oligos (filled = non-EMS, hollow = EMS), or a `listVariants` call over `template.start…end` |
| per-set lane | REF and ALT arrows stacked on the same footprint with 3′ tips coloured by allele, plus the common arrow facing the other way; hatch a `blocked`/`skipped` orientation |

Genomic footprints come from each oligo's `genomic {region, start, end, strand, blocks[]}`; an ALT primer spanning a
deletion has **two blocks** (deletion capture, `as_alt`: `[{11257,11282},{11284,11286}]`).

### 3.5 Allele matrix (results tab "Alleles")

Rows: `results.genotyping.genomes[]` — the reference first (`is_reference: true`), then the finished pan-genome genomes in
`request.genomes` order (sorted by `system_name`). Columns: the allele, then one per submitted set.

| Cell | Values | Extras to show |
| --- | --- | --- |
| allele (`genomes[].allele`) | `ref` ● · `alt` ▲ · `other` ◆ · `ambiguous` ? · `missing` – · `unavailable` × | `observed` core (compare against `variant.core.ref` / `.alt`), `orthologous_copies` as a superscript, `paralog_copies`, `source` (`amplicon` / `megablast`), and `reason` for `missing`/`unavailable` |
| prediction (`sets[].genomes[]`) | `ref` ● · `alt` ▲ · `both` ◐ · `none` ○ · `no_call` ⊘ · `unknown` · | `agrees` ✓/✗/– , `strength: "weak"` marker, ⚠ when `off_locus_products > 0`, tooltip = `reasons[]`, and the three primer calls (`status`, `likelihood`, `mm_pos`, `residual_mm_pos`) |
| detail panel | `copies[]` table | `region:start–end`, `strand`, `identity` **and** `gap_compressed_identity` (the ortholog test uses the gap-compressed one, ≥ 95 %), `aligned_length`, `observed`, `flank_edits`, `call`, `anchors`, `source` |

- **Glyph + colour + text**, never colour alone; reuse the grid keyboard handler of `PangenomeMatrix.tsx:307` and the row
  model idea of `buildMatrixRows` (`:49`).
- Filters: "disagreements only" (`agrees === false`) and "issues only" (`allele` not `ref`/`alt`, or `predicted` not
  `ref`/`alt`).
- Summaries: `genotyping.summary` (`ref + alt + other + ambiguous + missing + unavailable = genomes_total`) and per set
  `predicted_ref + predicted_alt + both + none + no_call + unknown = genomes_total` and
  `agree + disagree + not_comparable = genomes_total`; `weak` counts `strength: "weak"`. **The reference is never counted
  in a summary.**
- `control` (`{status, allele, reasons}`) belongs above the matrix: `pass` / `warn` / `fail`. A `fail` also raises results
  warning `REFERENCE_CONTROL_FAILED {allele, sets}` — say plainly that that set's other predictions are unreliable.
- **Off-locus threshold:** an off-locus product changes a prediction only when each primer has ≤ 2 mismatches
  (`check.genotype_offlocus_max_mismatches`). Weaker ones go to the job-level warning `WEAK_OFF_TARGETS
  {count, max_mismatches, examples}` (≤ 5 examples: `system_name, set_id, pair_id, allele, region, start, end, size,
  orientation, left_mm, right_mm`) and **still appear** in `specificity.off_targets` and the pan-genome results. Explain
  that the warning is informational, not a failure ([docs](primer_design_api.md#off-locus-products-and-weak_off_targets)).
- **Partial jobs:** `results.genotyping` first appears after the reference stage (reference entry + each set's
  `reference`, `control`, `specificity`), then grows genome by genome; the example job showed 1, 2, 9, 10 and finally 12
  entries. Render growth, not an empty grid.

### 3.6 Order sheet

| Piece | Source | Note |
| --- | --- | --- |
| order rows | `sets[].order[]`, three rows per set in REF, ALT, common order | `name = {label}_{set id}_{REF\|ALT\|COM}` plus `_FAM`/`_HEX` for a tailed AS primer (`rs871475760_S1_REF_FAM`); `notes` carries "deliberate mismatch −2 A→G (Little 1995 Table 9.8.1)" or "tail hairpin 51.5 °C". Ship `order_seq` (what is synthesized) *and* `target_seq` |
| KASP mix | `assay.kasp_mix` | 12 µL REF + 12 µL ALT + 30 µL common at 100 µM + 46 µL water = 100 µL, with its `source` string; `null` for `as_pcr` |
| submission string | `variant.submission_sequence` | 50 bp of reference each side of `[REF/ALT]`, minimal alleles with `-` written as nothing (`[A/]`, `[/GT]`); biallelic non-EMS SNV neighbours appear as IUPAC codes. Warning `SUBMISSION_NEIGHBOURS_OMITTED {ids}` means indel or multi-allelic neighbours were left as reference bases |
| exports | order sheet TSV/CSV (with `#` note lines for the mix and the submission string), genotype calls TSV (`genomes × sets`), the submission line, primers FASTA named from `order[].name` | all cells through `tsvCell`; `order[].name` is unique, which also fixes the duplicate `G_P1_F` FASTA names of the pairs exporter |

---

## 4. Check integration

- **The design hands you a ready body.** `check.request` is a complete `POST /primers/check` body:
  `{system_name, mode: "region", checks: ["specificity","pangenome"], pairs, genotyping}` with no `genomes` (so the whole
  panel is searched) and no `params`. Posting it unchanged is the simplest correct path; add `genomes` to narrow it.
- **Pair ids** are `S{n}_REF` and `S{n}_ALT`; `genotyping.sets[] = {id, ref_pair, alt_pair}`; `expected` is the **REF**
  product and is identical on both pairs of a set; `genotyping.variant` is `{region, ...variant.vcf}` (VCF style, no `-`).
- **Match results to sets by uppercase `target_seq`, never by id or rank.** Build the triple
  `pairKey(ref.left, ref.right) + '#' + pairKey(alt.left, alt.right)` from `job.request.pairs` grouped through
  `job.request.genotyping.sets`, and look up each design set's identical triple (`matchGenotypingResults`). Never match on
  `order_seq` — the tail is not in the job. Design sets with no match are `notChecked`; submitted triples with no design
  set are orphans ("results for primers not in the current design"). This mirrors `matchCheckResults` (`results.ts:44-81`)
  and survives a re-design that renumbers sets.
- **Caps.**

  | Cap | Value | Where |
  | --- | --- | --- |
  | sets in the server-built `check.request` | 5 | `genotyping.check_max_sets`; the rest are listed in `check.omitted_set_ids` |
  | pairs per check job | 10 | `check.max_pairs`; `genotyping.sets` is 1–5 |
  | distinct primers in the server-built request | 13 | `genotyping.check_max_unique_primers` — the most that fit under 6,000 CPU-s over the **full panel** |
  | distinct primers accepted by the check endpoint | 20 | `400 TOO_MANY_PRIMERS` |

  Gate the Submit button on the mirrored cost plus 10 pairs / 20 primers, not on a fixed primer count: over three genomes
  all 20 primers fit; over the full panel the limit bites at 14.
- **Cost mirror** (`src/cost.ts`): add `genotyping_cpu_s = input.genotyping ? (1 + pan.length) * 0.2 : 0`
  (`check.genotype_cpu_s_per_genome`, charged for the reference plus each pan-genome genome) to the existing terms. Real
  server values to pin ([docs](primer_design_api.md#post-primerscheck-the-genotyping-block)):

  | Workload (region mode, sorghum) | Without `genotyping` | With |
  | --- | --- | --- |
  | 1 set (3 primers), 119 assemblies | 1,333 | 1,357 |
  | 2 sets (6 primers), 119 | 2,666 | **2,690** |
  | 4 sets (12 primers), 119 | 5,332 | 5,356 |
  | 13 distinct primers, 119 | 5,776 | 5,800 (the most accepted) |
  | 14 distinct primers, 119 | 6,221 | 6,245 → `422 JOB_TOO_LARGE` |
  | 2 sets (6 primers), 3 genomes | 95 | 95 |

  Use the real sizes in parity tests: reference 708,735,318 bases; the 119 others 83,774,241,795.
- **The estimate is deliberately conservative.** The documented example job estimated **273 CPU-s** and measured
  **176.8** — about **1.5×**. The allele caller itself measured 0.03–0.10 CPU-s per genome (0.005–0.057 in the full-panel
  run) against the 0.2 charged. Label the number "estimated", never "expected runtime".
- **Job ids and polling are unchanged:** deterministic ids, `202` newly queued vs `200` existing, back off 1 s ×1.5 up to
  10 s, `404 UNKNOWN_JOB` ⇒ "results expired, re-run". Genotyping jobs hash with `"2+g1"`, so a caller-version change
  invalidates only genotyping jobs.
- **Submit-time rejections are user errors, not bugs.** `400 GENOTYPING_SET_INVALID` with
  `details {set_id, reason, pair_ids, primer?, sequence?, mm_pos?}` — reasons `unknown_pair`, `same_pair`, `pair_reused`,
  `expected_required`, `expected_differs`, `no_shared_common`, `not_at_variant`, `alleles_swapped`, `too_many_edits`,
  `common_in_zone`, `expected_mismatch`. They only happen if the UI edits the sets, so surface the message verbatim.
  `genotyping` outside gene/region mode is `400 INVALID_REQUEST {field: "genotyping", reason: "mode"}`.
- **Defensive state:** a `done` job whose `request.genotyping` is set but whose `results.genotyping` is missing means an
  older worker served it — show "Allele check unavailable (server update in progress)" with a Re-run button, not an empty
  matrix.

---

## 5. Behaviour the UI must not get wrong

| # | Behaviour | Why it is like that |
| --- | --- | --- |
| 1 | **`strength` stays `null` when a prediction comes only from off-locus products** — and also for every `unknown`, `no_call` and `none` row. Do not render "normal" as a default | `strength` describes the on-locus signal; with no on-locus product there is nothing to grade. Rows 1–4 and 8 of the prediction table have no strength ([docs](primer_design_api.md#amplification-prediction)) |
| 2 | **A genome can be `predicted: "both"` although its allele is a clean single call** | `sorghum_pi180348` carries a paralog at 1:72,907,129–72,907,210 that S2's REF pair amplifies with 1 + 2 mismatches. Its allele is `alt` (two orthologous copies), the prediction is `both` with `reasons: ["ref_signal_off_locus"]`, `off_locus_products: 1`, `agrees: false`, `paralog_copies: 1`. Present it as "this set also lights the other dye from a paralog", not as a data error |
| 3 | **Deletion ALT genomes can be `unknown (shift_tract_uncertain)` by design** | when an AS primer's 3′ base lies inside an indel's shift tract (`discrimination.in_shift_tract: true`), the wrong-allele primer may still prime through a bulge, so the prediction is `unknown`, `strength` `null`, `agrees` `null`, and the genome is excluded from concordance. In the M11 deletion panel all 5 ALT genomes are `unknown`. Never count those as disagreements |
| 4 | **`missing` carries a reason** — `no_orthologous_copy`, `variant_not_covered`, `fallback_failed`, `fallback_budget` — and the prediction is then `none` with reason `no_orthologous_copy` | they mean different things: "no comparable locus in that assembly" vs "the megablast fallback failed or ran out of budget" (results warnings `GENOTYPE_FALLBACK_FAILED` / `GENOTYPE_FALLBACK_BUDGET`). In the full panel, `is36143` and `pi525695` are `no_orthologous_copy` everywhere. Conversely, a `disagree` row may carry an **empty** `reasons` array, so never explain a disagreement from `reasons` alone |
| 5 | **`specificity` and `pangenome` cannot answer allele questions** | both S1 pairs "amplify" in all 11 assemblies (`summary.amplifies: 11`), because a single 3′-terminal mismatch still gives a `likely_weak` product. Only `results.genotyping` says which allele a genome carries and which dye would light |
| 6 | **The Ensembl feature flag and outages degrade, they do not fail the feature** | `PRIMERS_VARIATION_ENABLED=0` ⇒ both variants endpoints and *design by id* answer `503 FEATURE_DISABLED` (retry after 300 s), `genomes` reports `variation.available: false`, and manual designs keep working with `neighbours.data: "none"` + warning `NO_VARIATION_DATA` (whose message says lookups are **disabled on this server**). An Ensembl outage gives `503 VARIATION_SOURCE_UNAVAILABLE {retry_after_s, reason}` (`timeout`, `transport`, `http_5xx`, `rate_limited`, `invalid_response`, `breaker_open`, `queue_full`) on the pickers and on design-by-id, while a *manual* design returns `200` with `neighbours.data: "unavailable"` and warning `NEIGHBOURS_UNAVAILABLE {reason}` — meaning "no neighbour screening", not "bad primers". Always offer manual entry as the fallback and honour `retry_after_s` with a countdown ([docs](primer_design_api.md#ensembl-and-what-happens-without-it)) |
| 7 | **`DESIGN_BUDGET_EXHAUSTED` is still a `200`** | the per-request budget (54 Primer3 runs, 272 `ntthal` calls) can leave candidates unscored: the response carries the sets found so far plus warning `DESIGN_BUDGET_EXHAUSTED {orientation, primer3_runs, thermo_calls, max_primer3_runs, max_thermo_calls, not_scored, sets_returned}`, and an orientation with nothing left reports `status: "no_sets", reason: "budget_exhausted"`. Never treat it as an error banner; show the sets and a note |
| 8 | **`sets: []` is a `200` too** | warning `NO_SETS`; the reason lives in `orientations.*.attempts[].explain` and in `blockers`. Show the orientation cards, not "request failed" |
| 9 | **Send `target_seq`, keep `order_seq` for ordering** | the tail is not part of the assay's specificity and fails the check's primer pattern |
| 10 | **Identify a set by `key`, not `id`** | `S1`/`S2` are positional; `key` is content-derived and stable across re-designs |
| 11 | **`results.genotyping` is absent, not `null`, for ordinary jobs** | so old code paths see exactly the pre-genotyping results |

---

## 6. Fixtures for offline work

### 6.1 Real captures (`test/primers/fixtures/docs/`)

Wrapper shape: `{source, request {method, path[, body]}, status, response}`. **Read-only**: `docs_examples.test.js`
compares the documented examples against these files and fails if a capture is not cited, so copy them into
`gramene-primers/test/fixtures/genotyping/` rather than editing them here.

| File | Size | What it shows |
| --- | --- | --- |
| `capture-variants-list-1_11180-11290.json` | 3.0 KB | 4 entries in `1:11180–11290`: 2 EMS SNVs, 1 EVA SNV, 1 shiftable deletion (`shift: 2`, zone 11282–11286). `include_ems=false` must leave 2 |
| `capture-variants-list-no-variation-data.json` | 0.5 KB | `422 NO_VARIATION_DATA` on `sorghum_rio` |
| `capture-variants-lookup-tmp_1_11502_C_CGT.json` | 1.3 KB | id lookup of an insertion, 2 merged ids, warning `DUPLICATE_VARIANT_IDS`, `alt_maps_to: null` |
| `capture-variants-lookup-unknown-variant.json` | 0.4 KB | `404 UNKNOWN_VARIANT` |
| `capture-genomes-sorghum_bicolor.json` | 36 KB | the additive `variation` block and 120 genomes, exactly 1 with `has_variation: true` |
| `capture-genotyping-design-rs871475760-kasp.json` | 18 KB | **the main design fixture**: KASP, 2 sets — S1 reverse 65 bp `usable` 9.18, S2 forward 92 bp `poor` 19.63; issues `ALT_PRIMER_SUBOPTIMAL`, `TAILED_STRUCTURE` (warn 51.5 °C and high 66.22 °C), `NEIGHBOUR_IN_PRIMER`; full `check.request` (4 pairs, 6 unique primers), order rows, `kasp_mix`, submission string |
| `capture-genotyping-design-tmp_1_11502_C_CGT.json` | 14 KB | insertion: forward orientation `blocked` (`neighbour_at_3p`, `distance_from_3p: 4`), 1 set, `AS_TM_IMBALANCE`, `WEAK_TERMINAL_CLASS`, response warnings `DUPLICATE_VARIANT_IDS`, `ORIENTATION_BLOCKED`, `DENSE_NEIGHBOURS`; ALT primer ends on an inserted base (`allele: "CGT"`, `inserted_bases: 1`) |
| `capture-genotyping-design-manual-deletion.json` | 22 KB | manual Ensembl-style deletion: `requested_id: null`, ids from Ensembl, both orientations at `relaxation_level: 1`, warnings `SHIFTABLE_INDEL`, `DENSE_NEIGHBOURS`, `RELAXED_CONSTRAINTS`; issues `SHIFT_TRACT_DISCRIMINATION`, `NEIGHBOUR_IN_PRIMER`, `WEAK_TERMINAL_CLASS`; ALT primer with two `genomic.blocks`; `alt_offset: -1` |
| `capture-genotyping-design-ref-mismatch.json` | 0.6 KB | `400 REF_MISMATCH {region, position, given, genome}` |
| `capture-check-genotyping-submit.json` | 1.7 KB | real `202`: `{job_id, status, kind: "pangenome", queue_position, progress {done:0,total:12}, estimate {cpu_s: 273}, created_at, warnings}` |
| `capture-check-genotyping-result.json` | 66 KB | **the main results fixture**: `done` job over the reference + 11 assemblies, 2 sets. 6 `ref` / 5 `alt`; 4 genomes with 2 copies; pi180348 `paralog_copies: 1`; S1 agrees 11/11, S2 10/11 (pi180348 `both`); `WEAK_OFF_TARGETS` count 28 with 5 examples; `timings_ms`; `transcriptome: null` |
| `capture-check-genotyping-no-shared-common.json` | 1.3 KB | `400 GENOTYPING_SET_INVALID` reason `no_shared_common` |
| `capture-check-genotyping-not-at-variant.json` | 1.4 KB | `400 GENOTYPING_SET_INVALID` reason `not_at_variant` |

### 6.2 What the captures do **not** cover — synthesize these in the mock

The result capture only exercises `allele` ∈ {`ref`,`alt`}, `source: "amplicon"`, `copy.call` ∈ {`ref`,`alt`},
`predicted` ∈ {`ref`,`alt`,`both`}, `strength: "normal"`, primer statuses {`match`,`terminal_mismatch`},
`control: "pass"`, and `reason: null` everywhere. The mock must add: `other` / `ambiguous` / `missing` (each reason) /
`unavailable`; `source: "megablast"`; `predicted: "none"` / `"no_call"` / `"unknown"`; `strength: null` and
`"weak"`; statuses `weak`, `uncertain`, `blocked`, `no_product`, `unknown`; `agrees: null`; `control: "warn"` / `"fail"`
with `REFERENCE_CONTROL_FAILED`; a partial job (reference only, then growing); `REPEAT_TOO_LONG` entries with
`zone: null`; `template_only` responses; `as_pcr` sets with `deliberate_mismatch` and `tm_method: "ntthal_duplex"`;
`DESIGN_BUDGET_EXHAUSTED`; `NO_SETS`; and the 503 paths of §5 row 6.

### 6.3 Driving the playground mock

- `examples/playground/mockClient.ts` (`MockPrimersClient`, `createMockClient`) gains `listVariants`, `getVariant` and
  `designGenotyping` replaying the captures (import the JSON and return `capture.response`), plus a genotyping job
  simulation: reuse `snapshot()`/`results()` and emit `results.genotyping` after the reference stage, then one genome per
  tick, recomputing the summaries. `MockClientOptions.failWith` gains `VARIATION_SOURCE_UNAVAILABLE`,
  `NO_VARIATION_DATA`, `REF_MISMATCH`, `FEATURE_DISABLED`.
- `examples/playground/pages.ts`: extend `MockVariant` and add pages `genotyping-rs871475760`,
  `genotyping-insertion-blocked`, `genotyping-manual-deletion`, `genotyping-no-variation`, `genotyping-ensembl-down`,
  each with an `initialState` carrying a `genotyping` slice. `App.tsx:120` must include `'genotyping'` in `modes`.
- Existing mock affordances to reuse: `?api=mock|live&page=<id>&mockError=<CODE>`, "Restore with an expired job",
  the controlled-state inspector, `delayScale` for tests.

### 6.4 Swagger contract fixtures

- **Already present** under `test/primers/fixtures/contract/requests/` (hand-written by the API session, validated by
  `contract.test.js`): `variants-list.json`, `variants-lookup.json`, `genotyping-design-{id,aspcr,manual-deletion,
  template-only,all-fields}.json`, `check-genotyping-{two-sets,gene-mode}.json`. Use them as the target shapes.
- `npm run fixtures` in gramene-primers (`test/fixtures.gen.test.ts`) must **produce the same names**, because
  `contract.test.js` routes fixtures by prefix:
  `design-*` → `POST /primers/design`, `check-*` → `POST /primers/check`, `check-genotyping-*` → the same with a
  `genotyping` block, `genotyping-design-*` → `POST /primers/genotyping/design`, `variants-*` → the GET endpoints.
- Two file shapes are accepted: a bare request body (the `design*`/`check*`/`genotyping-design*` files) or a wrapper
  `{"method", "path", "body"|"query", "expect": "valid"|"invalid"}` (the `variants-*` files — GETs have no body).
  Extend the generator's `Case` union with `kind: 'genotyping' | 'variants'` and its `manifest.json` entries
  (`definition: 'PrimerGenotypingRequest'`, and the two variant operations).
- **Caveat:** the compatibility test `compat_ids.test.js` stores a sha256 of each `check-*.json` body and fails if one
  changes ("changed since its id was recorded"); it deliberately **skips `check-genotyping-*.json`**. So regenerating the
  genotyping fixtures is safe, but do not let the generator rewrite the existing `check-*.json` bodies.

---

## 7. Every error and warning code the UI should surface

### 7.1 Errors (handler bodies `{message, code, details}`)

| Status | Code | Show |
| --- | --- | --- |
| 400 | `VARIANT_WINDOW_TOO_LONG {length, max}` | "Window too large — at most 50,000 bp" |
| 400 | `REGION_OUT_OF_BOUNDS` | "That position is outside the sequence" (`details.length` is the region length) |
| 400 | `INVALID_VARIANT {reason}` | `id_or_manual` → "give either an id or region/position/ref/alt"; `alleles` → "REF and ALT must differ and be A/C/G/T or `-`"; `allele_too_long` → "at most 50 nt" |
| 400 | `REF_MISMATCH {region, position, given, genome}` | "The genome has **{genome}** at {region}:{position}" — mark the ref field |
| 400 | `ALT_REQUIRED {id, alts}` / `ALT_NOT_AT_SITE {id, alt, alleles}` | offer the listed alleles |
| 400 | `UNSUPPORTED_ALLELE {allele}` | "`*` and N alleles cannot be designed" |
| 400 | `VARIANT_TOO_CLOSE_TO_END {region_length, needed}` | "Too close to the end of the sequence for any product" |
| 400 | `VARIANT_TOO_REPETITIVE {shift, max}` | "The indel slides more than 1,000 bp"; `shift: null` means it slid past the read window |
| 400 | `INVALID_PARAMS {param[, level]}` | name the parameter; `level` means a relaxed ladder level breaks it |
| 400 | `PRIMER3_INPUT_ERROR {primer3_error}` | "Primer3 rejected the request" |
| 400 | `GENOTYPING_SET_INVALID {set_id, reason, …}` | §4; the server message is already user-readable |
| 400 | `INVALID_REQUEST {field[, reason]}` / validator `errors[]` | use `flattenValidationErrors` (`errors.ts:52`) |
| 404 | `UNKNOWN_VARIANT {id, system_name}` | "Ensembl does not know that id" |
| 404 | `UNKNOWN_GENOME`, `UNKNOWN_REGION`, `UNKNOWN_JOB` | `UNKNOWN_JOB` ⇒ "results expired, re-run" |
| 422 | `NO_VARIATION_DATA {system_name}` | "No known variants for this genome — enter the variant by hand" |
| 422 | `VARIANT_NOT_ON_ASSEMBLY {id}` / `AMBIGUOUS_VARIANT_MAPPING {id, mappings}` | "That id does not map to this assembly" / "maps to several sequences" |
| 422 | `NO_SEQUENCE`, `AMBIGUOUS_ASSEMBLY {candidates}`, `JOB_TOO_LARGE {estimate_cpu_s, limit}` | `JOB_TOO_LARGE` already handled by `ErrorBanner.tsx:82` |
| 500 | `THERMO_FAILED {binary}`, `PRIMER3_FAILED`, `INTERNAL` | "Server error — try again" |
| 503 | `THERMO_UNAVAILABLE {binary, retry_after_s}` (60 s or 5 s), `PRIMER3_UNAVAILABLE`, `BUSY` (5 s), `MONGO_UNAVAILABLE` (30 s), `QUEUE_FULL` (60 s), `JOB_STORE_UNAVAILABLE` (5 s) | countdown from `retry_after_s`; `BUSY` already auto-retries (`useDesign.ts` `MAX_BUSY_RETRIES`) |
| 503 | `VARIATION_SOURCE_UNAVAILABLE {retry_after_s, reason}` | "Variant lookups are temporarily unavailable ({reason}) — retry in Ns, or enter the variant by hand" |
| 503 | `FEATURE_DISABLED {retry_after_s: 300}` | "Variant lookups are switched off on this server" (manual entry still works) |
| 504 | `DEADLINE_EXCEEDED` | "The design took too long — try fewer sets" |

### 7.2 Warnings — variants endpoints

`DUPLICATE_VARIANT_IDS {key, ids}` (ids merged) · `VARIATION_RECORDS_SKIPPED {count, reasons}` (malformed Ensembl records
dropped) · `REF_MISMATCHES {count}` (entries disagreeing with our FASTA) · `VARIANTS_TRUNCATED {returned, total, limit}`.

### 7.3 Warnings — design response level (in this order)

| Code | User-facing meaning |
| --- | --- |
| `DUPLICATE_VARIANT_IDS` | several Ensembl ids describe this variant; they were merged |
| `EMS_TARGET {sources}` | an EMS mutation, private to one BTx623 line; natural neighbours warn instead of blocking |
| `MULTIALLELIC_SITE {key, other_alts}` | other alleles exist here; assemblies carrying them mismatch both AS primers |
| `SHIFTABLE_INDEL {shift, forward_position, reverse_position, zone}` | the indel can slide; the wrong-allele primer may still prime through a bulge |
| `ORIENTATION_BLOCKED {orientation, ids, distances}` | a known variant sits in the last 5 nt of that orientation's AS primer |
| `NO_VARIATION_DATA {system_name}` | neighbours were not screened (no data, or lookups disabled) |
| `NEIGHBOURS_UNAVAILABLE {reason}` | Ensembl failed during this design; neighbours were not screened |
| `DENSE_NEIGHBOURS {count, window}` | more than 2 known non-EMS variants within 30 bp |
| `SUBMISSION_NEIGHBOURS_OMITTED {ids}` | indel/multi-allelic neighbours left as reference bases in the submission string |
| `VARIANT_IN_REPEAT {orientation, masked_bases}` | the repeat mask reached the allele-specific window |
| `REPEAT_MASK_FAILED`, `NO_REPEAT_MASK`, `BLAST_DEPTH_MASK`, `MOSTLY_REPEAT`, `ASSEMBLY_MISMATCH`, `AMBIGUOUS_ASSEMBLY`, `PRIMER3_WARNING` | as in `/primers/design` (already in `WARNING_TEXT`) |
| `ORIENTATION_SKIPPED {orientation, reason}` | too close to the region end, or an `N` in the AS window |
| `ORIENTATION_NO_SETS {orientation, levels_tried}` | nothing above the floors in that orientation |
| `RELAXED_CONSTRAINTS {orientation, level, changes}` | these sets needed relaxed Primer3 constraints |
| `NO_SETS {}` | no set in either orientation — read the orientation cards |
| `DESIGN_BUDGET_EXHAUSTED {…}` | the design stopped scoring candidates; the sets shown are still valid |

### 7.4 Set issues (`sets[].issues[]`, with `severity`; warn/high repeated in `sets[].warnings`)

| Code | Severity | Meaning |
| --- | --- | --- |
| `ALT_PRIMER_SUBOPTIMAL {oligo, problems}` | warn | a derived oligo Primer3 would not have picked, but it clears the hard floors |
| `AS_TM_IMBALANCE {diff}` | warn (high > 2.0 °C) | the two AS primers' Tm differ by more than 1.0 °C |
| `COMMON_TM_OUT_OF_RANGE {value}` | warn | `common_minus_as` outside −1…+3 °C |
| `TAILED_STRUCTURE {oligo\|pair, metric, value}` | warn (high ≥ 55 °C hairpin) | a tailed hairpin/dimer at or above 47 °C; 55 °C is the KASP annealing temperature |
| `MISMATCH_NOT_APPLICABLE {position}` | warn | the base at −k is the same in both AS primers, so no mismatch was applied |
| `MISMATCH_STRUCTURE {oligo, metric, value}` | warn (high ≥ 55 °C) | a deliberate-mismatch primer's own structure |
| `NEIGHBOUR_IN_PRIMER {role, ids, distances}` | warn | known variants inside a primer, outside its last 5 nt |
| `NEIGHBOUR_AT_3P {role, ids, distances}` | **high** | a known variant in the last 5 nt that was not allowed to block |
| `WEAK_DISCRIMINATION {primer, mm_pos}` | **high** | the primer is predicted to amplify the *other* allele |
| `SHIFT_TRACT_DISCRIMINATION {primer, shift}` | warn | the 3′ base lies inside the indel's shift tract |
| `WEAK_TERMINAL_CLASS {}` | info | both 3′ mismatches are in Little's weak class; not scored |

### 7.5 Check warnings

- **Submit:** `MAX_PRODUCT_SIZE_RAISED`, `EXPECTED_IGNORED`, `GENOMES_IGNORED`, plus the reference assembly's
  `ASSEMBLY_MISMATCH` / `AMBIGUOUS_ASSEMBLY`.
- **Results:** `REFERENCE_CONTROL_FAILED {allele, sets}` (that set's predictions are unreliable) ·
  `WEAK_OFF_TARGETS {count, max_mismatches, examples}` (informational; the products are still listed) ·
  `GENOTYPE_FALLBACK_FAILED` / `GENOTYPE_FALLBACK_BUDGET` (those genomes are `missing`) · `GENOTYPE_FAILED` (that genome
  is `unavailable`) · and the pre-existing `NO_FASTA_FOR_REALIGN`, `ANNOTATION_UNAVAILABLE`,
  `PANGENOME_TRANSCRIPT_MODELS_ONLY`, `TRANSCRIPT_GENE_UNMAPPED`.

Add every new code to `WARNING_TEXT` (`src/components/Warnings.tsx:5`) so a message-less warning still reads well.

---

## 8. Acceptance checklist

**Against fixtures (unit / component tests)**

1. `GET /primers/variants` on `1:11180–11290` lists 4 entries; "hide EMS" leaves 2; the deletion row shows `shift: 2`.
2. An id lookup of `tmp_1_11502_C_CGT` shows 2 merged ids and the `DUPLICATE_VARIANT_IDS` notice; an unknown id shows the
   404 message; an id containing `,`/`*` is URL-encoded and still resolves.
3. rs871475760 KASP with `num_sets: 2` renders S1 (reverse, 65 bp, `usable`, score 9.18) and S2 (forward, 92 bp, `poor`,
   score 19.63) with the `TAILED_STRUCTURE` high chip on S2's ALT primer, plus `ALT_PRIMER_SUBOPTIMAL` and
   `NEIGHBOUR_IN_PRIMER`; the map shows two lanes, the variant marker and the zone band.
4. `tmp_1_11502_C_CGT` shows the forward orientation **blocked** with the blocker text and distance 4, and its ALT primer
   marked as ending on an inserted base.
5. The manual deletion shows `requested_id: null`, both orientations at level 1, the `SHIFTABLE_INDEL` and
   `SHIFT_TRACT_DISCRIMINATION` notices, and the ALT primer's two genomic blocks.
6. Checking S1 + S2 posts exactly `check.request` (4 pairs, `genotyping.sets` S1/S2, `genotyping.variant`
   `{region:"1", position:11109, ref:"C", alt:"A"}`); the cost line reads ≈ 95 CPU-s for 3 genomes and ≈ 2,690 for the
   full panel; 14 distinct primers over the full panel disables Submit.
7. With `capture-check-genotyping-result.json`: the Alleles tab shows 12 rows (reference first), 6 `ref` / 5 `alt`,
   pi180348 with 2 copies and `paralog_copies: 1`, S1 11/11 agree, S2 10/11 with pi180348 `both`, `ref_signal_off_locus`
   and ⚠; "disagreements only" leaves exactly that one row; the `WEAK_OFF_TARGETS` notice explains itself.
8. Synthetic rows render correctly: `no_call` ⊘ with `common_primer_3p_mismatch`, `unknown` with
   `shift_tract_uncertain`, `missing` with its reason, `unavailable` ×, `strength: null` (no "normal" label),
   `agrees: null` as "–", and a `fail` control banner.
9. Results match after a re-design that renumbers sets (match by target-sequence triple), and orphan submitted sets are
   labelled.
10. Order-sheet TSV for S1 equals the three `order` rows plus the `#` mix and submission lines; FASTA names are unique;
    `tsvCell` escapes a leading `-`/`=`.
11. A saved-view round trip keeps variant, assay, `checkedSetKeys` and `jobId`; state stays `v: 1`; an old build ignores
    the slice.
12. `DESIGN_BUDGET_EXHAUSTED`, `NO_SETS`, `503 VARIATION_SOURCE_UNAVAILABLE`, `503 FEATURE_DISABLED` and
    `400 REF_MISMATCH` each render their intended state, and the 503s show a countdown.
13. `npm run typecheck && npm test && npm run build && npm run lint:pkg` pass; new components are axe-clean;
    `npm run fixtures` writes the names of §6.4 and the copied files validate in gramene-swagger's `contract.test.js`.

**Manual browser checks (the ones that matter here)**

14. Playground `?api=mock`: every new page renders; keyboard navigation works across the allele matrix (arrows, Home/End,
    Enter to open a cell) and the variant table (up/down/enter); focus is visible throughout.
15. Playground `?api=live` against a swagger that serves `/primers`: pick a variant, design, submit a **3-genome** check
    and watch a real job go `queued → running → done` with the matrix filling in genome by genome; then reload the page
    and confirm the job re-attaches from `jobId`.
16. Light, dark and `auto` themes; the page never scrolls horizontally; wide tables scroll inside their own container.
17. Copy buttons copy `target_seq` and `order_seq` separately; the exports download and open in a spreadsheet with the
    columns intact.
18. gramene-search with `PRIMERS_API=<v11 swagger>`: the Primers tab offers the Genotyping mode, a design plus check
    survives switching detail tabs and re-selecting the gene, and a saved view restores the genotyping slice.
19. An expired job ("Restore with an expired job") shows "Re-run check", not an error.

---

## 9. Known mismatches in the surrounding code

| Where | Mismatch | What to do |
| --- | --- | --- |
| `gramene-search/src/demo.js:146` | the sorghum site points `ensemblRest` at `https://data.gramene.org/pansite-ensembl-108`, while this API reads variants from **Ensembl REST release 115** (`https://data.gramene.org/pansite-ensembl-115`) | harmless for this feature — variants only ever come through swagger — but do **not** use `config.ensemblRest` for anything variant-related, and expect the Location tab's Dalliance browser (`dallianceBrowser.js:60,66`) to show release-108 annotation beside your release-115 variants. The same 108 URL appears at `demo.js:98` and `:249` |
| `gramene-search/src/demo.js:147` | `grameneData: 'https://data.sorghumbase.org/sorghum_v10b'`, so the demo loads **gene documents from v10b** while the primers API is **v11** | `PRIMERS_API` already points the Primers tab at its own base (`demo.js:167`), so the two APIs coexist, but the `gene` prop comes from v10b: gene ids, `system_name` and especially `taxon_id` can differ (`src/bundles/docs.js:11` records `sorghum_v10b`'s 4558001 vs `sorghum_v11`'s 4558006). Resolve genomes through `listGenomes(gene.system_name)` and treat a missing/mismatched genome as "choose a genome", never as an error; for release-sensitive work point `grameneData` at the v11 API |
| `gramene-search/src/components/results/details/Primers.js:22` | `MODES` is `['gene','transcript','region','sequence']` | add `'genotyping'`; no new props are needed — the component already passes `systemName={gene.system_name}` (`:55`) and the designer reads `variation` from `listGenomes`. `uiViewState.js` and `viewSnapshot.js` need no change (`snapshotPrimers` only strips `sequence`), but keep the genotyping slice small |
| `gramene-primers/src/components/ParamsPanel.tsx:188` | the preset radio list is hard-coded `['pcr','qpcr']` | parameterize it for `['kasp','as_pcr']` |
| `gramene-primers/src/types.ts:656` | `ResultsTab` has no allele tab | add the genotyping tab ids (`'alleles'`, `'order'`) as a separate union used by `GenotypingState.view.tab`, so `PrimerDesignerState.view.resultsTab` keeps validating in old builds |
| spec §8 vs reality | spec §8.8 promised a **partial** job fixture and a live endpoint | no partial capture exists (only the finished job); simulate partials in the mock. The pm2 repoint that gives the UI session a live genotyping endpoint is a separate, user-gated step |

---

## Appendix — numbers worth pinning in tests

| Value | Number |
| --- | --- |
| Design budget per request | 54 Primer3 runs, 272 `ntthal` calls, 8 candidates scored per orientation, 20 pairs per Primer3 run |
| Hard floors | AS Tm ≥ 52 °C; derived-primer GC ≥ 15 % |
| Structure thresholds | warn 47 °C, high 55 °C; `AS_TM_IMBALANCE` warn 1.0 °C, high 2.0 °C; `common_minus_as` window −1…+3 °C |
| Variant limits | window ≤ 50,000 bp; `limit` default 2000, max 5000; allele ≤ 50 nt; shift ≤ 1,000 bp |
| Check caps | 5 sets, 10 pairs, 13 primers (server-built request); 20 distinct primers, 6,000 CPU-s (endpoint) |
| Off-locus threshold | ≤ 2 mismatches per primer (`check.genotype_offlocus_max_mismatches`) |
| Ortholog test | gap-compressed identity ≥ 95 %, product size within ±20 %; paralogs counted at ≥ 80 % |
| Example job | estimate 273 CPU-s, measured 176.8 CPU-s (≈ 1.5×), 36.5 s wall, 12 genome tasks, `WEAK_OFF_TARGETS` 28 |
| Set keys | rs871475760 S1 `f9df650ad116`, S2 `1accc54c262d`; insertion S1 `53942cb55348`; deletion S1 `a5277232d8ab`, S2 `cb3ef66afd37` |
