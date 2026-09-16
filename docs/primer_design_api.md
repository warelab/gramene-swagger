# Primer Design API

PCR and qPCR primer design on the assemblies this site serves, Primer-BLAST-style specificity checks against the
reference genome (and its transcriptome), and pan-genome coverage checks across every assembly of the same species.
This documents the `/primers*` endpoints served by this gramene-swagger instance.

- **Base path:** `https://data.sorghumbase.org/sorghum_v11` (local dev: `http://localhost:50111/sorghum_v11`)
- **Controller:** `api/controllers/primers.js`; helpers under `api/helpers/primers/` (check algorithm in `check/`,
  job store and worker in `jobs/`)
- **Engines:** Primer3 2.6.1 (`primer3_core`, spawned per request) and BLAST+ 2.13.0 (`blastn`, `blastdbcmd`)
- **Data:** bgzipped FASTA and BLAST databases under `/scratch/olson/fasta/<system_name>/`; gene models from mongo
- **Design spec:** `docs/primer_design_spec.md` (the approved plan's Overrides take precedence over it)

| Endpoint | Purpose | Auth |
| --- | --- | --- |
| `POST /primers/design` | Design primers (synchronous, ≤ 45 s) | none |
| `GET /primers/genomes?system_name=` | Same-species genomes and what each can do | none |
| `POST /primers/check` | Queue a specificity / pan-genome check | none |
| `GET /primers/check/{job_id}` | Poll a check: progress, partial and final results | none |

> **POST bodies must be JSON with `Content-Type: application/json`**, otherwise swagger rejects them with
> `400` and `errors[].code` `INVALID_CONTENT_TYPE` before the handler runs. Bodies are limited to **100 kb** (`413`).

- Every `/primers` response carries `Cache-Control: no-store` (errors, validator 400s and 413s included).
- CORS is `*`; `Retry-After` is exposed to browsers (`Access-Control-Expose-Headers: Retry-After`).
- Unknown paths under `/primers` are a JSON `404 NOT_FOUND`; a wrong method is a JSON `405 METHOD_NOT_ALLOWED` with an
  `Allow` header — never express's HTML page.
- No authentication (decided: anonymous, like `/gene_lists/validate`). The host is protected by concurrency caps and
  the check cost guard instead.

---

## Quick start

```bash
BASE=https://data.sorghumbase.org/sorghum_v11
post() { curl -s -w '\n%{http_code}\n' -H 'Content-Type: application/json' -d "$2" "$BASE$1"; }

post /primers/design '{"mode":"transcript","gene_id":"SORBI_3001G000200"}'            # qPCR, junction-spanning
post /primers/design '{"mode":"gene","gene_id":"SORBI_3001G000200","flank_up":200,"flank_down":100}'
curl -s "$BASE/primers/genomes?system_name=sorghum_bicolor"

post /primers/check '{"system_name":"sorghum_bicolor","mode":"gene","gene_id":"SORBI_3004G087700",
  "pairs":[{"id":"P2","left":"GGACAGCTCCACAACATATCAG","right":"GGACATTTGAAGCCCATGGCC",
            "expected":{"region":"4","start":7423537,"end":7423746}}]}'                   # -> 202 {job_id,...}
curl -s "$BASE/primers/check/<job_id>"                                                   # poll until done|error
```

---

## Errors

Two shapes exist, depending on which layer rejected the request.

**Swagger validator** (types, required fields, `additionalProperties: false`, patterns, maxima, Content-Type):

```jsonc
{ "message": "Validation errors",
  "errors": [ { "code": "INVALID_REQUEST_PARAMETER", "in": "body", "name": "body",
                "errors": [ { "code": "OBJECT_ADDITIONAL_PROPERTIES", "params": ["bogus"], "message": "..." } ] } ] }
```

Inner codes you will meet: `INVALID_CONTENT_TYPE`, `OBJECT_ADDITIONAL_PROPERTIES`, `OBJECT_MISSING_REQUIRED_PROPERTY`,
`PATTERN`, `ENUM_MISMATCH`, `INVALID_TYPE`, `MINIMUM`/`MAXIMUM`, `MAX_LENGTH`, `ARRAY_LENGTH_SHORT`/`_LONG`, `REQUIRED`.

**Handler** (everything that needs the request's meaning, the catalog or the filesystem):

```jsonc
{ "message": "unknown gene NOPE", "code": "UNKNOWN_GENE", "details": { "gene_id": "NOPE" } }
```

- `code` is UPPER_SNAKE; `details` is always an object.
- **Every `503` carries `details.retry_after_s` and a matching `Retry-After` header.** Honour it.
- A `500` is `INTERNAL` with the message `internal error`; stack traces and filesystem paths never reach the client
  (paths in messages are reduced to basenames).
- **Send JSON types exactly and omit optional fields instead of sending `null`**: `null` fails swagger's type check.

| Status | Codes |
| --- | --- |
| 400 | validator errors; `INVALID_REQUEST` (`details.field`), `INVALID_PARAMS` (design: `details.param`; check: `details.field`), `INVALID_SEQUENCE`, `TEMPLATE_TOO_LONG`, `REGION_OUT_OF_BOUNDS`, `INTERVAL_OUT_OF_BOUNDS`, `SYSTEM_NAME_MISMATCH`, `PRIMER3_INPUT_ERROR` (`details.primer3_error`), `DUPLICATE_PAIR_ID`, `TOO_MANY_PRIMERS`, `PRODUCT_TOO_LONG_TO_CHECK`, `GENOME_NOT_CHECKABLE` |
| 404 | `UNKNOWN_GENE`, `UNKNOWN_TRANSCRIPT`, `UNKNOWN_GENOME`, `UNKNOWN_REGION`, `UNKNOWN_JOB`, `NOT_FOUND` |
| 405 | `METHOD_NOT_ALLOWED` (`details.allowed_methods`) |
| 413 | body over 100 kb (body-parser shape, below; no `code`) |
| 422 | `NO_SEQUENCE`, `AMBIGUOUS_ASSEMBLY` (`details.candidates`), `NO_BLASTDB`, `JOB_TOO_LARGE` |
| 500 | `PRIMER3_FAILED`, `GENE_STRUCTURE_MISMATCH`, `INTERNAL` |
| 503 | `BUSY` (5 s), `PRIMER3_UNAVAILABLE`, `MONGO_UNAVAILABLE` (30 s), `QUEUE_FULL` (60 s), `JOB_STORE_UNAVAILABLE` (5 s), `FEATURE_DISABLED` (300 s) |
| 504 | `DEADLINE_EXCEEDED` |

**Body parser.** An oversized or malformed JSON body is rejected while swagger parses it, before validation, and
comes back in body-parser's own shape — check `type`, there is no `code`:

```jsonc
// 413
{ "message": "request entity too large", "expected": 110034, "length": 110034, "limit": 102400, "type": "entity.too.large" }
// 400
{ "message": "Unexpected end of JSON input", "expose": true, "status": 400, "body": "{\"mode\":", "type": "entity.parse.failed" }
```

---

## Coordinate conventions

These are part of the contract; every number in requests and responses follows them.

- **Template coordinates** are 1-based inclusive, in template orientation (Primer3 runs with
  `PRIMER_FIRST_BASE_INDEX=1`). **Intervals** (`target`, `included`, `excluded`, `mask` runs) are `[start, length]`.
- **Primer footprints:** left primer `[start, end]` reads 5′→3′ along the template; the right primer's `end` is its
  5′ base on the template (Primer3's `PRIMER_RIGHT_i=pos,len` names the rightmost base, so the footprint is
  `[pos-len+1, pos]`). `product_size = right.end - left.start + 1`. Primer sequences are always **uppercase** 5′→3′.
- **Junction value `p`** is the boundary between cDNA bases `p` and `p+1`. A left primer `[a,b]` spans `p` when
  `p-a+1 ≥ min_5_prime_overlap` and `b-p ≥ min_3_prime_overlap`; a right primer when `p-a+1 ≥ min_3` and `b-p ≥ min_5`.
- **Genomic strand of a primer** is the strand whose 5′→3′ sequence equals the primer: the left primer is on
  `template.strand`, the right primer on the opposite strand.
- **`genomic.blocks`** are ascending `{start, end}` pieces; a cDNA primer spanning an exon–exon junction has two. Fetch
  the blocks on `genomic.strand` and join them 5′→3′ (for strand −1, in descending order) to get the primer back.
  Verified: SORBI_3001G000200.1 cDNA 856–875 → `1:13395–13400 + 13637–13650 (−)`.
- **`product.genomic`** is the envelope of both primers; `product.genomic_size` (transcript templates only) is its length.
- **Gene-relative position `p`** (mongo `gene_structure.exons`) maps to genomic as `strand == 1 ? start+p-1 : end-p+1`.
- **Check amplicons** (`start`, `end`) are the genomic 5′ ends of the forward-facing and reverse-facing primer sites,
  which is what `product.genomic` gives for a designed pair — pass those as `expected`.

---

## `POST /primers/design`

Synchronous. One Primer3 run per request, on one template.

### Modes

| `mode` | Required | Optional | Template |
| --- | --- | --- | --- |
| `gene` | `gene_id` | `system_name` (must match the gene), `transcript_id` (feature overlay), `flank_up`, `flank_down` | Genomic gene span ± flanks (clamped to the region), in transcription orientation |
| `transcript` | `gene_id` | `transcript_id` (default: canonical), `junction_spanning` (default `true`) | Spliced cDNA of the transcript |
| `region` | `system_name`, `region` | | Genomic region on the given strand (`strand` defaults to 1) |
| `sequence` | `sequence` | `system_name` (resolved; used for repeat masking and later checks) | Cleaned pasted sequence |

### Request fields

| field | type | notes |
| --- | --- | --- |
| `mode` | `gene` \| `transcript` \| `region` \| `sequence` | required |
| `gene_id`, `transcript_id` | string ≤ 255 | equality lookups only |
| `system_name` | `^[a-z0-9_]+$`, ≤ 128 | must be in the maps catalog |
| `region` | `{region ≤ 255, start ≥ 1, end ≥ 1, strand 1\|-1}` | `end - start + 1 ≤ 50,000` |
| `sequence` | string ≤ 60,000 chars | FASTA headers allowed; see [Sequence mode](#sequence-mode) |
| `flank_up`, `flank_down` | int 0–10,000 | gene mode |
| `target`, `included` | `[start, length]` | template coordinates; must lie inside the template |
| `excluded` | ≤ 50 × `[start, length]` | |
| `avoid_repeats` | bool | see [Repeat avoidance](#repeat-avoidance) |
| `repeat_mask_mode` | `n_mask` (default) \| `three_prime` | only with `avoid_repeats` |
| `junction_spanning` | bool | transcript mode; default `true` |
| `template_only` | bool | return the template (with its mask) and no pairs; Primer3 is not run, and `product_size_ranges` are not checked against the template length (every other rule still applies) |
| `params` | object | closed set, below; unknown keys are rejected |

### Params, presets and Primer3 tags

Unset params come from the preset: `pcr` for gene/region/sequence, `qpcr` for transcript mode. The effective values
are echoed in `settings.params`; the client only needs to send what it changed.

| param | limits | Primer3 tag | `pcr` | `qpcr` |
| --- | --- | --- | --- | --- |
| `opt_size` / `min_size` / `max_size` | int 15–36 | `PRIMER_OPT_SIZE` / `_MIN_SIZE` / `_MAX_SIZE` | 20 / 18 / 25 | 20 / 18 / 24 |
| `opt_tm` / `min_tm` / `max_tm` | 30–90 | `PRIMER_OPT_TM` / `_MIN_TM` / `_MAX_TM` | 60 / 57 / 63 | 60 / 58 / 62 |
| `opt_gc` / `min_gc` / `max_gc` | 0–100 | `PRIMER_OPT_GC_PERCENT` / `_MIN_GC` / `_MAX_GC` | – / 30 / 70 | – / 35 / 65 |
| `max_tm_diff` | 0–30 | `PRIMER_PAIR_MAX_DIFF_TM` | 3 | 2 |
| `max_poly_x` | int 0–10 | `PRIMER_MAX_POLY_X` | 4 | 4 |
| `gc_clamp` | int 0–5 | `PRIMER_GC_CLAMP` | – | – |
| `max_end_stability` | 0–100 | `PRIMER_MAX_END_STABILITY` | – | – |
| `max_ns` | int 0–5 | `PRIMER_MAX_NS_ACCEPTED` | – | – |
| `salt_monovalent`, `salt_divalent`, `dntp_conc`, `dna_conc` | 0–1000, 0–100, 0–100, 0–10000 | `PRIMER_SALT_MONOVALENT`, `_SALT_DIVALENT`, `_DNTP_CONC`, `_DNA_CONC` | – | – |
| `num_return` | int 1–20 | `PRIMER_NUM_RETURN` | 5 | 5 |
| `min_3_prime_overlap_of_junction`, `min_5_prime_overlap_of_junction` | int 1–20 | `PRIMER_MIN_3/5_PRIME_OVERLAP_OF_JUNCTION` (+ the `INTERNAL` twins) | – | 4, 7 |
| `product_size_ranges` | 1–10 × `[a, b]`, ints 20–50,000 | `PRIMER_PRODUCT_SIZE_RANGE` | `[[100,1000]]` | `[[70,150]]` |

"–" means Primer3's own default. Always sent by the server: `PRIMER_TASK=generic`, `PRIMER_EXPLAIN_FLAG=1`,
`PRIMER_LIBERAL_BASE=1`, `PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT=1`, `PRIMER_PRODUCT_MIN_TM=0`/`MAX_TM=150` (forces
`product_tm` to be printed without constraining it), `PRIMER_FIRST_BASE_INDEX=1`. Conditional: `SEQUENCE_TARGET`,
`SEQUENCE_INCLUDED_REGION`, `SEQUENCE_EXCLUDED_REGION`, `SEQUENCE_OVERLAP_JUNCTION_LIST` (transcript mode with
`junction_spanning`), and the masking tags below.

**Cross-field rules** (`400 INVALID_PARAMS`), mirrored by the UI:

- `min ≤ opt ≤ max` for size, Tm and GC; every product range has `a < b`.
- `max_size` ≤ the smallest product-range start; `gc_clamp` ≤ `min_size`.
- Transcript mode with `junction_spanning`: `min_5_prime_overlap_of_junction` and `min_3_prime_overlap_of_junction`
  (defaults 7 and 4) must each be ≤ `floor(max_size / 2)` — Primer3 otherwise aborts the whole run.
- Intervals must lie inside the template (`400 INTERVAL_OUT_OF_BOUNDS`), and at least one product range must start
  within the template length.

### Sequence mode

1. Lines starting with `>` are dropped; whitespace and digits are removed. If more than one FASTA record contains
   sequence, the records are joined into one template and warning `MULTIPLE_RECORDS` is added (primers may span the
   joins).
2. What remains must be IUPAC nucleotides (`ACGTRYKMSWBDHVN`, any case), 20–50,000 nt, else `400 INVALID_SEQUENCE`
   (longer: `TEMPLATE_TOO_LONG`).
3. Ambiguity codes other than N become N, with warning `IUPAC_CONVERTED`.
4. Case is kept only with `avoid_repeats`, where the user's lowercase **is** the mask (`mask_source: user_lowercase`).
5. Every `genomic` field in the response is `null`.

### Examples

```jsonc
{"mode":"gene","gene_id":"SORBI_3001G000200","flank_up":200,"flank_down":100,"params":{"product_size_ranges":[[300,800]]}}
{"mode":"transcript","gene_id":"SORBI_3001G000200","junction_spanning":true}
{"mode":"region","system_name":"sorghum_bicolor","region":{"region":"1","start":11080,"end":15099,"strand":-1},
 "target":[499,50],"included":[100,3000],"excluded":[[1000,40]],
 "params":{"min_size":19,"max_size":22,"min_tm":58,"max_tm":61,"min_gc":40,"max_gc":60,"max_tm_diff":2,"product_size_ranges":[[200,300]]}}
{"mode":"sequence","sequence":">amp\nATGGCCRYT...","system_name":"sorghum_bicolor","avoid_repeats":true}
```

### Response (`200`)

Coordinates and sequences below are real (SORBI_3001G000200.1); Tm and penalties are illustrative.

```jsonc
{
  "template": {
    "mode": "transcript", "system_name": "sorghum_bicolor", "gene_id": "SORBI_3001G000200", "transcript_id": "SORBI_3001G000200.1",
    "region": "1", "start": 11180, "end": 14899, "strand": -1,       // transcript: genomic envelope of its exons
    "length": 1982, "seq": "…",                                         // UPPERCASE, unmasked
    "masked": false, "mask_source": null, "mask": [], "masked_fraction": 0,
    "features": {
      "gene": null,                                                     // gene mode: {start, end} in template coordinates
      "exons": [ { "id": "EER93047-1", "start": 1, "end": 397, "genomic": { "start": 14503, "end": 14899 } }, "…" ],
      "cds": { "start": 299, "end": 1651 },                              // null for non-coding transcripts
      "junctions": [397, 493, 597, 869, 992, 1081, 1187, 1285, 1546, 1630]
    }
  },
  "pairs": [ {
    "rank": 0, "penalty": 0.41, "product_size": 124, "product_tm": 81.2, "compl_any_th": 0.0, "compl_end_th": 0.0,
    "left":  { "seq": "ATTACATCAAATAGGCCTTG", "start": 856, "end": 875, "len": 20, "tm": 59.1, "gc": 35.0,
               "self_any_th": 0, "self_end_th": 0, "hairpin_th": 0, "end_stability": 3.1, "penalty": 0.2,
               "junction": { "position": 869, "overlap_5p": 14, "overlap_3p": 6 },
               "genomic": { "region": "1", "start": 13395, "end": 13650, "strand": -1,
                            "blocks": [ { "start": 13395, "end": 13400 }, { "start": 13637, "end": 13650 } ] } },
    "right": { "seq": "AACTTCTTTGTCGATCCATG", "start": 960, "end": 979, "len": 20, "…": "…", "junction": null,
               "genomic": { "region": "1", "start": 13291, "end": 13310, "strand": 1, "blocks": [ { "start": 13291, "end": 13310 } ] } },
    "product": { "start": 856, "end": 979, "genomic": { "region": "1", "start": 13291, "end": 13650, "strand": -1 }, "genomic_size": 360 }
  } ],
  "explain": { "left": { "raw": "considered 9120, low tm 3280, ok 1022", "considered": 9120, "low tm": 3280, "ok": 1022 },
               "right": { "raw": "…" }, "pair": { "raw": "considered 91, ok 5" } },
  "settings": { "preset": "qpcr", "junction_spanning": true, "avoid_repeats": false, "repeat_mask_mode": null,
                "params": { "opt_size": 20, "…": "…" } },
  "engine": { "primer3": "2.6.1" },
  "warnings": []
}
```

- `pairs` may be empty: that is a `200` with warning `NO_PAIRS`; read `explain` to see why candidates were rejected.
- `product_tm` is `null` when Primer3 did not report it. `explain` is `null` with `template_only`.
- With `repeat_mask_mode: n_mask`, the `product_tm` of every product covering masked bases is recomputed on the
  **unmasked** template with Primer3's `long_seq_tm` and the same salts (Primer3 itself would count the `N` bases as
  A/T and report a Tm several degrees too low).
- `settings.params.max_ns` is forced to `0` when `avoid_repeats` uses `n_mask`.

**Warnings** (`[{code, message}]`, in this order): `ASSEMBLY_MISMATCH`, `AMBIGUOUS_ASSEMBLY` (from the assembly
resolver), `MULTIPLE_RECORDS` (sequence mode: several FASTA records joined), `IUPAC_CONVERTED`, `SINGLE_EXON_TRANSCRIPT` (junction spanning requested on a one-exon transcript: designed
without it), `JUNCTIONS_TRUNCATED` (more than 200 junctions), `REPEAT_MASK_FAILED`, `NO_REPEAT_MASK`,
`BLAST_DEPTH_MASK`, `MOSTLY_REPEAT` (more than 80% masked), `PRIMER3_WARNING`, `NO_PAIRS`.

### Limits

| Limit | Value |
| --- | --- |
| Template length | ≤ 50,000 nt (`TEMPLATE_TOO_LONG`) |
| Flank | ≤ 10,000 per side |
| Internal genomic fetch | ≤ 2 Mb |
| `num_return` / product ranges / `excluded` | ≤ 20 / ≤ 10 / ≤ 50 |
| Concurrency per API process | 4 running + 16 waiting, 10 s wait, then `503 BUSY` |
| Deadline | 45 s from the start of work (template + mask + Primer3); Primer3 itself ≤ 30 s → `504 DEADLINE_EXCEEDED` |
| Primer3 output | ≤ 2 MB (`500 PRIMER3_FAILED`) |

A client that disconnects aborts its design (Primer3 and any megablast are killed). Measured through the dev
instance's HTTP API on squam (warm, 2026-09-12): 30–100 ms for a gene, transcript, region or pasted-sequence design,
~0.13 s with a real soft-mask (sorghum_tx436pac, 17 kb), ~0.35 s with a BLAST depth mask; a 50 kb template with 20
pairs takes ~0.5 s in Primer3.

---

## Repeat avoidance

With `avoid_repeats: true` the server builds a mask for the template and keeps primers out of it:

- `repeat_mask_mode: n_mask` (default): masked bases become `N` and `PRIMER_MAX_NS_ACCEPTED=0`, so no primer can
  overlap them.
- `repeat_mask_mode: three_prime`: masked bases are lowercased and `PRIMER_LOWERCASE_MASKING=1`, which only keeps
  primer **3′ ends** out of the mask.

`template.mask_source` says where the mask came from:

| `mask_source` | When | What it is |
| --- | --- | --- |
| `softmask` | the assembly's `repeat_masking` is `soft_masked` | The assembly's own soft-masked FASTA (`dna_sm`) at the same coordinates — the real RepeatMasker annotation |
| `blast_depth` | otherwise, when the genome has a BLAST DB | Megablast of the template (20 kb chunks; exons ± 100 bp in transcript mode) against its own genome: bases covered by ≥ 3 HSPs of ≥ 50 bp at ≥ 85% identity. Warning `BLAST_DEPTH_MASK` |
| `user_lowercase` | sequence mode with lowercase input | The caller's own lowercase |
| `null` + `NO_REPEAT_MASK` | nothing is available | Unmasked design |

> **Honest note on masking quality.** Only **four sorghum assemblies have a real soft-mask:** `sorghum_rio`,
> `sorghum_tx2783pac`, `sorghum_tx430nano` and `sorghum_tx436pac`. Of the other sorghum assemblies, 76 ship a
> `dna_sm` file that is just an unmasked copy (`repeat_masking: unmasked_copy`, detected by sampling it for lowercase —
> file size alone is not a reliable signal) and 40 have none (`absent`). Everything except those four gets the
> **`blast_depth` mask, which masks copy number, not repeat families**: it also masks multi-copy gene families and
> recent duplications (37% of the SORBI_3004G087700.3 cDNA is masked), and misses low-copy repeats (against
> RepeatMasker on tx436pac: precision 0.91, recall 0.60). If megablast fails or runs out of its 20 s budget the design
> proceeds unmasked with `REPEAT_MASK_FAILED`. Specificity is what the check endpoint is for; masking only lowers the
> odds of a repetitive primer, it does not prove uniqueness.

`sorghum_pi534133`, `sorghum_pi576434` and `sorghum_pi656057` do have a genuinely soft-masked `dna_sm`, but only as an
uncompressed, unindexed file, so they currently resolve as `absent`; bgzipping and indexing those files would make them
`soft_masked`.

---

## `GET /primers/genomes`

The query genome and every other genome of its species, with what each one supports. The species is the first
taxonomy ancestor of rank `species` (never `taxon_id / 1000`: rice genomes hang under 39946, 1736656 and 1736659, all
under species 4530).

| param | in | type | required |
| --- | --- | --- | --- |
| `system_name` | query | `^[a-z0-9_]+$`, ≤ 128 | yes |

```jsonc
{ "system_name": "sorghum_bicolor",
  "species": { "taxon_id": 4558, "name": "Sorghum bicolor" },
  "counts": { "total": 120, "with_blastdb": 120, "with_cdna_blastdb": 120 },
  "genomes": [
    { "system_name": "sorghum_bicolor", "display_name": "Sb bicolor BTx623 v3", "taxon_id": 4558006, "map_id": "GCA_000003195.3",
      "is_query": true, "has_sequence": true, "has_blastdb": true, "has_cdna_blastdb": true,
      "repeat_masking": "unmasked_copy", "total_bases": 708735318, "warnings": [] },
    { "system_name": "sorghum_rio", "…": "…", "repeat_masking": "soft_masked" } ] }
```

- The query genome is first, then the rest by `display_name`. No filesystem paths are ever returned.
- A genome whose assembly cannot be resolved is still listed, with `has_*` false, `repeat_masking: absent`,
  `total_bases: null` and a `warnings` entry (e.g. `AMBIGUOUS_ASSEMBLY`).
- `ASSEMBLY_MISMATCH` means fewer than 80% of the map's **real** regions were found (by name and length) in the FASTA
  index. Synthetic bins (`UNANCHORED`, `UNPLACED`, `UNASSIGNED`, `UNLOCALIZED`, any case) are excluded from both the
  count and the total, and a map with only bins (e.g. `sorghum_tx430nano`, `selaginella_moellendorffii`) gets no
  warning. The message reads `only X of Y map regions match …` with Y = real regions.
- `404 UNKNOWN_GENOME`; `503 MONGO_UNAVAILABLE` only if the catalog has never loaded (a stale catalog keeps being
  served while mongo is down). Results are cached for 10 minutes.

---

## `POST /primers/check`

Queues a check of up to 10 primer pairs. Specificity against the reference genome always runs (plus its
transcriptome in transcript mode); `pangenome` adds every other checkable assembly of the species.

### Request (frozen contract)

```jsonc
{
  "system_name": "sorghum_bicolor",               // required; the reference genome
  "mode": "gene",                                 // gene | transcript | region | sequence; default region
  "gene_id": "SORBI_3004G087700",                 // required in gene and transcript modes
  "transcript_id": "SORBI_3004G087700.1",         // optional; defaults to the canonical transcript
  "checks": ["specificity", "pangenome"],         // default [specificity]; specificity always runs
  "genomes": ["sorghum_353", "sorghum_leoti"],    // pangenome only; ≤ 150; omitted = all checkable same-species genomes
  "params": { "max_product_size": 4000, "ignore_mismatches": 6, "max_amplifying_mismatches": 3, "min_total_mismatches": 2,
              "min_3p_mismatches": 2, "three_prime_window": 5, "include_unlikely": false, "repeat_site_threshold": 5 },
  "pairs": [ { "id": "P3",                        // ^[A-Za-z0-9_.:-]+$, ≤ 64, unique
               "left": "GATATCAGTGGAATCATAAGACCG", "right": "CATCGATATCAGGATCTGGCTT",   // ^[ACGTacgt]{15,36}$
               "expected": { "region": "4", "start": 7422482, "end": 7423061 } } ]     // gene/region modes only
}
```

| param | range | default | meaning |
| --- | --- | --- | --- |
| `max_product_size` | 50–10,000 | 4000 | largest product called (raised automatically, below) |
| `max_amplifying_mismatches` | 0–5 | 3 | a product amplifies only if **each** primer site has at most this many edits; products over it are `unlikely`. Must be below `ignore_mismatches` |
| `ignore_mismatches` | 3–6 | 6 | a primer site with this many edits cannot prime (product ignored, not even `unlikely`) |
| `min_total_mismatches` | 0–6 | 2 | with `min_3p_mismatches`: a site this bad is blocked (product `unlikely`) |
| `min_3p_mismatches` | 1–5 | 2 | mismatches required inside the 3′ window for blocking |
| `three_prime_window` | 3–10 | 5 | bases from the 3′ end counted as the 3′ window |
| `include_unlikely` | bool | false | also list `unlikely` products |
| `repeat_site_threshold` | 1–100 | 5 | a primer with more near-perfect genome sites is flagged `repetitive` |

**Handler rules**

- ≤ 10 pairs and ≤ 20 distinct primers (`TOO_MANY_PRIMERS`); duplicate pair ids → `DUPLICATE_PAIR_ID`.
- `params.max_amplifying_mismatches` must be below `ignore_mismatches`: an explicit value that is not is
  `400 INVALID_PARAMS` with `details {field: "params.max_amplifying_mismatches", max_amplifying_mismatches,
  ignore_mismatches}`. When it is omitted, the default 3 is lowered to `ignore_mismatches − 1` if needed
  (`{"ignore_mismatches": 3}` checks with a cap of 2); the filled value is part of the job id and of `results.params`.
  A value outside 0–5 or not an integer is a validator `400`.
- Gene/transcript modes look the gene up: its `system_name` must equal `system_name` (`SYSTEM_NAME_MISMATCH`); an
  unknown transcript is `404 UNKNOWN_TRANSCRIPT`.
- `expected` is honoured only in gene and region modes; elsewhere it is dropped with warning `EXPECTED_IGNORED`.
- **`max_product_size` is raised automatically** to `min(10000, max(current, ceil(1.2 × largest expected size)))` with
  warning `MAX_PRODUCT_SIZE_RAISED`, so a 3.6 kb on-target is not silently missed. An expected product over 10 kb is
  `400 PRODUCT_TOO_LONG_TO_CHECK`.
- The reference needs a dna BLAST DB (and a cdna DB in transcript mode), else `422 NO_BLASTDB`.
- `genomes` must be in the catalog (`404 UNKNOWN_GENOME`), of the same species and have the needed DB, else
  `400 GENOME_NOT_CHECKABLE` with `details.genomes` and `details.reasons` (`other_species`, `no_blastdb`,
  `no_cdna_blastdb`, `ambiguous_assembly`, `assembly_unavailable`). `genomes` without the pangenome check is ignored
  (`GENOMES_IGNORED`).
- The cost guard runs last (below).

### Responses

```jsonc
// 202: newly queued (also when a job that ended in error is re-queued); 200: an identical job already exists
{ "job_id": "9f2c0a4be1d34c7a8e5f00112233aabb", "status": "queued", "kind": "pangenome", "queue_position": 0,
  "progress": { "done": 0, "total": 4, "stage": "queued", "running": [] },
  "estimate": { "cpu_s": 16 }, "created_at": "2026-09-12T20:01:02.000Z", "warnings": [] }
```

## `GET /primers/check/{job_id}`

`job_id` is 32 lowercase hex characters (else `400` from the validator). Poll it; back off (start at 1 s, ×1.5,
cap 10 s) and honour `retry_after_s` on a `503`.

```jsonc
{ "job_id": "…", "status": "running", "kind": "pangenome", "partial": true, "queue_position": null,
  "progress": { "done": 37, "total": 121, "stage": "pangenome", "running": ["sorghum_is929", "sorghum_rio"] },
  "created_at": "…", "started_at": "…", "finished_at": null, "attempts": 1,
  "request": { "…normalized request…": "" }, "warnings": [], "estimate": { "cpu_s": 1880 },
  "results": { "…partial results, same schema…": "" }, "error": null }

{ "job_id": "…", "status": "done",  "partial": false, "progress": { "done": 121, "total": 121, "stage": "done" }, "results": { "…": "" } }
{ "job_id": "…", "status": "error", "error": { "code": "REFERENCE_BLAST_FAILED", "message": "…" } }
```

`404 UNKNOWN_JOB` when the id is unknown or expired — show "results expired, re-run".

### Job lifecycle

```
POST ─► queued ─► running ─► done   (results kept 24 h)
                     │  └──► error  (kept 1 h; re-POSTing the same request re-queues it)
                     └─ worker stopped ─► queued again (front of the queue), then running with attempts+1
```

- **Deterministic ids.** `job_id = sha256(canonicalJSON({v, algorithm_version, normalized request, DB fingerprints}))`
  truncated to 32 hex. The normalized request has primers uppercased, `checks` sorted with `specificity`, `genomes`
  expanded and sorted, all params filled and the transcript defaulted, so key order, case and defaults do not matter.
  Identical requests share one job and its results; a rebuilt BLAST DB or FASTA changes the id.
- **No cancel** in v1 (ids are shared across identical requests).
- **TTLs:** queued/running jobs 6 h (refreshed while running), done 24 h, error 1 h; each site keeps at most 500
  finished jobs. Results are gzipped and at most 5 MB (else error `RESULT_TOO_LARGE`).
- **Queue limit:** 50 queued jobs per site (both queues) → `503 QUEUE_FULL` (retry after 60 s).
- **Host-wide caps** (shared by every site and dev instance on the host): ≤ 2 running jobs, of which ≤ 1 pan-genome;
  each worker runs ≤ 2 jobs. A specificity job uses 4 BLAST threads, a pan-genome job ≤ 8 single-thread BLAST
  processes, all at `nice 10`. Job timeout 30 min (`JOB_TIMEOUT`).
- **Restarts.** Stopping the worker (SIGTERM/SIGINT, e.g. `pm2 restart`) puts its running jobs back at the front of
  their queue within a second; clients see `running → queued → running → done` without re-posting. After a hard kill
  the next worker requeues them when it starts (and a running entry whose heartbeat is older than 60 s is reclaimed).
  A job interrupted twice ends with error `WORKER_LOST`.
- **Redis down:** `POST /primers/check` and `GET /primers/check/{id}` answer `503 JOB_STORE_UNAVAILABLE` within 3 s;
  design and genomes keep working. The worker runs nothing without a claimed slot, so checks fail closed.
- **Job error codes:** `REFERENCE_BLAST_FAILED` (BLAST failed twice, or timed out once: `details {target, cause}`;
  stderr tail with paths reduced to basenames), `NO_BLASTDB`, `JOB_TIMEOUT`, `RESULT_TOO_LARGE {bytes, limit}`,
  `WORKER_LOST`, `MONGO_UNAVAILABLE` (mongo failed during the run on the job's last attempt; see Operations),
  `CHECK_FAILED`.

### Cost guard

Estimated before queueing, in CPU-seconds: single-thread coefficients measured on sorghum, plus a factor for the
concurrent pan-genome searches.

```
cpu_s = ceil( P × [ ref_Gb × c(5)  +  (transcript ? 0.15 × c(5) : 0)
                    + f_pan × Σ over pan-genome genomes of (transcript ? 0.15 : total_bases/1e9) × c(6)
                    + genome_tasks × 0.6 ] )
c(5) = 5.2, c(6) = 2.2, c(7) = 1.2 CPU-s per primer·Gb        P = distinct primers
f_pan = 2.0 (check.pangenome_cpu_factor)                       unknown genome sizes count as 1 Gb
genome_tasks = 1 + (transcript ? 0 : number of pan-genome genomes)
```

- **Re-alignment term:** the last term is FASTA re-alignment in the worker, at 0.6 CPU-s per primer per genome task (measured 0.26–0.58 on sorghum). cDNA tasks re-align from BLAST's own alignment and are not charged.
- **Pan-genome factor:** `f_pan` covers pan-genome assemblies being searched by up to 8 concurrent single-thread `blastn` processes.
  - The full-panel P3 check (2 primers, 119 assemblies, 83.8 Gb) used **758 CPU-s**: 733 in `blastn` and 25 in the worker.
  - Without the factor the formula estimated 520: word size 6 measured 4.33 CPU-s per primer·Gb, about twice the single-thread 2.2.
  - With `f_pan = 2.0` that job estimates 889 CPU-s.
- **Limit:** `cpu_s > 6000` → `422 JOB_TOO_LARGE {estimate_cpu_s, limit}`.
- **UI:** the UI mirrors this formula for display.

| Workload (sorghum) | Estimate |
| --- | --- |
| Specificity, 10 primers | 43 CPU-s (37 BLAST + 6 re-alignment) → 15–20 s wall on 4 threads |
| + transcriptome (transcript mode) | + 8 CPU-s |
| Pan-genome, all 119 other assemblies (83.8 Gb), 2 primers | ≈ 890 CPU-s (measured 758 CPU-s, 108 s wall on 8 processes) |
| Pan-genome, all 119 other assemblies, 10 primers | ≈ 4,440 CPU-s (3,690 BLAST + 720 re-alignment + 37 reference) |
| Pan-genome, all 119 other assemblies, 20 primers | ≈ 8,890 CPU-s → **refused** (`JOB_TOO_LARGE`); up to 13 distinct primers fit under 6,000 |
| Pan-genome in transcript mode (119 cDNA DBs at 0.15 Gb), 10 primers | ≈ 840 CPU-s |

Sites whose lower bound already exceeds `max_amplifying_mismatches` are not re-aligned (step 4 below), which keeps
re-alignment to about 1–1.6 s per pan-genome assembly for the 2-primer presence/absence pair (it was 5–7 s when every
candidate site was re-aligned).

Measured on the dev instance with check algorithm version 2 (squam, 2026-09-13; wall time from submit to `done`, one
job at a time, `PRIMERS_GLOBAL_MAX_JOBS=1`):

| Check (sorghum_bicolor reference) | Estimate | Wall | `timings_ms` |
| --- | --- | --- | --- |
| P2, specificity, 2 primers | 9 CPU-s | 8 s | reference 5,671 |
| 5 pairs / 8 primers, specificity, `include_unlikely` | 35 CPU-s | 16 s | reference 13,849 |
| qPCR, transcript mode (genome + cDNA), 2 primers | 11 CPU-s | 10 s | reference 6,197; transcriptome 1,130 |
| P3, pan-genome over 3 assemblies | 22 CPU-s | 12 s | reference 4,258; 353 6,377; grassl 6,454; leoti 7,291 |
| SORBI_3001G046200 pair, pan-genome over 4 assemblies (max size 5,932) | 27 CPU-s | 16 s | reference 5,145; per genome 6,961–8,934 |

### How a check works

Check algorithm version **2** (`results.engine.algorithm_version`; a version change changes every job id, so results
of an older algorithm are never reused).

1. **BLAST** each distinct primer: `blastn -task blastn-short -reward 1 -penalty -1 -ungapped -evalue 30000
   -searchsp 15000000000 -dust no -soft_masking false -max_hsps 100000`, word size **5** for the reference genome
   and transcriptome, **6** for pan-genome assemblies; `-max_target_seqs` = max(5000, sequences in the DB). A failed
   BLAST is retried once after 2 s, except a `BLAST_TIMEOUT` (10 min), which is not repeated: the reference then fails
   the job with `REFERENCE_BLAST_FAILED` (`details.cause: "BLAST_TIMEOUT"`), a pan-genome genome gets status `error`
   with `error.code: "BLAST_TIMEOUT"`.
2. **Sites:** each hit is extended over its unaligned tails to a full primer footprint; hits that cannot be within
   5 edits are dropped. A primer with an indel usually comes back as two hits, one per side of the gap: a hit that
   fails that filter but is within 5 edits of its own part and aligns at least 11 bases is held, and after BLAST two
   held halves of the same primer on the same sequence and strand that line up (colinear, 5′ ends within 3 bp) are
   joined into one gapped site (edits = both halves plus the gap).
3. **Candidates:** a forward-facing and a reverse-facing site on the same sequence within the maximum product size
   form a product: `LR`/`RL` for a proper pair, `LL`/`RR` for a single primer priming both ways. A pair whose left
   and right primers are identical labels its products `LR`. Pan-genome searches allow
   `max(max_product_size, 1.5 × reference_size + 500)`. A product must be at least as long as its two primers
   together, and after re-alignment **products whose two primer footprints overlap are discarded**.
4. **Re-alignment:** candidate sites are re-aligned to their genome window by dynamic programming (a gap counts as a
   mismatch), giving `mm`, mismatches in the 3′ window and whether the 3′-terminal base mismatches. A site whose lower
   bound (`mm` of the hit + 1 per unaligned tail) already exceeds `max_amplifying_mismatches` cannot amplify and is
   not re-aligned: it keeps that bound as `mm` with `approx: true` and `mm_pos: null`, so its products are `unlikely`.
   cDNA hits use BLAST's own alignment (`approx: true` when tails had to be estimated). Published coordinates are
   clamped to the sequence.
5. **Classification** (Primer-BLAST 3′ rule plus a per-primer cap), per product, first match wins:

   | likelihood | rule |
   | --- | --- |
   | ignored (dropped) | either primer site has ≥ `ignore_mismatches` edits |
   | `unlikely` | either site has more than `max_amplifying_mismatches` edits (default 3) |
   | `unlikely` | either site is 3′-blocked: ≥ `min_total_mismatches` edits with ≥ `min_3p_mismatches` of them in the 3′ window |
   | `likely_weak` | a 3′-terminal mismatch |
   | `likely` | everything else |

   `likely` and `likely_weak` count as amplicons: **a product amplifies only if each primer has at most
   `max_amplifying_mismatches` edits and neither primer is 3′-blocked.** Products with a 4- or 5-edit primer are
   `unlikely`: counted in `unlikely_count`, listed only with `include_unlikely`, and usable as a pan-genome `nearest`.
   Mismatch positions (`left_mm_pos`, `right_mm_pos`) are distances from the primer's 3′ end (1 = terminal base),
   listed 5′→3′.
6. **Verdicts, pan-genome statuses and gene annotation** (below).

With these defaults the §10.4 fixtures on sorghum_bicolor give: P2 two off-targets (4:7437317–7437526 `LR` 0/0 and
5:66890615–66890824 `RL` 2/2); P1 two off-targets, while its 5:66891530–66892235 product (a 4-edit left site,
mm positions 20, 17, 14, 11) is `unlikely`; P3 `specific`, its 4:7436221–7436845 product (4-edit left site) being
`unlikely`; P3 across sorghum_353, grassl and leoti `single_mismatch` in all three; the SORBI_3001G046200 pair
`no_amplicon` in sorghum_353 and sorghum_is12661 (their best products need 4+ edits in a primer) and `single_perfect`
in ji2731 and leoti. Lower `max_amplifying_mismatches` (e.g. 1) for only near-perfect products; raising it to 4–5
also calls products that are unlikely to amplify.

**Whole sorghum panel** (P3 with `genomes` omitted: all 119 other assemblies, default params; squam, 2026-09-13):
- **Run:** 108 s wall on 8 BLAST processes, 758 CPU-s.
- **Statuses:** `single_perfect` 66, `single_mismatch` 51, `multiple` 1, `no_amplicon` 1; `amplifies` 118, `truncated` 0.
- **Primary product sizes:** 580 bp ×94, 581 ×23, 579 ×1.
- **The `multiple` genome**, sorghum_pi601552a8087c, has two perfect 580 bp products (4:7940574 `LR` and contig_5132:70605
  `RL`), most likely an assembly duplication.
- **The `no_amplicon` genome**, sorghum_pi536008, has only an `unlikely` L4/R4 639 bp product at 1:55838383 as its
  `nearest`.

### Verdicts

| mode | target | on-target | verdicts |
| --- | --- | --- | --- |
| gene, region with `expected` | genome | an `LR`/`RL` amplicon within ±2 bp of `expected` | `specific`, `off_targets`, `on_target_missing`, `truncated`, `error` |
| gene/region without `expected`, sequence | genome | exactly one perfect `LR`/`RL` amplicon → `on_target_inferred: true` | `specific`, `off_targets`, `unverified_target`, `truncated`, `error` |
| transcript | genome | none; amplicons inside the gene span go to `gdna_products` | `specific`, `off_targets`, `truncated`, `error` — never `on_target_missing` |
| transcript | cDNA (`results.transcriptome`) | amplicons on transcripts of `gene_id`, grouped by gene (isoforms collapsed) | `specific`, `off_targets`, `on_target_missing`, `truncated`, `error` |

- `off_targets` includes single-primer (`LL`/`RR`) products. `off_target_count` is exact; at most 100 are listed.
- `truncated` means a cap was hit (5,000 candidates per pair, 20,000 re-aligned sites per genome, a primer with too
  many sites, or BLAST's `-max_target_seqs`): the list may be incomplete.

### Pan-genome statuses

| status | meaning (genome target: distinct amplicons; cDNA target: distinct gene groups) |
| --- | --- |
| `single_perfect` | exactly one amplicon, no mismatches |
| `single_mismatch` | exactly one amplicon with mismatches; flags `mismatch_in_3p_window`, `terminal_mismatch` |
| `multiple` | two or more amplicons: `primary`, `others` (≤ 10), `other_amplicons` |
| `no_amplicon` | none; `nearest` is the best unlikely product, if any (it can be approximate: `approx: true`, mismatch counts are lower bounds) |
| `db_unavailable` | no usable BLAST DB (with `reason`) |
| `error` | BLAST failed twice, or timed out once (with `error {code, message}`, e.g. `BLAST_TIMEOUT`) |

- The **primary** amplicon is chosen by: overlapping an annotated ortholog of the design gene, then fewest mismatches,
  then size closest to the reference product; `size_delta = size - reference_size`.
- `ortholog_annotated` says whether the design gene has any annotated ortholog in that genome — with `no_amplicon` it
  separates "no ortholog annotated there" from "the primers fail there". `ortholog` is `null` in region/sequence modes.
- `summary` counts always add up to `genomes_total`, and `amplifies = single_perfect + single_mismatch + multiple`.
  While a job runs, partial documents cover the genomes finished so far.
- `summary.truncated` is the number of genomes with `truncated: true` (a site, candidate or re-alignment cap was hit).
  It is not a status and not part of the `genomes_total` sum: such a genome's status is a lower bound (a
  `no_amplicon` or `single_*` there may have missed products).
- In transcript mode the pan-genome search runs against each assembly's **annotated transcripts only**
  (`PANGENOME_TRANSCRIPT_MODELS_ONLY`): an unannotated copy cannot be found.

### Sensitivity guarantee

BLAST finds a substitution-only site only if it contains an exact match of at least the word size and scores about
11 or more. Each primer reports the number of mismatches that is **guaranteed** to be found:

`guaranteed_max_mismatches(L, W)` = the largest `k` with `ceil((L-k)/(k+1)) ≥ W` and `L - 2k ≥ 12`.

| primer length | W = 5 (reference) | W = 6 (pan-genome) |
| --- | --- | --- |
| 20 | 3 | 2 |
| 22 | 3 | 2 |
| 25 | 4 | 3 |

Sites with more mismatches are often found too (word size 5 found a 22-mer with 3 mismatches that word size 7 missed),
but not guaranteed. **Mid-primer insertions or deletions in short primers can be missed**: an indel is found when the
primer part on at least one side of it scores about 11 (both sides are then joined into one gapped site), so indels
near the middle of primers shorter than about 22 nt may be missed.

**The guarantee and the amplification cap can differ in the pan-genome.** A product amplifies with up to
`max_amplifying_mismatches` (default 3) mismatches per primer.
- **Reference search (word size 5):** guarantees 3 mismatches for 20–24 nt primers.
- **Pan-genome search (word size 6):** guarantees only 2 for 20–23 nt primers, and 3 from 24 nt.

A 3-mismatch product in another assembly can therefore be missed, and that genome reported as `no_amplicon`.

`results.sensitivity_note` is written for each request. It names the actual cap, the word sizes and the request's
primer lengths, and it adds the pan-genome sentence only when a pan-genome search runs. For P3 (22 and 24 nt) with
`checks: [specificity, pangenome]` and the default cap:

> Sites are detected only if they contain an exact match of at least the word size and an ungapped score of about 11
> or more. An indel is found when the primer part on at least one side of it scores that much (both sides are joined
> into one gapped site); indels near the middle of primers shorter than about 22 nt may be missed. A product amplifies
> only when each primer has at most 3 mismatches (max_amplifying_mismatches). The reference search (word size 5) finds
> every site with up to 3 mismatches for these 22-24 nt primers, which covers the cap. The pan-genome search (word
> size 6) finds every site with up to 2-3 mismatches for these 22-24 nt primers, so products with 3 mismatches per
> primer can be missed in other genomes for the shorter primers. See primers[].sensitivity.

### Results (`status: done`)

```jsonc
{
  "engine": { "algorithm_version": "2", "blast": "2.13.0",
              "reference": "blastn-short r1 p-1 ws5 ungapped e30000 searchsp1.5e10", "pangenome": "… ws6 …" },
  "params": { "max_product_size": 4000, "ignore_mismatches": 6, "max_amplifying_mismatches": 3, "…": "…" },   // all 8, defaults filled
  "reference": { "system_name": "sorghum_bicolor", "map_id": "GCA_000003195.3", "total_bases": 708735318 },
  "sensitivity_note": "Sites are detected only if …",
  "primers": { "GGACAGCTCCACAACATATCAG": { "len": 22, "near_perfect_sites": 2, "repetitive": false, "truncated": false,
               "sensitivity": { "reference": { "word_size": 5, "guaranteed_max_mismatches": 3 } } } },
  "specificity": { "target": "genome", "pairs": [ {
     "id": "P2", "verdict": "off_targets", "on_target_inferred": false, "truncated": false,
     "on_target": { "region": "4", "start": 7423537, "end": 7423746, "size": 210, "strand": 1, "orientation": "LR",
                    "likelihood": "likely", "left_mm": 0, "right_mm": 0, "left_3p_mm": 0, "right_3p_mm": 0,
                    "left_mm_pos": [], "right_mm_pos": [], "terminal_mismatch": false, "approx": false,
                    "genes": [ { "id": "SORBI_3004G087700", "name": "…", "biotype": "protein_coding", "strand": 1 } ] },
     "off_target_count": 2, "unlikely_count": 0,
     "off_targets": [
       { "region": "4", "start": 7437317, "end": 7437526, "orientation": "LR", "left_mm": 0, "right_mm": 0, "…": "…",
         "genes": [ { "id": "SORBI_3004G087800" } ] },
       { "region": "5", "start": 66890615, "end": 66890824, "orientation": "RL", "left_mm": 2, "right_mm": 2,
         "left_mm_pos": [10, 7], "right_mm_pos": [15, 10], "likelihood": "likely", "genes": [ { "id": "SORBI_3005G183900" } ] },
       "…" ],
     "gdna_products": [] } ] },
  "transcriptome": null,        // transcript mode: {target: "cdna", gene_id, transcript_id, pairs: [{id, verdict, on_target: group, off_targets: [group], …}]}
  "pangenome": { "target": "genome", "pairs": [ {
     "id": "P3", "reference_size": 580, "max_size": 4000,
     "summary": { "genomes_total": 3, "single_perfect": 0, "single_mismatch": 3, "multiple": 0, "no_amplicon": 0, "db_unavailable": 0, "error": 0, "amplifies": 3, "truncated": 0 },
     "genomes": [ { "system_name": "sorghum_353", "display_name": "…", "status": "single_mismatch", "ortholog_annotated": true, "truncated": false,
        "primary": { "region": "4", "start": 7499931, "end": 7500510, "size": 580, "size_delta": 0, "orientation": "LR",
                     "likelihood": "likely_weak", "left_mm": 1, "left_mm_pos": [1], "mismatch_in_3p_window": true,
                     "terminal_mismatch": true, "genes": [ { "id": "353.004G093500" } ], "ortholog": true },
        "other_amplicons": 0, "others": [], "nearest": null } ] } ] },
  "warnings": [ { "code": "…", "message": "…" } ],
  "timings_ms": { "reference": 11000, "sorghum_353": 1500, "total": 18000 }
}
```

A cDNA group is `{gene_id, orientation, isoforms: [{transcript_id, start, end, size, likelihood, left_mm, right_mm}],
size_min, size_max, likelihood, left_mm, right_mm, left_3p_mm, right_3p_mm, left_mm_pos, right_mm_pos,
terminal_mismatch, approx}`. `unlikely[]` appears on a pair only with `include_unlikely`.

Result warnings: `NO_FASTA_FOR_REALIGN` (mismatch counts are lower bounds, `approx: true`), `ANNOTATION_UNAVAILABLE`
(mongo was unreachable: `genes: null`, `ortholog: null`; the worker instead fails the attempt with `MONGO_UNAVAILABLE`
and retries it after a restart, so a finished job carries this warning only when run outside the worker),
`TRANSCRIPT_GENE_UNMAPPED`, `PANGENOME_TRANSCRIPT_MODELS_ONLY`,
and assembly warnings. Submit warnings: `EXPECTED_IGNORED`, `MAX_PRODUCT_SIZE_RAISED`, `GENOMES_IGNORED`, and the
reference assembly's `ASSEMBLY_MISMATCH` / `AMBIGUOUS_ASSEMBLY` prefixed with its system name.

---

## Operations

### Processes

| pm2 app | Entry | Role |
| --- | --- | --- |
| `sorghum_swagger11` | `app.js` | The API, including `/primers/design`, `/primers/genomes` and job submit/status |
| `sorghum_primers11` | `api/helpers/primers/jobs/worker_main.js` | **The check worker, a separate pm2 app.** `app.js` does not start it |

```bash
cd <checkout>
pm2 start api/helpers/primers/jobs/worker_main.js --name sorghum_primers11 --cwd <checkout> --exp-backoff-restart-delay 5000
pm2 save
```

- **Both processes log `primers site_key=<key>` at startup and the two values must match**, otherwise the worker
  watches a different Redis namespace than the API writes to:
  `pm2 logs sorghum_swagger11 --lines 500 --nostream | grep 'primers site_key='` and the same for `sorghum_primers11`.
  The key defaults to `<swagger.yaml basePath without '/'>:<mongo db>`, e.g. `sorghum_v11:sorghum11`; both processes read
  `basePath` from `api/swagger/swagger.yaml`.
- The worker holds a per-site lock in Redis (a second worker for the same site waits), exits with **75** when mongo is
  unavailable (gramene-mongodb-config never reconnects) and lets pm2 restart it, and on SIGTERM/SIGINT requeues its
  running jobs and exits within a second (pm2's default kill timeout is enough).
- Mongo failing **during** a check (no collection, a query error or a 15 s timeout, in any stage) also ends in exit
  75: the job fails internally with `MONGO_UNAVAILABLE`, is requeued at the front of its queue, and the restarted
  worker runs it again, so no check finishes with `ANNOTATION_UNAVAILABLE` because of an outage. A job that is still
  hit on its last attempt (`check.max_attempts`, 2) ends with error `MONGO_UNAVAILABLE` (1 h error TTL; a re-POST runs
  it again) instead of being requeued forever.
- **BLAST+ and Primer3 isolation.** `blastn`, `blastdbcmd` (worker and design-time megablast mask) and `primer3_core`
  run with `PATH=/usr/bin:/bin` and `NCBI_DONT_USE_NCBIRC=1`, in a private `0700` `mkdtemp` directory under
  `primers.tmp_dir` that is removed afterwards, never with `/tmp` (or any shared directory) as their working directory.
  BLAST+ would otherwise read a `.ncbirc` from its working directory, `$HOME`, `/etc` or the binary directory, so any
  local user could plant one in `/tmp`. If that directory cannot be created, nothing is spawned: the check's BLAST
  fails, and a design answers `500 PRIMER3_FAILED` (`503 PRIMER3_UNAVAILABLE` for the version probe).
- **Rollback:** `pm2 stop sorghum_primers11`, then `PRIMERS_ENABLED=0 pm2 restart sorghum_swagger11 --update-env`
  (every `/primers` endpoint answers `503 FEATURE_DISABLED`; the rest of the API is unaffected).
- **Apache:** the proxy timeout must be ≥ 60 s (design runs up to 45 s) and `Cache-Control: no-store` must pass
  through mod_cache.

### Configuration

The `primers:` block of `config/default.yaml` (read by `api/helpers/primers/config.js`, merged over built-in defaults)
holds every limit quoted above. Environment overrides (read by `config.js`; a malformed value is ignored with a
warning, never a startup crash):

| Variable | Config key | Notes |
| --- | --- | --- |
| `PRIMERS_ENABLED` | `enabled` | `1`/`true`/`yes`; anything else disables `/primers` |
| `PRIMERS_SITE_KEY` | `site_key` | Redis namespace; set the same value for the API and the worker |
| `PRIMER3_CORE`, `BLASTN`, `BLASTDBCMD` | `primer3_core`, `blastn`, `blastdbcmd` | binaries (defaults under `/home/olson/bin`) |
| `PRIMERS_FASTA_ROOT` | `fasta_root` | default `/scratch/olson/fasta` |
| `PRIMERS_JOB_STORE` | `check.store` | `redis` (default) or `memory` (tests only; the worker needs redis) |
| `PRIMERS_REDIS_URL` | `check.redis_url` | default `redis://localhost:6380/1` |
| `PRIMERS_GLOBAL_PREFIX` | `check.global_prefix` | host-wide slot keys; default `primers:global:` |
| `PRIMERS_GLOBAL_MAX_JOBS` | `check.global_max_jobs` | host-wide running jobs (default 2); read by the worker |
| `PRIMERS_MAX_QUEUED` | `check.max_queued` | per-site queue limit (default 50) |
| `SWAGGER_HOST`, `SWAGGER_SCHEMES` | (swagger.yaml `host`, `schemes`) | `app.js` only: rewrite the documented host for a dev instance |
| `HOST`, `PORT` | | `app.js` listen address (default all interfaces) and port (default 50011) |

Other keys (`assembly_overrides`, `repeat_masking_overrides`, `check.max_job_cpu_s`, …) are set in `default.yaml`, or
for a one-off run through the config library's `NODE_CONFIG` variable, e.g.
`NODE_CONFIG='{"primers":{"check":{"max_job_cpu_s":3000}}}'`.

- `assembly_overrides: {system_name: prefix}` picks the assembly in a directory with several (`422 AMBIGUOUS_ASSEMBLY`
  lists the candidate prefixes).
- `repeat_masking_overrides: {system_name: soft_masked|unmasked_copy|absent}` overrides the sampled masking state.

### Redis

Redis on port 6380, database 1, shared with gramene-blast. All keys carry a TTL except queues and sorted sets of ids.

| Key | Content |
| --- | --- |
| `primers:<site_key>:job:<id>` | job document (JSON; includes the resolved assemblies while active) |
| `primers:<site_key>:partial:<id>`, `…:result:<id>` | gzipped partial / final results |
| `primers:<site_key>:queue:spec`, `…:queue:pan` | queued job ids |
| `primers:<site_key>:running`, `…:finished` | sorted sets: heartbeat / finish time |
| `primers:<site_key>:worker` | worker lock (30 s, refreshed) |
| `primers:global:slots`, `primers:global:pan_slots` | host-wide running-job slots (`<site_key>:<id>`), pruned after 60 s without heartbeat |

```bash
redis-cli -p 6380 -n 1 --scan --pattern 'primers:sorghum_v11:sorghum11:*' | head
redis-cli -p 6380 -n 1 LRANGE primers:sorghum_v11:sorghum11:queue:pan 0 -1
redis-cli -p 6380 -n 1 ZRANGE primers:global:slots 0 -1 WITHSCORES
redis-cli -p 6380 -n 1 --scan --pattern 'primers:sorghum_v11_dev:*' | xargs -r redis-cli -p 6380 -n 1 DEL   # dev cleanup only
redis-cli -p 6380 INFO memory | grep used_memory_human
```

### Primer3 and BLAST

```bash
cd /home/olson/src && curl -fL -o primer3-2.6.1.tar.gz https://github.com/primer3-org/primer3/archive/refs/tags/v2.6.1.tar.gz
tar xzf primer3-2.6.1.tar.gz && cd primer3-2.6.1/src && make -j4
mkdir -p /home/olson/primer3-2.6.1/bin && cp primer3_core ntthal oligotm /home/olson/primer3-2.6.1/bin/
ln -s /home/olson/primer3-2.6.1/bin/primer3_core /home/olson/bin/primer3_core && primer3_core -about   # 2.6.1
```

The thermodynamic tables are compiled in; do not create `primer3_config/`. BLAST+ 2.13.0 lives in `/home/olson/bin`.
Genome data per `system_name`: `dna/<Prefix>.dna[_sm].toplevel.fa.gz` (bgzip, with `.fai` and `.gzi`),
`<Prefix>.dna.toplevel.*` and `<Prefix>.cdna.all.*` BLAST DBs (`.nal` for multi-volume). The prefix is resolved from the
directory listing (map `_id` suffix, then region matching), never derived from the system name.

### Dev instance

```bash
cd /usr/local/gramene/subsites/sorghum/v11/gramene-swagger-primers
HOST=127.0.0.1 PORT=50111 SWAGGER_HOST=localhost:50111 SWAGGER_SCHEMES=http PRIMERS_SITE_KEY=sorghum_v11_dev PRIMERS_GLOBAL_MAX_JOBS=1 node app.js
PRIMERS_SITE_KEY=sorghum_v11_dev PRIMERS_GLOBAL_MAX_JOBS=1 node api/helpers/primers/jobs/worker_main.js
# workstation: ssh -L 50111:localhost:50111 squam.cshl.edu, then http://localhost:50111/sorghum_v11/docs/?url=/sorghum_v11/swagger
```

Dev jobs use the same Redis database and the same host-wide slots as production (the site key separates the
namespaces), so they count against the host caps. Never bind port 50011 with a dev instance.

### Tests

```bash
npm run test:primers                                                   # node --test "test/primers/unit/**/*.test.js"; offline
PRIMERS_REALDATA=1 node --test --test-concurrency=1 "test/primers/unit/**/*.test.js"   # + mongo, /scratch FASTA, BLAST (serially)
PRIMERS_IT_BASE=http://127.0.0.1:50111/sorghum_v11 node --test --test-concurrency=1 "test/primers/integration/**/*.test.js"
```

- Unit tests never require `app.js`. `contract.test.js` validates `api/swagger/swagger.yaml` with sway, checks the
  request definitions against the handlers' closed sets, and validates every request fixture in
  `test/primers/fixtures/contract/requests/` (from gramene-primers `npm run fixtures`; another directory with
  `PRIMERS_CONTRACT_FIXTURES`). A fixture is either a bare request body (`design*.json` → `POST /primers/design`,
  `check*.json` → `POST /primers/check`) or `{"method", "path", "body", "query", "expect": "valid"|"invalid"}`.
- Integration tests need the dev API **and** the dev worker; `check_jobs.test.js` runs real BLAST jobs one at a time
  and expects check algorithm version 2 (the §10.4 cases as listed under [How a check works](#how-a-check-works)).
  `PRIMERS_IT_EXPECT_JOB_TOO_LARGE=1` enables the cost-guard test, only against an API started with a limit well below
  the ~8,900 CPU-s of its request (already over the default 6,000), e.g. `NODE_CONFIG='{"primers":{"check":{"max_job_cpu_s":1000}}}'`. Against a
  default API that request would be queued.
- `contract.test.js` also runs every fixture through the handlers' own request rules (`check/normalize.js` up to the
  catalog lookup, `design.normalize`), and requires the check fixtures together to send every check param.
- Manual checks (not automated): stop the worker with Ctrl-C during a pan-genome job and start it again — the job ends
  `done` with `attempts: 2`; start the API with `PRIMERS_REDIS_URL=redis://localhost:6399` — `POST /primers/check`
  answers `503 JOB_STORE_UNAVAILABLE` within 3 s while `/primers/design` works; start it with
  `PRIMERS_FASTA_ROOT=/nonexistent` — designs answer JSON errors (`422 NO_SEQUENCE`) and the process stays up.
