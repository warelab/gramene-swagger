# Primer Design API

PCR and qPCR primer design on the assemblies this site serves, Primer-BLAST-style specificity checks against the
reference genome (and its transcriptome), and pan-genome coverage checks across every assembly of the same species.
This documents the `/primers*` endpoints served by this gramene-swagger instance.

- **Base path:** `https://data.sorghumbase.org/sorghum_v11` (local dev: `http://localhost:50111/sorghum_v11`)
- **Controller:** `api/controllers/primers.js`; helpers under `api/helpers/primers/` (check algorithm in `check/`,
  job store and worker in `jobs/`)
- **Engines:** Primer3 2.6.1 (`primer3_core`, spawned per request; `ntthal` for genotyping designs) and BLAST+ 2.13.0
  (`blastn`, `blastdbcmd`)
- **Data:** bgzipped FASTA and BLAST databases under `/scratch/olson/fasta/<system_name>/`; gene models from mongo; known
  variants from Ensembl REST release 115
- **Design spec:** `docs/primer_design_spec.md` (the approved plan's Overrides take precedence over it)
- **Genotyping spec:** `docs/genotyping_design_spec.md`; where it differs, this document describes what is implemented
  ([Genotyping primers](#genotyping-primers-kasp--allele-specific-pcr))

| Endpoint | Purpose | Auth |
| --- | --- | --- |
| `POST /primers/design` | Design primers (synchronous, ≤ 45 s) | none |
| `GET /primers/genomes?system_name=` | Same-species genomes and what each can do | none |
| `POST /primers/check` | Queue a specificity / pan-genome check | none |
| `GET /primers/check/{job_id}` | Poll a check: progress, partial and final results | none |
| `GET /primers/variants?system_name=&region=&start=&end=` | Known variants (Ensembl) in a window, normalized | none |
| `GET /primers/variants/{variant_id}?system_name=` | One Ensembl variation id | none |
| `POST /primers/genotyping/design` | KASP / allele-specific PCR primer sets for one variant (synchronous, ≤ 45 s) | none |

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

curl -s "$BASE/primers/variants?system_name=sorghum_bicolor&region=1&start=11180&end=11290"   # known variants (Ensembl)
post /primers/genotyping/design '{"system_name":"sorghum_bicolor","variant":{"id":"rs871475760","alt":"A"},
  "assay":{"type":"kasp","num_sets":2}}'                                                  # KASP sets + a check request
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
| 400 | validator errors; `INVALID_REQUEST` (`details.field`), `INVALID_PARAMS` (design: `details.param`; check: `details.field`), `INVALID_SEQUENCE`, `TEMPLATE_TOO_LONG`, `REGION_OUT_OF_BOUNDS`, `INTERVAL_OUT_OF_BOUNDS`, `SYSTEM_NAME_MISMATCH`, `PRIMER3_INPUT_ERROR` (`details.primer3_error`), `DUPLICATE_PAIR_ID`, `TOO_MANY_PRIMERS`, `PRODUCT_TOO_LONG_TO_CHECK`, `GENOME_NOT_CHECKABLE`; genotyping: `VARIANT_WINDOW_TOO_LONG`, `INVALID_VARIANT`, `REF_MISMATCH`, `ALT_REQUIRED`, `ALT_NOT_AT_SITE`, `UNSUPPORTED_ALLELE`, `VARIANT_TOO_CLOSE_TO_END`, `VARIANT_TOO_REPETITIVE`, `GENOTYPING_SET_INVALID` |
| 404 | `UNKNOWN_GENE`, `UNKNOWN_TRANSCRIPT`, `UNKNOWN_GENOME`, `UNKNOWN_REGION`, `UNKNOWN_JOB`, `NOT_FOUND`, `UNKNOWN_VARIANT` |
| 405 | `METHOD_NOT_ALLOWED` (`details.allowed_methods`) |
| 413 | body over 100 kb (body-parser shape, below; no `code`) |
| 422 | `NO_SEQUENCE`, `AMBIGUOUS_ASSEMBLY` (`details.candidates`), `NO_BLASTDB`, `JOB_TOO_LARGE`, `NO_VARIATION_DATA`, `AMBIGUOUS_VARIANT_MAPPING`, `VARIANT_NOT_ON_ASSEMBLY` |
| 500 | `PRIMER3_FAILED`, `GENE_STRUCTURE_MISMATCH`, `INTERNAL`, `THERMO_FAILED` |
| 503 | `BUSY` (5 s), `PRIMER3_UNAVAILABLE`, `MONGO_UNAVAILABLE` (30 s), `QUEUE_FULL` (60 s), `JOB_STORE_UNAVAILABLE` (5 s), `FEATURE_DISABLED` (300 s), `VARIATION_SOURCE_UNAVAILABLE` (30 s; 5 s for `queue_full`), `THERMO_UNAVAILABLE` (60 s, or 5 s) |
| 504 | `DEADLINE_EXCEEDED` |

The genotyping codes are described under [Genotyping error codes](#genotyping-error-codes).

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
  "variation": { "available": true, "source": "ensembl", "release": "115" },
  "counts": { "total": 120, "with_blastdb": 120, "with_cdna_blastdb": 120 },
  "genomes": [
    { "system_name": "sorghum_bicolor", "display_name": "Sb bicolor BTx623 v3", "taxon_id": 4558006, "map_id": "GCA_000003195.3",
      "is_query": true, "has_sequence": true, "has_blastdb": true, "has_cdna_blastdb": true, "has_variation": true,
      "repeat_masking": "unmasked_copy", "total_bases": 708735318, "warnings": [] },
    { "system_name": "sorghum_rio", "…": "…", "has_variation": false, "repeat_masking": "soft_masked" } ] }
```

- The query genome is first, then the rest by `display_name`. No filesystem paths are ever returned.
- `variation` (the query genome) and `genomes[].has_variation` say which genomes have known variants from Ensembl; see
  [variation fields](#get-primersgenomes-variation-fields).
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
  `WORKER_LOST`, `MONGO_UNAVAILABLE` (on the job's last attempt: a mongo connection failure during the run, or genes
  queries that the worker found hung or timing out job after job; a genes query timeout alone is not fatal; see
  Operations),
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
  "genotyping": { "…": "…" },   // only for a request with a genotyping block: see "Results: results.genotyping"
  "warnings": [ { "code": "…", "message": "…" } ],
  "timings_ms": { "reference": 11000, "sorghum_353": 1500, "total": 18000 }
}
```

A cDNA group is `{gene_id, orientation, isoforms: [{transcript_id, start, end, size, likelihood, left_mm, right_mm}],
size_min, size_max, likelihood, left_mm, right_mm, left_3p_mm, right_3p_mm, left_mm_pos, right_mm_pos,
terminal_mismatch, approx}`. `unlikely[]` appears on a pair only with `include_unlikely`.

Result warnings: `NO_FASTA_FOR_REALIGN` (mismatch counts are lower bounds, `approx: true`), `ANNOTATION_UNAVAILABLE`
(gene annotation could not be read from mongo from the failure to the end of the job: the reference or genomes whose
annotation finished before it keep their genes, the genome products annotated from then on have `genes: null` and
`ortholog: null`, and cDNA products whose transcripts were not mapped yet are grouped by transcript id without the
isoform suffix. With `details: {cause: "timeout"}` a genes
query timed out twice, so mongo did not answer it within 45 s: it was slow, saturated or hung (see Operations), and
the worker finishes the job with this warning, unless its self-heal valve takes the timeout for a hung mongo and
restarts the worker instead (see Operations). A mongo connection failure makes the worker fail the attempt with
`MONGO_UNAVAILABLE` and retry it after a restart instead, so a finished job carries the warning without `details`
only when run outside the worker),
`TRANSCRIPT_GENE_UNMAPPED`, `PANGENOME_TRANSCRIPT_MODELS_ONLY`,
and assembly warnings; genotyping jobs can add the
[genotyping results warnings](#results-warnings). Submit warnings: `EXPECTED_IGNORED`, `MAX_PRODUCT_SIZE_RAISED`, `GENOMES_IGNORED`, and the
reference assembly's `ASSEMBLY_MISMATCH` / `AMBIGUOUS_ASSEMBLY` prefixed with its system name.

---

## Genotyping primers (KASP / allele-specific PCR)

Design primers that tell the two alleles of one variant apart, then check them across the pan-genome. A **set** is two
allele-specific (AS) primers, whose 3′ base is the base that differs between the REF and the ALT allele, plus one
**common** primer on the other side of the variant. KASP sets carry the FAM and HEX 5′ tails of the KASP chemistry; sets
for gel-based allele-specific PCR (`as_pcr`) carry a deliberate mismatch near the 3′ end instead.

- **Variants** come from Ensembl REST release 115, always through this API (a browser never calls Ensembl), or are
  entered by hand on any genome that has sequence. On this site only `sorghum_bicolor` has Ensembl variation data.
- **Identity:** a variant is its `key`, `region:position:REF:ALT` in left-aligned VCF form (`1:11109:C:A`, `1:11282:CA:C`,
  `1:11502:C:CGT`). Designs and check jobs use the key; Ensembl ids are for lookup and display.
- **Design spec:** `docs/genotyping_design_spec.md`; the § numbers below refer to it. Where the spec and this section
  differ, this section describes what is implemented, and the spec's last section lists the known differences.

| Endpoint | Purpose |
| --- | --- |
| `GET /primers/genomes?system_name=` | Gains `variation` and `genomes[].has_variation`: which genomes have Ensembl variants |
| `GET /primers/variants?system_name=&region=&start=&end=` | Known variants in a window, normalized, merged and checked against the genome |
| `GET /primers/variants/{variant_id}?system_name=` | One Ensembl variation id |
| `POST /primers/genotyping/design` | KASP or AS-PCR sets for one variant (synchronous, ≤ 45 s) |
| `POST /primers/check` | Gains an optional `genotyping` block: which allele each assembly carries |

**Typical flow**

1. `GET /primers/genomes`: offer the Ensembl variant picker when `variation.available` is true; otherwise ask for region,
   position, ref and alt.
2. `GET /primers/variants` over a window (the gene ± 2 kb, say), or `GET /primers/variants/{variant_id}` for a deep link.
3. `POST /primers/genotyping/design` with the chosen id (or the entry's `vcf`), or with a manual variant.
4. Export the order sheet from `sets[].order`, `assay.kasp_mix` and `variant.submission_sequence`.
5. `POST /primers/check` with the design's `check.request` (or the sets you choose), then poll `GET /primers/check/{job_id}`
   for `results.genotyping`.

### Variants

Every variant object (in list, lookup and design responses) carries these representations. The values are real:

| Field | Meaning | rs871475760 | rs5413864115 | tmp_1_11502_C_CGT |
| --- | --- | --- | --- | --- |
| `key` | `region:position:REF:ALT`, left-aligned VCF with an anchor base for indels | `1:11109:C:A` | `1:11282:CA:C` | `1:11502:C:CGT` |
| `vcf` | `{position, ref, alt}` of the key; what the check takes | `{11109, C, A}` | `{11282, CA, C}` | `{11502, C, CGT}` |
| `minimal` | Ensembl style `{start, end, ref, alt}`, `-` for an empty allele; an insertion has `start = end + 1` | `{11109, 11109, C, A}` | `{11283, 11283, A, -}` | `{11503, 11502, -, GT}` |
| `label` | display text | `1:11109 C/A` | `1:11283-11283 A/-` | `1:11502^11503 -/GT` |
| `kind` | `snv`, `mnv`, `insertion`, `deletion` or `complex` | `snv` | `deletion` | `insertion` |
| `shift` | how many bases the event can slide right and still give the same haplotype | 0 | 2 | 0 |
| `zone` | `{start, end}`: the VCF span, both discriminating positions and both ALT 3′ anchor mappings. Only an allele-specific primer may touch it | `{11109, 11109}` | `{11282, 11286}` | `{11502, 11503}` |
| `discriminating.forward` | `{position, ref_base, alt_base, alt_maps_to}`: the 3′ base of a forward AS primer, the first base (left to right) where the haplotypes differ. `alt_maps_to` is the reference coordinate of the ALT base, `null` for an inserted base | `{11109, C, A, 11109}` | `{11285, A, G, 11286}` | `{11503, A, G, null}` |
| `discriminating.reverse` | the same, scanning right to left (plus-strand bases) | `{11109, C, A, 11109}` | `{11283, A, C, 11282}` | `{11502, C, T, null}` |

- `ids`: every Ensembl id of the event, the requested id first, then `rs` ids by number, then the rest alphabetically.
  Records describing the same event (an EVA and a SAP record of one insertion, say) merge into one entry.
- `synonyms`: Ensembl's synonyms of the requested id (lookup and design by id; `.` is dropped); `[]` in listings.
- `alleles`, `multiallelic`: the site's alleles, reference first; `multiallelic` is `{alleles, other_alts}` or `null`.
  A multi-allelic record gives one entry per alternative allele.
- `records[]`: `{id, source, ems}`. `source` is Ensembl's short source code (`EVA`, `EMS_PMID38100514_Jiao`,
  `SAP_PMID35653240_Boatwri`), or `null` for an id known only from its lookup. `ems` is true when the source matches
  `variation.ems_source_pattern` (`^EMS_`). A record is never dropped because of its id; an id Ensembl reports in an
  unusable form becomes `null`.
- `ems`: every merged record is EMS. EMS mutations are private to BTx623-background mutant lines.
- `consequence`, `ref_verified` (REF matches this site's FASTA), `designable` and `issues[]` (`{code, message, details}`,
  codes `REF_MISMATCH`, `STAR_ALLELE`, `ALLELE_TOO_LONG`, `UNSUPPORTED_ALLELE`, `REPEAT_TOO_LONG`). Only `STAR_ALLELE` (a
  `*` sibling allele) leaves an entry designable. With `REPEAT_TOO_LONG` (the event slides more than `variation.max_shift`,
  1,000 bp, or up to the region end) `zone` and `discriminating` are `null`.
- Design responses add `requested_id` (the id, or `null` for a manual variant) right after `key`, and
  `submission_sequence`.

**Manual input** (`variant` in a design request): alleles are `A`, `C`, `G`, `T` in any case, at most 50 nt, or `-` for
an empty allele.

| Style | When | `position` is | rs5413864115 entered this way |
| --- | --- | --- | --- |
| VCF | neither allele is `-` | the first base of `ref` | `{"region": "1", "position": 11282, "ref": "CA", "alt": "C"}` |
| Ensembl | exactly one allele is `-` | deletion: the first deleted base; insertion: the base after the insertion point (`minimal.start`) | `{"region": "1", "position": 11283, "ref": "A", "alt": "-"}` |

The server left-aligns every input, so both rows give key `1:11282:CA:C`. The check's `genotyping.variant` takes VCF style
only.

**Orientation.** Coordinates are 1-based inclusive on the reference plus strand, and the genotyping template is always
the plus strand.

| `orientation` | allele-specific primer | its 3′ end | common primer |
| --- | --- | --- | --- |
| `forward` | left primer, plus strand | `discriminating.forward.position` | right primer, after `zone.end` |
| `reverse` | right primer, minus strand | `discriminating.reverse.position` | left primer, before `zone.start` |

### Ensembl, and what happens without it

- **Calls.** `GET {variation.base_url}/overlap/region/{species}/{region}:{start}-{end}?feature=variation;content-type=application/json`
  in fixed 10 kb chunks (`variation.chunk_bp`), and `GET {variation.base_url}/variation/{species}/{id}?content-type=application/json`.
  - Only the species in `variation.species` are ever called (`sorghum_bicolor` here).
  - The region must be a sequence of the assembly and the id must match the id pattern before either goes into a URL,
    and both are URL-encoded.
  - Redirects are refused. Each call has an 8 s timeout and a body cap: 1 MB per chunk, 256 KB per lookup.
- **Resilience,** per API process:
  - At most 4 outbound calls run at a time (`variation.max_concurrent`). A call that cannot start within 5 s, or finds 256
    calls already waiting, fails with reason `queue_full` (`retry_after_s` 5).
  - There is one in-flight request per chunk or id, shared by every caller. A client that disconnects stops waiting
    without cancelling it.
  - A circuit breaker opens after 3 failures within 60 s and stays open for 60 s (reason `breaker_open`, `retry_after_s` =
    the seconds left). A success closes it.
  - An LRU cache of 500 entries holds validated records only. Results are kept for 1 h (release 115 is static), unknown
    ids for 5 min, and failures (`timeout`, `transport`, `http_5xx`, `rate_limited`) for 30 s.
  - A request repeated inside those 30 s answers 503 at once, with the seconds left, and does not count again. The breaker
    therefore opens only after failures on 3 different chunks or ids.

| Ensembl answers | Result |
| --- | --- |
| 200 with JSON of the expected shape | the records, cached |
| 400 on `/variation/…` with a JSON `error` matching `not found` | `404 UNKNOWN_VARIANT`; cached 5 min, not a breaker failure |
| another 4xx, a redirect, a non-JSON body (e.g. an HTML 404 from a wrong base URL), a body over the cap, an unexpected shape | `503 VARIATION_SOURCE_UNAVAILABLE`, reason `invalid_response`, `retry_after_s` 30; a breaker failure, not cached |
| 429, 5xx, a timeout, a network error | reason `rate_limited`, `http_5xx`, `timeout` or `transport`, `retry_after_s` 30; cached 30 s; a breaker failure |

Records are kept when they are variation features of the requested region with integer coordinates and 2–10 alleles made
of `ACGTN` bases, `-` or `*`. Other records are dropped and counted in warning `VARIATION_RECORDS_SKIPPED {count, reasons}`.

**When Ensembl is unavailable:**

| Call | Answer |
| --- | --- |
| `GET /primers/variants`, `GET /primers/variants/{variant_id}` | `503 VARIATION_SOURCE_UNAVAILABLE {retry_after_s, reason}`; offer manual entry |
| `POST /primers/genotyping/design` with `variant.id` | the same 503 |
| `POST /primers/genotyping/design` with a manual variant, on a genome with variation data | `200`: the sets are designed, with `neighbours.data: "unavailable"` and warning `NEIGHBOURS_UNAVAILABLE {reason}`. No known variant was screened, and `variant.ids` stays empty |
| `POST /primers/genotyping/design` on a genome without variation data | `200` with `neighbours.data: "none"` and warning `NO_VARIATION_DATA`; Ensembl is never called |
| `POST /primers/check`, the check worker, `POST /primers/design`, `GET /primers/genomes` | never call Ensembl |

All the Ensembl work of a design (resolving the id, fetching the neighbours) happens **before** the design takes a Primer3
slot, inside its 45 s deadline. A slow or failing Ensembl therefore never holds one of the 4 slots, and never turns an
ordinary `POST /primers/design` into `503 BUSY`.

### Feature flags

| Setting | Effect |
| --- | --- |
| `PRIMERS_ENABLED=0` (`enabled`) | Every `/primers` endpoint, the new ones included, answers `503 FEATURE_DISABLED` (retry after 300 s). |
| `PRIMERS_VARIATION_ENABLED=0` (`variation.enabled`) | Ensembl is switched off. `GET /primers/variants` and `GET /primers/variants/{variant_id}` answer `503 FEATURE_DISABLED` (300 s), and so does a design with `variant.id` on a genome listed in `variation.species` (on other genomes that is `422 NO_VARIATION_DATA`). Manual designs keep working on every genome, as on a genome without data: `neighbours.data: "none"` and warning `NO_VARIATION_DATA`, whose message says that known-variant lookups are disabled on this server (not that the genome lacks data). `GET /primers/genomes` reports `variation.available: false` and `has_variation: false` everywhere. Checks are unaffected. |
| `PRIMERS_VARIATION_URL` (`variation.base_url`) | Another Ensembl REST base URL. Accepted only as `https://…`, or as `http://127.0.0.1[:port]` for a loopback test server, without credentials, query or fragment. Any other value is ignored with a startup warning. |

### `GET /primers/genomes`: variation fields

Two additive fields; nothing else changed.

- `variation` describes the query genome as `{available, source, release}`. It is never `null`: a genome without data has
  `{available: false, source: null, release: null}`.
- `genomes[].has_variation` is true when `variation.enabled` is on, the genome's `system_name` is a key of
  `variation.species`, and its assembly resolves with a FASTA. On sorghum_v11 that is `sorghum_bicolor` alone (1 of 120).
  No Ensembl call is made.

<!-- example: request GET /primers/genomes capture=capture-genomes-sorghum_bicolor.json#/request/path -->
```http
GET /sorghum_v11/primers/genomes?system_name=sorghum_bicolor
```

<!-- example: response GET /primers/genomes 200 capture=capture-genomes-sorghum_bicolor.json#/response -->
```json
{
  "system_name": "sorghum_bicolor", "species": { "taxon_id": 4558, "name": "Sorghum bicolor" },
  "variation": { "available": true, "source": "ensembl", "release": "115" },
  "counts": { "total": 120, "with_blastdb": 120, "with_cdna_blastdb": 120 },
  "genomes": [
    {
      "system_name": "sorghum_bicolor", "display_name": "Sb bicolor BTx623 v3", "taxon_id": 4558006,
      "map_id": "GCA_000003195.3", "is_query": true, "has_sequence": true, "has_blastdb": true,
      "has_cdna_blastdb": true, "has_variation": true, "repeat_masking": "unmasked_copy", "total_bases": 708735318,
      "warnings": []
    },
    {
      "system_name": "sorghum_rio", "display_name": "Sb bicolor PI651496 Rio", "taxon_id": 4558116,
      "map_id": "GCA_015952705.1", "is_query": false, "has_sequence": true, "has_blastdb": true,
      "has_cdna_blastdb": true, "has_variation": false, "repeat_masking": "soft_masked", "total_bases": 729379862,
      "warnings": []
    },
    "…"
  ]
}
```

### `GET /primers/variants`

The known variants of a window on a genome with variation data, one entry per `key`.

| param | in | type | required | rule |
| --- | --- | --- | --- | --- |
| `system_name` | query | `^[a-z0-9_]+$`, ≤ 128 | yes | a genome with variation data, else `422 NO_VARIATION_DATA` |
| `region` | query | string ≤ 255 | yes | a sequence of the assembly, else `404 UNKNOWN_REGION` |
| `start`, `end` | query | integer ≥ 1 | yes | `start ≤ end ≤` region length, else `400 REGION_OUT_OF_BOUNDS`; at most 50,000 bp (`variation.max_window`), else `400 VARIANT_WINDOW_TOO_LONG {length, max}` |
| `types` | query | comma-separated `snv`, `mnv`, `insertion`, `deletion`, `complex` | no | default all |
| `include_ems` | query | boolean | no | default `true`; `false` drops entries whose records are all EMS |
| `limit` | query | integer 1–5000 | no | default 2000 |

- An entry is in the window when its `minimal` span overlaps `[start, end]`. Entries are sorted by `vcf.position`, then
  `key`.
- A multi-allelic record gives one entry per alternative allele, and records of the same event merge into one entry
  whatever their ids. `synonyms` is always `[]` here.
- Entries that cannot be designed are listed too, with `designable: false` and their `issues`, so a picker can say why.
- `total` counts the entries that pass `types` and `include_ems`, and `returned` is at most `limit`. When more matched,
  `truncated` is `true` and warning `VARIANTS_TRUNCATED {returned, total, limit}` is added.
- REF is checked against this site's FASTA (`ref_verified`); disagreeing entries add warning `REF_MISMATCHES {count}`.
- The window is read from the FASTA once (± 1,250 bp). The 10 kb Ensembl chunks are cached for an hour and shared with
  later designs and lookups at the same locus. The example below took 0.32 s with cold caches.

<!-- example: request GET /primers/variants capture=capture-variants-list-1_11180-11290.json#/request/path -->
```http
GET /sorghum_v11/primers/variants?system_name=sorghum_bicolor&region=1&start=11180&end=11290
```

<!-- example: response GET /primers/variants 200 capture=capture-variants-list-1_11180-11290.json#/response -->
```json
{
  "system_name": "sorghum_bicolor", "region": "1", "start": 11180, "end": 11290,
  "source": { "name": "ensembl", "release": "115" }, "total": 4, "returned": 4, "truncated": false,
  "variants": [
    {
      "key": "1:11182:A:G", "ids": ["rs873026643"], "synonyms": [], "label": "1:11182 A/G", "kind": "snv",
      "region": "1", "vcf": { "position": 11182, "ref": "A", "alt": "G" },
      "minimal": { "start": 11182, "end": 11182, "ref": "A", "alt": "G" }, "alleles": ["A", "G"],
      "multiallelic": null, "shift": 0, "zone": { "start": 11182, "end": 11182 },
      "discriminating": {
        "forward": { "position": 11182, "ref_base": "A", "alt_base": "G", "alt_maps_to": 11182 },
        "reverse": { "position": 11182, "ref_base": "A", "alt_base": "G", "alt_maps_to": 11182 }
      },
      "records": [{ "id": "rs873026643", "source": "EVA", "ems": false }], "ems": false,
      "consequence": "3_prime_UTR_variant", "ref_verified": true, "designable": true, "issues": []
    },
    {
      "key": "1:11193:C:T", "ids": ["tmp_1_11193_C_T"], "synonyms": [], "label": "1:11193 C/T", "kind": "snv",
      "region": "1", "vcf": { "position": 11193, "ref": "C", "alt": "T" },
      "minimal": { "start": 11193, "end": 11193, "ref": "C", "alt": "T" }, "alleles": ["C", "T"],
      "multiallelic": null, "shift": 0, "zone": { "start": 11193, "end": 11193 },
      "discriminating": {
        "forward": { "position": 11193, "ref_base": "C", "alt_base": "T", "alt_maps_to": 11193 },
        "reverse": { "position": 11193, "ref_base": "C", "alt_base": "T", "alt_maps_to": 11193 }
      },
      "records": [{ "id": "tmp_1_11193_C_T", "source": "EMS_PMID38100514_Jiao", "ems": true }], "ems": true,
      "consequence": "3_prime_UTR_variant", "ref_verified": true, "designable": true, "issues": []
    },
    {
      "key": "1:11203:C:T", "ids": ["tmp_1_11203_C_T"], "synonyms": [], "label": "1:11203 C/T", "kind": "snv",
      "region": "1", "vcf": { "position": 11203, "ref": "C", "alt": "T" },
      "minimal": { "start": 11203, "end": 11203, "ref": "C", "alt": "T" }, "alleles": ["C", "T"],
      "multiallelic": null, "shift": 0, "zone": { "start": 11203, "end": 11203 },
      "discriminating": {
        "forward": { "position": 11203, "ref_base": "C", "alt_base": "T", "alt_maps_to": 11203 },
        "reverse": { "position": 11203, "ref_base": "C", "alt_base": "T", "alt_maps_to": 11203 }
      },
      "records": [{ "id": "tmp_1_11203_C_T", "source": "EMS_PMID29378822_Addo-Qu", "ems": true }], "ems": true,
      "consequence": "3_prime_UTR_variant", "ref_verified": true, "designable": true, "issues": []
    },
    {
      "key": "1:11282:CA:C", "ids": ["rs5413864115"], "synonyms": [], "label": "1:11283-11283 A/-",
      "kind": "deletion", "region": "1", "vcf": { "position": 11282, "ref": "CA", "alt": "C" },
      "minimal": { "start": 11283, "end": 11283, "ref": "A", "alt": "-" }, "alleles": ["A", "-"],
      "multiallelic": null, "shift": 2, "zone": { "start": 11282, "end": 11286 },
      "discriminating": {
        "forward": { "position": 11285, "ref_base": "A", "alt_base": "G", "alt_maps_to": 11286 },
        "reverse": { "position": 11283, "ref_base": "A", "alt_base": "C", "alt_maps_to": 11282 }
      },
      "records": [{ "id": "rs5413864115", "source": "EVA", "ems": false }], "ems": false,
      "consequence": "3_prime_UTR_variant", "ref_verified": true, "designable": true, "issues": []
    }
  ],
  "warnings": []
}
```

A genome without variation data:

<!-- example: request GET /primers/variants capture=capture-variants-list-no-variation-data.json#/request/path -->
```http
GET /sorghum_v11/primers/variants?system_name=sorghum_rio&region=1&start=1&end=100
```

<!-- example: response GET /primers/variants 422 capture=capture-variants-list-no-variation-data.json#/response -->
```json
{
  "message": "sorghum_rio has no known-variant data; enter the variant as region, position, ref and alt",
  "code": "NO_VARIATION_DATA", "details": { "system_name": "sorghum_rio" }
}
```

**Errors:** validator `400`; `400 VARIANT_WINDOW_TOO_LONG`, `REGION_OUT_OF_BOUNDS`; `404 UNKNOWN_GENOME`, `UNKNOWN_REGION`;
`422 NO_VARIATION_DATA`, `NO_SEQUENCE`, `AMBIGUOUS_ASSEMBLY`; `503 VARIATION_SOURCE_UNAVAILABLE`, `FEATURE_DISABLED`,
`MONGO_UNAVAILABLE`.
**Warnings:** `VARIATION_RECORDS_SKIPPED {count, reasons}`, `REF_MISMATCHES {count}`, `VARIANTS_TRUNCATED {returned, total, limit}`.

### `GET /primers/variants/{variant_id}`

Resolves one Ensembl variation id.

| param | in | type | required | notes |
| --- | --- | --- | --- | --- |
| `variant_id` | path | `^[A-Za-z0-9][A-Za-z0-9_.:,*-]{0,254}$` | yes | URL-encode it: real ids contain `,` and `*` (`tmp_1_13549_TTA_T%2C%2A`) and run to 255 characters |
| `system_name` | query | `^[a-z0-9_]+$`, ≤ 128 | yes | a genome with variation data |

- Returns every entry of the id, one per alternative allele, including the non-designable ones with their `issues`.
- Each entry carries every Ensembl id of the same event; a merge adds warning `DUPLICATE_VARIANT_IDS {key, ids}`.
- `requested_id` is a top-level member; the entries do not repeat it.
- `synonyms` holds Ensembl's synonyms of the id: `GET /primers/variants/rs871475760?system_name=sorghum_bicolor` gives key
  `1:11109:C:A` with `synonyms: ["tmp_1_11109_C_A"]`.
- `records[].source` comes from the overlap records at the id's position, because the lookup itself reports a long title.
  An id that the overlap does not list is reported from its lookup mapping, with `source: null`.
- The id must map to exactly one sequence of this assembly: no mapping is `422 VARIANT_NOT_ON_ASSEMBLY {id}`, several are
  `422 AMBIGUOUS_VARIANT_MAPPING {id, mappings}`.

<!-- example: request GET /primers/variants/{variant_id} capture=capture-variants-lookup-tmp_1_11502_C_CGT.json#/request/path -->
```http
GET /sorghum_v11/primers/variants/tmp_1_11502_C_CGT?system_name=sorghum_bicolor
```

<!-- example: response GET /primers/variants/{variant_id} 200 capture=capture-variants-lookup-tmp_1_11502_C_CGT.json#/response -->
```json
{
  "requested_id": "tmp_1_11502_C_CGT", "system_name": "sorghum_bicolor",
  "source": { "name": "ensembl", "release": "115" },
  "variants": [
    {
      "key": "1:11502:C:CGT", "ids": ["tmp_1_11502_C_CGT", "rs5413863549"], "synonyms": [],
      "label": "1:11502^11503 -/GT", "kind": "insertion", "region": "1",
      "vcf": { "position": 11502, "ref": "C", "alt": "CGT" },
      "minimal": { "start": 11503, "end": 11502, "ref": "-", "alt": "GT" }, "alleles": ["-", "GT"],
      "multiallelic": null, "shift": 0, "zone": { "start": 11502, "end": 11503 },
      "discriminating": {
        "forward": { "position": 11503, "ref_base": "A", "alt_base": "G", "alt_maps_to": null },
        "reverse": { "position": 11502, "ref_base": "C", "alt_base": "T", "alt_maps_to": null }
      },
      "records": [
        { "id": "tmp_1_11502_C_CGT", "source": "SAP_PMID35653240_Boatwri", "ems": false },
        { "id": "rs5413863549", "source": "EVA", "ems": false }
      ],
      "ems": false, "consequence": "3_prime_UTR_variant", "ref_verified": true, "designable": true, "issues": []
    }
  ],
  "warnings": [
    {
      "code": "DUPLICATE_VARIANT_IDS",
      "message": "2 Ensembl ids describe the same event 1:11502:C:CGT; they were merged",
      "details": { "key": "1:11502:C:CGT", "ids": ["tmp_1_11502_C_CGT", "rs5413863549"] }
    }
  ]
}
```

An unknown id (Ensembl answers 400 "not found", which is cached for 5 minutes):

<!-- example: request GET /primers/variants/{variant_id} capture=capture-variants-lookup-unknown-variant.json#/request/path -->
```http
GET /sorghum_v11/primers/variants/rs0000000001?system_name=sorghum_bicolor
```

<!-- example: response GET /primers/variants/{variant_id} 404 capture=capture-variants-lookup-unknown-variant.json#/response -->
```json
{
  "message": "unknown variant rs0000000001 for sorghum_bicolor", "code": "UNKNOWN_VARIANT",
  "details": { "id": "rs0000000001", "system_name": "sorghum_bicolor" }
}
```

**Errors:** validator `400` (the id pattern); `404 UNKNOWN_VARIANT {id, system_name}`, `UNKNOWN_GENOME`;
`422 NO_VARIATION_DATA`, `NO_SEQUENCE`, `AMBIGUOUS_ASSEMBLY`, `AMBIGUOUS_VARIANT_MAPPING`, `VARIANT_NOT_ON_ASSEMBLY`;
`503 VARIATION_SOURCE_UNAVAILABLE`, `FEATURE_DISABLED`, `MONGO_UNAVAILABLE`.
**Warnings:** `DUPLICATE_VARIANT_IDS {key, ids}`, `VARIATION_RECORDS_SKIPPED`, `REF_MISMATCHES`.

### `POST /primers/genotyping/design`

Synchronous. Designs up to `num_sets` sets for one variant, in both orientations, and returns them ranked with order rows
and a ready check request.

#### Request

| field | type | notes |
| --- | --- | --- |
| `system_name` | `^[a-z0-9_]+$`, ≤ 128 | required |
| `variant` | object | required: exactly one of `{id [, alt]}` or `{region, position, ref, alt}` |
| `variant.id` | the Ensembl id pattern | the genome needs variation data; `alt` picks the allele at a multi-allelic site, where it is required |
| `variant.region`, `variant.position` | string ≤ 255, integer ≥ 1 | manual input |
| `variant.ref`, `variant.alt` | `^([ACGTacgt]{1,50}\|-)$` | manual input, in [VCF or Ensembl style](#variants) |
| `assay` | object | below; unset values take the defaults of `assay.type` |
| `params` | object | a closed subset of the design params, below |
| `avoid_repeats` | bool | default `false`. When true the template is repeat-masked as in `/primers/design`, except the allele-specific primer window, which is always exempt (`template.features.exempt`) |
| `repeat_mask_mode` | `n_mask` (default) \| `three_prime` | only with `avoid_repeats` |
| `template_only` | bool | resolve and verify the variant only. The response has `variant`, `template`, `assay`, `neighbours`, `settings`, `engine` and `warnings`, with `orientations: null`, `sets: []` and `check: null`. Primer3 is not run, and no Primer3 slot is taken unless `avoid_repeats` is set |
| `label` | `^[A-Za-z0-9_.-]{1,40}$` | prefix of the order-sheet oligo names; default: the first variant id, else the key (other characters become `_`) |

**Assay.** Any combination is accepted (tails on `as_pcr`, a deliberate mismatch on `kasp`); `assay` in the response echoes
the effective values.

| field | values | `kasp` default | `as_pcr` default |
| --- | --- | --- | --- |
| `type` | `kasp` \| `as_pcr` | `kasp` when omitted | |
| `orientation` | `both` \| `forward` \| `reverse` | `both` | `both` |
| `tails` | `none` \| `ref_fam_alt_hex` \| `ref_hex_alt_fam` | `ref_fam_alt_hex` | `none` |
| `deliberate_mismatch` | `none` \| `auto` | `none` | `auto` |
| `mismatch_position` | `2` \| `3`, the distance from the 3′ end (Little 1995 Table 9.8.1 is defined for 2; 3 is extrapolated) | `2` | `2` |
| `num_sets` | 1–10 (`genotyping.max_sets`) | 6 (`genotyping.num_sets_default`) | 6 |
| `max_relaxation` | 0–2 | 2 | 2 |
| `neighbour_policy` | `avoid_3p` \| `ignore` | `avoid_3p` | `avoid_3p` |

The KASP tails are FAM `GAAGGTGACCAAGTTCATGCT` and HEX `GAAGGTCGGAGTCAACGGATT`.

**Params, presets and the relaxation ladder.**

- **Accepted params:** `opt_size`, `min_size`, `max_size`, `opt_tm`, `min_tm`, `max_tm`, `opt_gc`, `min_gc`, `max_gc`,
  `max_tm_diff`, `max_poly_x`, `gc_clamp`, `max_end_stability`, `salt_monovalent`, `salt_divalent`, `dntp_conc` and
  `dna_conc`, with the bounds of [`POST /primers/design`](#params-presets-and-primer3-tags). `product_size_ranges` takes
  1–4 ranges of integers 20–1000. Anything else (`num_return`, `max_ns`, the junction params) is rejected.
- **Defaults and levels:** unset params come from the preset. Each ladder level applies its changes on top of the previous
  level.
- **Pinned params:** a param the client set is never changed by any level (`settings.pinned`).
- **Product minimum:** at every level the effective product-range minimum is raised to `2 × max_size + 1`.

| param | `kasp` L0 | L1 | L2 | `as_pcr` L0 | L1 | L2 |
| --- | --- | --- | --- | --- | --- | --- |
| `opt_size` / `min_size` / `max_size` | 22 / 18 / 30 | max 32 | – | 24 / 18 / 30 | max 32 | – |
| `opt_tm` / `min_tm` / `max_tm` | 60 / 57 / 63 | min 55, max 65 | min 52 | 60 / 57 / 63 | min 55, max 65 | min 52 |
| `min_gc` / `max_gc` | 30 / 70 | 20 / 80 | – | 30 / 70 | 20 / 80 | – |
| `max_tm_diff` | 3 | – | 6 | 3 | – | 6 |
| `max_poly_x` | 5 | – | – | 4 | – | – |
| `product_size_ranges` | `[[50,120]]`, effective `[[61,120]]` | `[[50,150]]`, effective `[[65,150]]` | – | `[[150,300]]` | – | – |

- **Hard floors, never relaxed.** Every allele-specific primer needs a Tm of at least 52 °C (`genotyping.as_min_tm`; for a
  deliberate-mismatch primer, the Tm of its perfect-match sequence). The derived primers (the ALT primer and any mismatch
  primer) also need GC ≥ 15 % (`genotyping.as_min_gc`). A set that fails is dropped (`rejected.below_floor`).
  Primer3's `PRIMER_PICK_ANYWAY` is never used to design.
- **Where the ladder stops.** It stops at the first level that keeps a set. Sets from level 1 or 2 carry that
  `relaxation_level`, and the orientation adds warning `RELAXED_CONSTRAINTS`. `settings.ladder` lists each level's changes,
  and `settings.floors` the floors.

#### How a design works

1. The request rules (below) are checked, and the 45 s deadline starts (`design.deadline_ms`; the semaphore wait counts).
2. The design resolves:
   - the assembly;
   - the variant: by id through Ensembl, or manual;
   - REF against the FASTA, the shift, both discriminating positions and the zone;
   - the known variants of the template window.

   **All Ensembl work happens here, before step 3.**
3. The design takes a Primer3 slot of the design semaphore. The slots are shared with `POST /primers/design`: 4 running and 16
   waiting per API process, a 10 s wait, then `503 BUSY`.
4. The template is built: the plus-strand reference from the lower discriminating position − F to the higher one + F,
   with `F = max(400, largest product-range bound + 40)` (`genotyping.template_flank`), plus its ALT haplotype. An optional
   repeat mask exempts the allele-specific window.
5. For each orientation:
   - Natural targets with `neighbour_policy: "avoid_3p"` get the neighbour check first: a known non-EMS variant in the
     last 5 nt of the allele-specific primer blocks the orientation.
   - Primer3 runs once per ladder level, with the allele-specific 3′ end forced (`SEQUENCE_FORCE_LEFT_END` or
     `SEQUENCE_FORCE_RIGHT_END`) and a 10 nt `SEQUENCE_TARGET` guard just outside the zone. The guard keeps the common
     primer off the variant on both haplotypes.
6. Each returned pair, in Primer3 order, is rejected when:
   - Primer3 ignored the forced end (`rejected.force`);
   - the common primer touches the zone or the allele-specific primer on either haplotype (`rejected.overlap`);
   - the common primer ends on a known non-EMS variant where neighbours block (`rejected.common_neighbour_3p`);
   - it repeats an earlier pair (`rejected.duplicate`).
7. Scoring, for at most 8 candidates per orientation:
   - The ALT primer is derived from the REF primer: same 5′ end, 3′ end at the first haplotype difference, which also
     covers indels.
   - A deliberate mismatch, when asked for, is looked up in Little 1995 Table 9.8.1.
   - Derived oligos are scored by Primer3 `check_primers` runs on their own haplotype (`rejected.alt_scoring_failed` when
     that fails), and the floors apply.
   - Tailed hairpins and dimers, and the duplex Tm of mismatch primers, come from `ntthal` with the same salts.
   - Discrimination is the check's own classification of each AS primer on the other allele, not a ΔTm.
8. Issues, quality, score and ranking; order rows; the check request. The slot is released.

#### Process budget

Each design request has a fixed budget:

| Limit | Value | Config |
| --- | --- | --- |
| Pairs returned by each Primer3 design run | 20 | `genotyping.num_return_per_run` |
| Candidates scored per orientation, over all levels | 8 | `genotyping.max_scored_per_orientation` |
| Primer3 runs per request (design and `check_primers` runs, each counted when it starts) | 54 | `genotyping.max_primer3_runs` |
| Distinct `ntthal` calls per request | 272 | `genotyping.max_thermo_calls` |
| `ntthal` calls at a time, and per-call timeout | 4 per request, 5 s | `genotyping.thermo_concurrency`, `genotyping.thermo_timeout_ms` |

- **Reservation.** Before scoring a candidate, the design reserves its worst case: 1 `check_primers` run (3 with a
  deliberate mismatch) and 15 `ntthal` calls with tails (plus 2 with a deliberate mismatch).
- **When a reservation does not fit:**
  - That pair and the rest of its level stay unscored (`attempts[].not_scored`), and the orientation's ladder stops.
  - **The response is still a `200`** with the sets scored so far, plus warning `DESIGN_BUDGET_EXHAUSTED {orientation,
    primer3_runs, thermo_calls, max_primer3_runs, max_thermo_calls, not_scored, sets_returned}`. The counters are as at
    the end of the request.
  - An orientation left without a set reports `status: "no_sets"` with `reason: "budget_exhausted"`.
  - It is never a `500`.
- **Timing.** The examples below stay within the default caps; they took 0.44–1.0 s over HTTP on the dev instance with
  live Ensembl. The largest case the caps allow, 54 Primer3 runs and 272 `ntthal` calls, measured about 1.15 s (spec §4.16).

#### Handler rules and errors

In this order; everything up to the Primer3 slot runs before a slot is taken.

| Rule | Error |
| --- | --- |
| Unknown keys (swagger rejects them first) | `400 INVALID_REQUEST {field}` |
| `variant` has `id` together with `region`, `position` or `ref`, or has neither `id` nor all four manual fields | `400 INVALID_VARIANT {reason: "id_or_manual"}` |
| Manual alleles both `-`, equal, or `-` against something other than bases | `400 INVALID_VARIANT {reason: "alleles"}` |
| An allele longer than 50 nt | `400 INVALID_VARIANT {reason: "allele_too_long"}` |
| The params' cross-field rules (`min ≤ opt ≤ max`, `gc_clamp ≤ min_size`, each range `a < b`), at every ladder level | `400 INVALID_PARAMS {param}`, plus `level` when a relaxed level breaks them |
| The genome and its assembly | `404 UNKNOWN_GENOME`; `422 AMBIGUOUS_ASSEMBLY`, `NO_SEQUENCE`; `503 MONGO_UNAVAILABLE` |
| By id: the genome has variation data | `422 NO_VARIATION_DATA {system_name}`; `503 FEATURE_DISABLED` when variation is switched off |
| By id: Ensembl resolves the id onto this assembly | `404 UNKNOWN_VARIANT {id, system_name}`; `422 VARIANT_NOT_ON_ASSEMBLY`, `AMBIGUOUS_VARIANT_MAPPING`; `503 VARIATION_SOURCE_UNAVAILABLE` |
| By id: `alt` omitted where the site has several designable alternative alleles | `400 ALT_REQUIRED {id, alts}` |
| By id: `alt` is not an allele of the site | `400 ALT_NOT_AT_SITE {id, alt, alleles}` |
| An allele `*`, or containing `N` | `400 UNSUPPORTED_ALLELE {allele}` |
| Manual: the region exists and the variant lies inside it | `404 UNKNOWN_REGION`; `400 REGION_OUT_OF_BOUNDS {region, position, length}` |
| REF is the genome's bases | `400 REF_MISMATCH {region, position, given, genome}` |
| The event slides at most 1,000 bp (`variation.max_shift`) | `400 VARIANT_TOO_REPETITIVE {region, position, shift, max}` |
| The known variants of the template window | by id, a 503 as above; a manual design degrades instead |
| A Primer3 slot | `503 BUSY` (5 s) |
| An orientation has room for the smallest product | `400 VARIANT_TOO_CLOSE_TO_END {region, position, region_length, needed}` |
| Primer3 and `ntthal` | `400 PRIMER3_INPUT_ERROR`; `500 PRIMER3_FAILED`, `THERMO_FAILED {binary}`; `503 PRIMER3_UNAVAILABLE`, `THERMO_UNAVAILABLE {binary, retry_after_s}` (60 s for a missing or non-executable binary, 5 s when it could not be started) |
| The 45 s deadline, at any step | `504 DEADLINE_EXCEEDED` |

`ntthal` is not checked at startup: a missing binary shows up as `503 THERMO_UNAVAILABLE` on the first design that needs it.

#### Response

| member | content |
| --- | --- |
| `variant` | The [variant](#variants), with `requested_id` and `submission_sequence`. A manual variant on a genome with variation data takes `ids`, `records`, `ems`, `consequence`, `alleles` and `multiallelic` from Ensembl (rs5413864115 below). |
| `template` | `{system_name, region, start, end, strand: 1, length, alt_length, seq, alt_seq, masked, mask_source, mask, masked_fraction, features}`. `seq` is the REF window and `alt_seq` the same window on the ALT haplotype. `features` are in template coordinates: the `variant` and `zone` spans, the `discriminating` indices, `alt_offset` (`len(alt) − len(ref)`), and the `exempt` intervals as `[start, length]`. |
| `assay` | The effective assay, plus `ems_target` (every record of the target is EMS) and `kasp_mix` (`kasp` only, else `null`). |
| `neighbours` | `{data: "ensembl" \| "none" \| "unavailable", window, variants, non_ems, ems, dense_non_ems}` for the template window, the target excluded. |
| `orientations` | `{forward, reverse}`, each `{status, reason, discriminating_position, relaxation_level, sets_found, blockers, attempts}`; `null` with `template_only`. |
| `sets` | The ranked sets, at most `num_sets`. |
| `check` | `{request, set_ids, unique_primers, omitted_set_ids}`; `null` when there is no set. |
| `settings` | `{preset, params, pinned, ladder, floors}`; `params` are level 0, with the product-range minimum raised. |
| `engine` | `{primer3, thermo, genotyping_design: "1", variation_source}`. `thermo` is `"ntthal "` plus the Primer3 version, because `ntthal` prints none. `variation_source` is `"ensembl 115"` or `null`. `primer3` and `thermo` are `null` with `template_only`. |
| `warnings` | Response-level [warnings](#genotyping-warning-codes). |

**Orientations**

| `status` | `reason` | meaning |
| --- | --- | --- |
| `ok` | `null` | Sets found; `relaxation_level` is the level that kept them. |
| `blocked` | `neighbour_at_3p` | A known non-EMS variant lies in the last 5 nt of the allele-specific primer (`genotyping.neighbour_3p_window`), listed in `blockers` with `distance_from_3p`. Primer3 is not run. Only natural targets with `neighbour_policy: "avoid_3p"` are blocked. |
| `skipped` | `too_close_to_end`, `n_in_primer_window`, `not_requested` | No room for the smallest product before the region end; an `N` in the allele-specific primer window; only the other orientation was requested. |
| `no_sets` | `null`, `below_floor`, `budget_exhausted` | Every allowed level was tried; `attempts[].explain` holds Primer3's reasons. |

`attempts[]` is `{level, changes, explain {left, right, pair}, pairs_returned, rejected {force, overlap,
common_neighbour_3p, alt_scoring_failed, below_floor, duplicate}, not_scored, sets}`. `pairs_returned` is always the sum of
`rejected`, `not_scored` and `sets`, and the `explain` lines are parsed as in `/primers/design`.

**Sets**

| member | content |
| --- | --- |
| `id`, `rank` | `S1`, `S2`, … by rank (`rank` is 0-based). |
| `key` | The first 12 hex characters of `sha256(orientation\|as_ref.target_seq\|as_alt.target_seq\|common.target_seq)`. It is stable across re-designs; match check results to sets by the target sequences, never by `id`. |
| `orientation`, `relaxation_level` | |
| `quality` | `poor` when any issue is high severity; `usable` when any issue is warn, or `relaxation_level ≥ 1`; else `good`. Quality does not change the order. |
| `score` | Lower is better. It sums: Primer3's pair penalty; 5 × the relaxation level; 2 × the AS Tm difference above 1.0 °C; 2 × the distance of `common_minus_as` outside −1…+3 °C; for `kasp`, `max(0, REF product − 100) / 25`; and 3 per high and 1 per warn issue, not counting `AS_TM_IMBALANCE`, `COMMON_TM_OUT_OF_RANGE` and info issues. The sum is exact, then rounded (spec §4.16). |
| `primers` | The `as_ref`, `as_alt` and `common` oligos, below. |
| `products` | `ref` and `alt`, each `{size, template, genomic, inserted_bases}`; an indel's ALT product differs in size by `alt_offset`. |
| `thermo` | `ref_common` and `alt_common` `{compl_any_th, compl_end_th}`, and `tailed` (the tailed cross-dimers, `null` without tails). |
| `tm_balance` | `{as_tm_diff, common_minus_as}`. |
| `neighbour_sites` | The distinct non-EMS known variants under the three primers. |
| `primer3_penalty` | The penalty of the Primer3 pair the set was built from, 4 decimals. |
| `issues`, `warnings` | `issues` holds every [set issue](#genotyping-warning-codes) with its `severity`; `warnings` holds the warn and high ones as `{code, message, details}`. |
| `check` | `{set: {id, ref_pair, alt_pair}, pairs}`: this set's part of a check request. |
| `order` | Three order rows, in the order REF, ALT, common. |

**Oligos** (`sets[].primers.*`)

| member | content |
| --- | --- |
| `role`, `allele`, `three_prime_base`, `haplotype` | `as_ref`, `as_alt` or `common`; the VCF allele the primer detects (`null` for the common primer, and `CGT` for the ALT primer of the insertion below); its 3′ base; `ref`, `alt` or `both`. |
| `target_seq` | The part that anneals, as ordered, including any deliberate mismatch. **It is the only sequence a check ever receives.** |
| `matched_seq` | The perfect match, before any deliberate mismatch. |
| `dye`, `tail_seq`, `order_seq` | `FAM`, `HEX` or `null`; the 5′ tail or `null`; `tail_seq + target_seq`, which is what the vendor synthesizes. |
| `len`, `order_len` | The lengths of `target_seq` and `order_seq`. |
| `tm`, `tm_method`, `matched_tm` | Primer3's Tm (`primer3`). For a deliberate-mismatch primer: the `ntthal` duplex Tm on its own allele (`ntthal_duplex`), with the perfect-match Tm in `matched_tm`. |
| `gc`, `hairpin_th`, `self_any_th`, `self_end_th`, `end_stability` | GC % computed from the sequence; the untailed structure Tm at 37 °C (negative values are reported as 0). |
| `primer3_problems` | Primer3's problem text for a derived oligo, e.g. `" Temperature too low;"`. It is informational; the floors decide. |
| `template`, `genomic`, `inserted_bases` | The span in `template.seq` or `alt_seq`; the genomic footprint (an ALT primer across a deletion has two `blocks`); the bases inside an insertion, which have no reference coordinate. |
| `deliberate_mismatch` | `null`, or `{position, original_base, new_base, template_base, terminal_pair, terminal_mismatch_class, added_mismatch_class, source}`. |
| `discrimination` | AS primers only. `{own_allele, other_allele}`, each `{mm_pos, likelihood}` as the check would classify the primer on that allele. `terminal_mismatch_class` is Little 1995's `max`, `strong`, `medium` or `weak`. `in_shift_tract` says whether the 3′ base lies in the indel's shift tract. |
| `tailed` | The tailed hairpin and self-dimer Tm; `null` without a tail. |
| `neighbours` | Known variants under the primer, EMS included, each with `distance_from_3p` (1 = the 3′ base). |

Temperatures and percentages have 2 decimals, rounded once, half away from zero, from exact decimals (spec §2.1).

**Neighbouring known variants.**

For natural targets with `neighbour_policy: "avoid_3p"`:
- In the last 5 nt of an allele-specific primer, a known variant blocks the orientation.
- In the last 5 nt of the common primer, it rejects that pair (`rejected.common_neighbour_3p`), and the orientation goes on.
- Elsewhere in a primer, it raises issue `NEIGHBOUR_IN_PRIMER` (warn). There is at most one such issue for the
  allele-specific pair and one for the common primer.

EMS neighbours never block and never raise an issue; they are listed in `neighbours` only.

For an EMS target (every record is EMS, warning `EMS_TARGET`), or with `neighbour_policy: "ignore"`, nothing blocks. A
non-EMS variant in a last-5-nt window becomes issue `NEIGHBOUR_AT_3P` (high) instead.

**Order sheet.**

- **`sets[].order[]`** rows are `{name, set_id, set_key, role, allele, dye, order_seq, target_seq, tail_seq, length, tm, gc,
  orientation, product_size_ref, product_size_alt, variant_key, notes}`.
  - `name` is `{label}_{set id}_{REF|ALT|COM}`, plus `_FAM` or `_HEX` for a tailed AS primer (`rs871475760_S1_REF_FAM`).
  - `notes` names a deliberate mismatch (`deliberate mismatch −2 A→G (Little 1995 Table 9.8.1)`) or a tailed hairpin issue
    (`tail hairpin 51.5 °C`).
- **`assay.kasp_mix`** is 12 µL of the REF primer, 12 µL of the ALT primer and 30 µL of the common primer, each at 100 µM,
  plus 46 µL of water. The ratio is from Makhoul et al. 2020; the water volume is inferred.
- **`variant.submission_sequence`** is 50 bp of reference on each side of `[REF/ALT]`, in minimal alleles with `-` written
  as nothing (`[/GT]`). Biallelic non-EMS SNV neighbours are written as IUPAC codes. Indel and multi-allelic neighbours
  stay reference bases and add warning `SUBMISSION_NEIGHBOURS_OMITTED {ids}`.

**The proposed check request.**

- **Body.** `check.request` is a ready `POST /primers/check` body: `{system_name, mode: "region", checks: ["specificity",
  "pangenome"], pairs, genotyping}`. It has no `genomes`, so every other assembly is checked, and no `params`.
- **Caps.** Sets are added in rank order while they fit 5 sets (`genotyping.check_max_sets`), 10 pairs (`check.max_pairs`)
  and 13 distinct primers (`genotyping.check_max_unique_primers`, the most that fit under the 6,000 CPU-s limit over the
  full panel). The rest are listed in `omitted_set_ids`.
- **Pairs.** Pair ids are `S1_REF`, `S1_ALT` and so on, `expected` is the REF product, and `genotyping.variant` is
  `variant.vcf` plus the region.

#### Example: rs871475760, KASP, two sets

The complete response. Only `template.seq` and `template.alt_seq` (801 nt each) are shortened, with `…`. It took 0.95 s.

<!-- example: request POST /primers/genotyping/design capture=capture-genotyping-design-rs871475760-kasp.json#/request/body -->
```json
{
  "system_name": "sorghum_bicolor", "variant": { "id": "rs871475760", "alt": "A" },
  "assay": { "type": "kasp", "num_sets": 2 }
}
```

<!-- example: response POST /primers/genotyping/design 200 capture=capture-genotyping-design-rs871475760-kasp.json#/response -->
```json
{
  "variant": {
    "key": "1:11109:C:A", "requested_id": "rs871475760", "ids": ["rs871475760"], "synonyms": ["tmp_1_11109_C_A"],
    "label": "1:11109 C/A", "kind": "snv", "region": "1", "vcf": { "position": 11109, "ref": "C", "alt": "A" },
    "minimal": { "start": 11109, "end": 11109, "ref": "C", "alt": "A" }, "alleles": ["C", "A"], "multiallelic": null,
    "shift": 0, "zone": { "start": 11109, "end": 11109 },
    "discriminating": {
      "forward": { "position": 11109, "ref_base": "C", "alt_base": "A", "alt_maps_to": 11109 },
      "reverse": { "position": 11109, "ref_base": "C", "alt_base": "A", "alt_maps_to": 11109 }
    },
    "records": [{ "id": "rs871475760", "source": "EVA", "ems": false }], "ems": false,
    "consequence": "downstream_gene_variant", "ref_verified": true, "designable": true, "issues": [],
    "submission_sequence": "YCCTCAAAAAGCTTCTCTAAGTGGTTATCCGAATATAGTCATACTCTATT[C/A]TGAATTTCTCGCTAGTCAAAGATAACAAAAATAGCATATTCTGGATTTCT"
  },
  "template": {
    "system_name": "sorghum_bicolor", "region": "1", "start": 10709, "end": 11509, "strand": 1, "length": 801,
    "alt_length": 801, "seq": "GCGAGTTCTCAAG…CATATGAT", "alt_seq": "GCGAGTTCTCAAG…CATATGAT", "masked": false,
    "mask_source": null, "mask": [], "masked_fraction": 0,
    "features": {
      "variant": { "start": 401, "end": 401 }, "zone": { "start": 401, "end": 401 },
      "discriminating": { "forward": 401, "reverse": 401 }, "alt_offset": 0, "exempt": [[365, 37], [401, 37]]
    }
  },
  "assay": {
    "type": "kasp", "orientation": "both", "tails": "ref_fam_alt_hex", "deliberate_mismatch": "none",
    "mismatch_position": 2, "num_sets": 2, "max_relaxation": 2, "neighbour_policy": "avoid_3p", "ems_target": false,
    "kasp_mix": {
      "stock_uM": 100, "as_ref_uL": 12, "as_alt_uL": 12, "common_uL": 30, "water_uL": 46, "total_uL": 100,
      "source": "Makhoul et al. 2020 (12:12:30 at 100 uM); water to 100 uL inferred"
    }
  },
  "neighbours": {
    "data": "ensembl", "window": { "start": 10709, "end": 11509 }, "variants": 57, "non_ems": 45, "ems": 12,
    "dense_non_ems": 0
  },
  "orientations": {
    "forward": {
      "status": "ok", "reason": null, "discriminating_position": 11109, "relaxation_level": 0, "sets_found": 8,
      "blockers": [],
      "attempts": [
        {
          "level": 0, "changes": {},
          "explain": {
            "left": {
              "raw": "considered 13, GC content failed 5, low tm 6, ok 2", "considered": 13, "GC content failed": 5,
              "low tm": 6, "ok": 2
            },
            "right": {
              "raw": "considered 4771, GC content failed 1718, low tm 1123, high tm 777, ok 1153", "considered": 4771,
              "GC content failed": 1718, "low tm": 1123, "high tm": 777, "ok": 1153
            },
            "pair": {
              "raw": "considered 1463, unacceptable product size 1443, ok 20", "considered": 1463,
              "unacceptable product size": 1443, "ok": 20
            }
          },
          "pairs_returned": 20,
          "rejected": {
            "force": 0, "overlap": 0, "common_neighbour_3p": 0, "alt_scoring_failed": 0, "below_floor": 0,
            "duplicate": 0
          },
          "not_scored": 12, "sets": 8
        }
      ]
    },
    "reverse": {
      "status": "ok", "reason": null, "discriminating_position": 11109, "relaxation_level": 0, "sets_found": 8,
      "blockers": [],
      "attempts": [
        {
          "level": 0, "changes": {},
          "explain": {
            "left": {
              "raw": "considered 4771, GC content failed 598, low tm 1550, high tm 1031, ok 1592", "considered": 4771,
              "GC content failed": 598, "low tm": 1550, "high tm": 1031, "ok": 1592
            },
            "right": { "raw": "considered 13, low tm 6, ok 7", "considered": 13, "low tm": 6, "ok": 7 },
            "pair": {
              "raw": "considered 905, unacceptable product size 883, tm diff too large 1, ok 21", "considered": 905,
              "unacceptable product size": 883, "tm diff too large": 1, "ok": 21
            }
          },
          "pairs_returned": 20,
          "rejected": {
            "force": 0, "overlap": 0, "common_neighbour_3p": 6, "alt_scoring_failed": 0, "below_floor": 0,
            "duplicate": 0
          },
          "not_scored": 6, "sets": 8
        }
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
          "target_seq": "ATCTTTGACTAGCGAGAAATTCAG", "matched_seq": "ATCTTTGACTAGCGAGAAATTCAG", "dye": "FAM",
          "tail_seq": "GAAGGTGACCAAGTTCATGCT", "order_seq": "GAAGGTGACCAAGTTCATGCTATCTTTGACTAGCGAGAAATTCAG",
          "len": 24, "order_len": 45, "tm": 57.1, "tm_method": "primer3", "matched_tm": null, "gc": 37.5,
          "hairpin_th": 0, "self_any_th": 0, "self_end_th": 0, "end_stability": 3.02, "primer3_problems": null,
          "template": { "start": 401, "end": 424, "sequence": "ref" },
          "genomic": {
            "region": "1", "start": 11109, "end": 11132, "strand": -1, "blocks": [{ "start": 11109, "end": 11132 }]
          },
          "inserted_bases": 0, "deliberate_mismatch": null,
          "discrimination": {
            "own_allele": { "mm_pos": [], "likelihood": "likely" },
            "other_allele": { "mm_pos": [1], "likelihood": "likely_weak" }, "terminal_mismatch_class": "max",
            "in_shift_tract": false
          },
          "tailed": { "hairpin_th": 42.62, "self_any_th": 3.63, "self_end_th": 11.64 }, "neighbours": []
        },
        "as_alt": {
          "role": "as_alt", "allele": "A", "three_prime_base": "T", "haplotype": "alt",
          "target_seq": "ATCTTTGACTAGCGAGAAATTCAT", "matched_seq": "ATCTTTGACTAGCGAGAAATTCAT", "dye": "HEX",
          "tail_seq": "GAAGGTCGGAGTCAACGGATT", "order_seq": "GAAGGTCGGAGTCAACGGATTATCTTTGACTAGCGAGAAATTCAT",
          "len": 24, "order_len": 45, "tm": 56.51, "tm_method": "primer3", "matched_tm": null, "gc": 33.33,
          "hairpin_th": 0, "self_any_th": 0, "self_end_th": 0, "end_stability": 2.57,
          "primer3_problems": " Temperature too low;", "template": { "start": 401, "end": 424, "sequence": "alt" },
          "genomic": {
            "region": "1", "start": 11109, "end": 11132, "strand": -1, "blocks": [{ "start": 11109, "end": 11132 }]
          },
          "inserted_bases": 0, "deliberate_mismatch": null,
          "discrimination": {
            "own_allele": { "mm_pos": [], "likelihood": "likely" },
            "other_allele": { "mm_pos": [1], "likelihood": "likely_weak" }, "terminal_mismatch_class": "max",
            "in_shift_tract": false
          },
          "tailed": { "hairpin_th": 51.5, "self_any_th": 26.85, "self_end_th": 11.47 }, "neighbours": []
        },
        "common": {
          "role": "common", "allele": null, "three_prime_base": "A", "haplotype": "both",
          "target_seq": "AGCTTCTCTAAGTGGTTATCCGA", "matched_seq": "AGCTTCTCTAAGTGGTTATCCGA", "dye": null,
          "tail_seq": null, "order_seq": "AGCTTCTCTAAGTGGTTATCCGA", "len": 23, "order_len": 23, "tm": 58.72,
          "tm_method": "primer3", "matched_tm": null, "gc": 43.48, "hairpin_th": 34.99, "self_any_th": 0,
          "self_end_th": 0, "end_stability": 4.55, "primer3_problems": null,
          "template": { "start": 360, "end": 382, "sequence": "ref" },
          "genomic": {
            "region": "1", "start": 11068, "end": 11090, "strand": 1, "blocks": [{ "start": 11068, "end": 11090 }]
          },
          "inserted_bases": 0, "deliberate_mismatch": null, "discrimination": null, "tailed": null,
          "neighbours": [
            {
              "key": "1:11069:G:A", "ids": ["tmp_1_11069_G_A"], "label": "1:11069 G/A", "start": 11069, "end": 11069,
              "alleles": "G/A", "ems": true, "distance_from_3p": 22
            }
          ]
        }
      },
      "products": {
        "ref": {
          "size": 65, "template": { "start": 360, "end": 424, "sequence": "ref" },
          "genomic": {
            "region": "1", "start": 11068, "end": 11132, "strand": 1, "blocks": [{ "start": 11068, "end": 11132 }]
          },
          "inserted_bases": 0
        },
        "alt": {
          "size": 65, "template": { "start": 360, "end": 424, "sequence": "alt" },
          "genomic": {
            "region": "1", "start": 11068, "end": 11132, "strand": 1, "blocks": [{ "start": 11068, "end": 11132 }]
          },
          "inserted_bases": 0
        }
      },
      "thermo": {
        "ref_common": { "compl_any_th": 0, "compl_end_th": 0 },
        "alt_common": { "compl_any_th": 0, "compl_end_th": 0 },
        "tailed": {
          "ref_alt_any_th": 6.89, "ref_alt_end_th": 8.21, "ref_common_any_th": 0, "ref_common_end_th": 0,
          "alt_common_any_th": 5.46, "alt_common_end_th": 0
        }
      },
      "tm_balance": { "as_tm_diff": 0.59, "common_minus_as": 1.62 }, "neighbour_sites": 0, "primer3_penalty": 7.1797,
      "warnings": [
        {
          "code": "ALT_PRIMER_SUBOPTIMAL",
          "message": "the derived as_alt primer would not have been chosen de novo (Primer3: temperature too low); it clears the hard floors",
          "details": { "oligo": "as_alt", "problems": " Temperature too low;" }
        },
        {
          "code": "TAILED_STRUCTURE", "message": "HEX-tailed as_alt: hairpin Tm 51.50 °C is above 47 °C",
          "details": { "oligo": "as_alt", "metric": "hairpin_th", "value": 51.5, "severity": "warn" }
        }
      ],
      "issues": [
        {
          "code": "ALT_PRIMER_SUBOPTIMAL", "severity": "warn",
          "message": "as_alt: Primer3 reports temperature too low",
          "details": { "oligo": "as_alt", "problems": " Temperature too low;" }
        },
        {
          "code": "TAILED_STRUCTURE", "severity": "warn", "message": "HEX-tailed as_alt hairpin 51.50 °C",
          "details": { "oligo": "as_alt", "metric": "hairpin_th", "value": 51.5 }
        }
      ],
      "check": {
        "set": { "id": "S1", "ref_pair": "S1_REF", "alt_pair": "S1_ALT" },
        "pairs": [
          {
            "id": "S1_REF", "left": "AGCTTCTCTAAGTGGTTATCCGA", "right": "ATCTTTGACTAGCGAGAAATTCAG",
            "expected": { "region": "1", "start": 11068, "end": 11132 }
          },
          {
            "id": "S1_ALT", "left": "AGCTTCTCTAAGTGGTTATCCGA", "right": "ATCTTTGACTAGCGAGAAATTCAT",
            "expected": { "region": "1", "start": 11068, "end": 11132 }
          }
        ]
      },
      "order": [
        {
          "name": "rs871475760_S1_REF_FAM", "set_id": "S1", "set_key": "f9df650ad116", "role": "as_ref",
          "allele": "C", "dye": "FAM", "order_seq": "GAAGGTGACCAAGTTCATGCTATCTTTGACTAGCGAGAAATTCAG",
          "target_seq": "ATCTTTGACTAGCGAGAAATTCAG", "tail_seq": "GAAGGTGACCAAGTTCATGCT", "length": 45, "tm": 57.1,
          "gc": 37.5, "orientation": "reverse", "product_size_ref": 65, "product_size_alt": 65,
          "variant_key": "1:11109:C:A", "notes": ""
        },
        {
          "name": "rs871475760_S1_ALT_HEX", "set_id": "S1", "set_key": "f9df650ad116", "role": "as_alt",
          "allele": "A", "dye": "HEX", "order_seq": "GAAGGTCGGAGTCAACGGATTATCTTTGACTAGCGAGAAATTCAT",
          "target_seq": "ATCTTTGACTAGCGAGAAATTCAT", "tail_seq": "GAAGGTCGGAGTCAACGGATT", "length": 45, "tm": 56.51,
          "gc": 33.33, "orientation": "reverse", "product_size_ref": 65, "product_size_alt": 65,
          "variant_key": "1:11109:C:A", "notes": "tail hairpin 51.5 °C"
        },
        {
          "name": "rs871475760_S1_COM", "set_id": "S1", "set_key": "f9df650ad116", "role": "common", "allele": null,
          "dye": null, "order_seq": "AGCTTCTCTAAGTGGTTATCCGA", "target_seq": "AGCTTCTCTAAGTGGTTATCCGA",
          "tail_seq": null, "length": 23, "tm": 58.72, "gc": 43.48, "orientation": "reverse", "product_size_ref": 65,
          "product_size_alt": 65, "variant_key": "1:11109:C:A", "notes": ""
        }
      ]
    },
    {
      "id": "S2", "key": "1accc54c262d", "rank": 1, "orientation": "forward", "relaxation_level": 0,
      "quality": "poor", "score": 19.63,
      "primers": {
        "as_ref": {
          "role": "as_ref", "allele": "C", "three_prime_base": "C", "haplotype": "ref",
          "target_seq": "GGTTATCCGAATATAGTCATACTCTATTC", "matched_seq": "GGTTATCCGAATATAGTCATACTCTATTC", "dye": "FAM",
          "tail_seq": "GAAGGTGACCAAGTTCATGCT", "order_seq": "GAAGGTGACCAAGTTCATGCTGGTTATCCGAATATAGTCATACTCTATTC",
          "len": 29, "order_len": 50, "tm": 57.28, "tm_method": "primer3", "matched_tm": null, "gc": 34.48,
          "hairpin_th": 39.2, "self_any_th": 12.01, "self_end_th": 12.01, "end_stability": 1.75,
          "primer3_problems": null, "template": { "start": 373, "end": 401, "sequence": "ref" },
          "genomic": {
            "region": "1", "start": 11081, "end": 11109, "strand": 1, "blocks": [{ "start": 11081, "end": 11109 }]
          },
          "inserted_bases": 0, "deliberate_mismatch": null,
          "discrimination": {
            "own_allele": { "mm_pos": [], "likelihood": "likely" },
            "other_allele": { "mm_pos": [1], "likelihood": "likely_weak" }, "terminal_mismatch_class": "max",
            "in_shift_tract": false
          },
          "tailed": { "hairpin_th": 36.16, "self_any_th": 3.63, "self_end_th": 12.01 }, "neighbours": []
        },
        "as_alt": {
          "role": "as_alt", "allele": "A", "three_prime_base": "A", "haplotype": "alt",
          "target_seq": "GGTTATCCGAATATAGTCATACTCTATTA", "matched_seq": "GGTTATCCGAATATAGTCATACTCTATTA", "dye": "HEX",
          "tail_seq": "GAAGGTCGGAGTCAACGGATT", "order_seq": "GAAGGTCGGAGTCAACGGATTGGTTATCCGAATATAGTCATACTCTATTA",
          "len": 29, "order_len": 50, "tm": 56.34, "tm_method": "primer3", "matched_tm": null, "gc": 31.03,
          "hairpin_th": 30.63, "self_any_th": 0.58, "self_end_th": 0, "end_stability": 0.98,
          "primer3_problems": " Temperature too low;", "template": { "start": 373, "end": 401, "sequence": "alt" },
          "genomic": {
            "region": "1", "start": 11081, "end": 11109, "strand": 1, "blocks": [{ "start": 11081, "end": 11109 }]
          },
          "inserted_bases": 0, "deliberate_mismatch": null,
          "discrimination": {
            "own_allele": { "mm_pos": [], "likelihood": "likely" },
            "other_allele": { "mm_pos": [1], "likelihood": "likely_weak" }, "terminal_mismatch_class": "max",
            "in_shift_tract": false
          },
          "tailed": { "hairpin_th": 66.22, "self_any_th": 27.34, "self_end_th": 10.26 }, "neighbours": []
        },
        "common": {
          "role": "common", "allele": null, "three_prime_base": "A", "haplotype": "both",
          "target_seq": "TCTTTGTCTACTGAGAAATCCAGA", "matched_seq": "TCTTTGTCTACTGAGAAATCCAGA", "dye": null,
          "tail_seq": null, "order_seq": "TCTTTGTCTACTGAGAAATCCAGA", "len": 24, "order_len": 24, "tm": 57.08,
          "tm_method": "primer3", "matched_tm": null, "gc": 37.5, "hairpin_th": 35.39, "self_any_th": 0,
          "self_end_th": 0, "end_stability": 3.86, "primer3_problems": null,
          "template": { "start": 441, "end": 464, "sequence": "ref" },
          "genomic": {
            "region": "1", "start": 11149, "end": 11172, "strand": -1, "blocks": [{ "start": 11149, "end": 11172 }]
          },
          "inserted_bases": 0, "deliberate_mismatch": null, "discrimination": null, "tailed": null,
          "neighbours": [
            {
              "key": "1:11161:A:T", "ids": ["rs872438201"], "label": "1:11161 A/T", "start": 11161, "end": 11161,
              "alleles": "A/T", "ems": false, "distance_from_3p": 13
            }
          ]
        }
      },
      "products": {
        "ref": {
          "size": 92, "template": { "start": 373, "end": 464, "sequence": "ref" },
          "genomic": {
            "region": "1", "start": 11081, "end": 11172, "strand": 1, "blocks": [{ "start": 11081, "end": 11172 }]
          },
          "inserted_bases": 0
        },
        "alt": {
          "size": 92, "template": { "start": 373, "end": 464, "sequence": "alt" },
          "genomic": {
            "region": "1", "start": 11081, "end": 11172, "strand": 1, "blocks": [{ "start": 11081, "end": 11172 }]
          },
          "inserted_bases": 0
        }
      },
      "thermo": {
        "ref_common": { "compl_any_th": 0, "compl_end_th": 0 },
        "alt_common": { "compl_any_th": 0, "compl_end_th": 0 },
        "tailed": {
          "ref_alt_any_th": 5.46, "ref_alt_end_th": 0, "ref_common_any_th": 0, "ref_common_end_th": 0,
          "alt_common_any_th": 0, "alt_common_end_th": 0
        }
      },
      "tm_balance": { "as_tm_diff": 0.95, "common_minus_as": -0.2 }, "neighbour_sites": 1, "primer3_penalty": 14.6344,
      "warnings": [
        {
          "code": "ALT_PRIMER_SUBOPTIMAL",
          "message": "the derived as_alt primer would not have been chosen de novo (Primer3: temperature too low); it clears the hard floors",
          "details": { "oligo": "as_alt", "problems": " Temperature too low;" }
        },
        {
          "code": "TAILED_STRUCTURE",
          "message": "HEX-tailed as_alt: hairpin Tm 66.22 °C is at or above 55 °C, the final KASP annealing temperature",
          "details": { "oligo": "as_alt", "metric": "hairpin_th", "value": 66.22, "severity": "high" }
        },
        {
          "code": "NEIGHBOUR_IN_PRIMER",
          "message": "known variant rs872438201 lies in the common primer, 13 nt from its 3′ end",
          "details": { "role": "common", "ids": ["rs872438201"], "distances": [13] }
        }
      ],
      "issues": [
        {
          "code": "ALT_PRIMER_SUBOPTIMAL", "severity": "warn",
          "message": "as_alt: Primer3 reports temperature too low",
          "details": { "oligo": "as_alt", "problems": " Temperature too low;" }
        },
        {
          "code": "TAILED_STRUCTURE", "severity": "high", "message": "HEX-tailed as_alt hairpin 66.22 °C",
          "details": { "oligo": "as_alt", "metric": "hairpin_th", "value": 66.22 }
        },
        {
          "code": "NEIGHBOUR_IN_PRIMER", "severity": "warn",
          "message": "rs872438201 in the common primer, 13 nt from its 3′ end",
          "details": { "role": "common", "ids": ["rs872438201"], "distances": [13] }
        }
      ],
      "check": {
        "set": { "id": "S2", "ref_pair": "S2_REF", "alt_pair": "S2_ALT" },
        "pairs": [
          {
            "id": "S2_REF", "left": "GGTTATCCGAATATAGTCATACTCTATTC", "right": "TCTTTGTCTACTGAGAAATCCAGA",
            "expected": { "region": "1", "start": 11081, "end": 11172 }
          },
          {
            "id": "S2_ALT", "left": "GGTTATCCGAATATAGTCATACTCTATTA", "right": "TCTTTGTCTACTGAGAAATCCAGA",
            "expected": { "region": "1", "start": 11081, "end": 11172 }
          }
        ]
      },
      "order": [
        {
          "name": "rs871475760_S2_REF_FAM", "set_id": "S2", "set_key": "1accc54c262d", "role": "as_ref",
          "allele": "C", "dye": "FAM", "order_seq": "GAAGGTGACCAAGTTCATGCTGGTTATCCGAATATAGTCATACTCTATTC",
          "target_seq": "GGTTATCCGAATATAGTCATACTCTATTC", "tail_seq": "GAAGGTGACCAAGTTCATGCT", "length": 50,
          "tm": 57.28, "gc": 34.48, "orientation": "forward", "product_size_ref": 92, "product_size_alt": 92,
          "variant_key": "1:11109:C:A", "notes": ""
        },
        {
          "name": "rs871475760_S2_ALT_HEX", "set_id": "S2", "set_key": "1accc54c262d", "role": "as_alt",
          "allele": "A", "dye": "HEX", "order_seq": "GAAGGTCGGAGTCAACGGATTGGTTATCCGAATATAGTCATACTCTATTA",
          "target_seq": "GGTTATCCGAATATAGTCATACTCTATTA", "tail_seq": "GAAGGTCGGAGTCAACGGATT", "length": 50,
          "tm": 56.34, "gc": 31.03, "orientation": "forward", "product_size_ref": 92, "product_size_alt": 92,
          "variant_key": "1:11109:C:A", "notes": "tail hairpin 66.2 °C"
        },
        {
          "name": "rs871475760_S2_COM", "set_id": "S2", "set_key": "1accc54c262d", "role": "common", "allele": null,
          "dye": null, "order_seq": "TCTTTGTCTACTGAGAAATCCAGA", "target_seq": "TCTTTGTCTACTGAGAAATCCAGA",
          "tail_seq": null, "length": 24, "tm": 57.08, "gc": 37.5, "orientation": "forward", "product_size_ref": 92,
          "product_size_alt": 92, "variant_key": "1:11109:C:A", "notes": ""
        }
      ]
    }
  ],
  "check": {
    "request": {
      "system_name": "sorghum_bicolor", "mode": "region", "checks": ["specificity", "pangenome"],
      "pairs": [
        {
          "id": "S1_REF", "left": "AGCTTCTCTAAGTGGTTATCCGA", "right": "ATCTTTGACTAGCGAGAAATTCAG",
          "expected": { "region": "1", "start": 11068, "end": 11132 }
        },
        {
          "id": "S1_ALT", "left": "AGCTTCTCTAAGTGGTTATCCGA", "right": "ATCTTTGACTAGCGAGAAATTCAT",
          "expected": { "region": "1", "start": 11068, "end": 11132 }
        },
        {
          "id": "S2_REF", "left": "GGTTATCCGAATATAGTCATACTCTATTC", "right": "TCTTTGTCTACTGAGAAATCCAGA",
          "expected": { "region": "1", "start": 11081, "end": 11172 }
        },
        {
          "id": "S2_ALT", "left": "GGTTATCCGAATATAGTCATACTCTATTA", "right": "TCTTTGTCTACTGAGAAATCCAGA",
          "expected": { "region": "1", "start": 11081, "end": 11172 }
        }
      ],
      "genotyping": {
        "variant": { "region": "1", "position": 11109, "ref": "C", "alt": "A" },
        "sets": [
          { "id": "S1", "ref_pair": "S1_REF", "alt_pair": "S1_ALT" },
          { "id": "S2", "ref_pair": "S2_REF", "alt_pair": "S2_ALT" }
        ]
      }
    },
    "set_ids": ["S1", "S2"], "unique_primers": 6, "omitted_set_ids": []
  },
  "settings": {
    "preset": "kasp",
    "params": {
      "opt_size": 22, "min_size": 18, "max_size": 30, "opt_tm": 60, "min_tm": 57, "max_tm": 63, "min_gc": 30,
      "max_gc": 70, "max_tm_diff": 3, "max_poly_x": 5, "product_size_ranges": [[61, 120]]
    },
    "pinned": [],
    "ladder": [
      {
        "level": 1,
        "changes": {
          "max_size": 32, "min_tm": 55, "max_tm": 65, "min_gc": 20, "max_gc": 80, "product_size_ranges": [[65, 150]]
        }
      },
      { "level": 2, "changes": { "min_tm": 52, "max_tm_diff": 6 } }
    ],
    "floors": { "as_min_tm": 52, "as_min_gc": 15 }
  },
  "engine": {
    "primer3": "2.6.1", "thermo": "ntthal 2.6.1", "genotyping_design": "1", "variation_source": "ensembl 115"
  },
  "warnings": []
}
```

#### Example: tmp_1_11502_C_CGT, forward orientation blocked

An insertion, designed by id:
- Two Ensembl ids describe the event, and they are merged.
- The forward orientation is blocked by rs5413863234, 4 nt from the allele-specific 3′ end, so the set comes from the
  reverse orientation.
- Its ALT primer ends on an inserted base (`inserted_bases: 1`, `allele: "CGT"`).

This is an excerpt: `"…": "…"` marks where members were left out.

<!-- example: request POST /primers/genotyping/design capture=capture-genotyping-design-tmp_1_11502_C_CGT.json#/request/body -->
```json
{ "system_name": "sorghum_bicolor", "variant": { "id": "tmp_1_11502_C_CGT" }, "assay": { "num_sets": 1 } }
```

<!-- example: response POST /primers/genotyping/design 200 capture=capture-genotyping-design-tmp_1_11502_C_CGT.json#/response -->
```json
{
  "variant": {
    "key": "1:11502:C:CGT", "requested_id": "tmp_1_11502_C_CGT", "ids": ["tmp_1_11502_C_CGT", "rs5413863549"],
    "label": "1:11502^11503 -/GT", "kind": "insertion", "vcf": { "position": 11502, "ref": "C", "alt": "CGT" },
    "minimal": { "start": 11503, "end": 11502, "ref": "-", "alt": "GT" }, "shift": 0,
    "zone": { "start": 11502, "end": 11503 },
    "discriminating": {
      "forward": { "position": 11503, "ref_base": "A", "alt_base": "G", "alt_maps_to": null },
      "reverse": { "position": 11502, "ref_base": "C", "alt_base": "T", "alt_maps_to": null }
    },
    "records": [
      { "id": "tmp_1_11502_C_CGT", "source": "SAP_PMID35653240_Boatwri", "ems": false },
      { "id": "rs5413863549", "source": "EVA", "ems": false }
    ],
    "…": "…"
  },
  "template": {
    "region": "1", "start": 11102, "end": 11903, "length": 802, "alt_length": 804,
    "features": {
      "variant": { "start": 401, "end": 401 }, "zone": { "start": 401, "end": 402 },
      "discriminating": { "forward": 402, "reverse": 401 }, "alt_offset": 2, "exempt": [[366, 37], [401, 37]]
    },
    "…": "…"
  },
  "neighbours": {
    "data": "ensembl", "window": { "start": 11102, "end": 11903 }, "variants": 59, "non_ems": 49, "ems": 10,
    "dense_non_ems": 5
  },
  "orientations": {
    "forward": {
      "status": "blocked", "reason": "neighbour_at_3p", "discriminating_position": 11503, "relaxation_level": null,
      "sets_found": 0,
      "blockers": [
        {
          "key": "1:11500:G:A", "ids": ["rs5413863234"], "label": "1:11500 G/A", "start": 11500, "end": 11500,
          "alleles": "G/A", "ems": false, "distance_from_3p": 4
        }
      ],
      "attempts": []
    },
    "reverse": {
      "status": "ok", "reason": null, "discriminating_position": 11502, "relaxation_level": 0, "sets_found": 8,
      "blockers": [], "…": "…"
    }
  },
  "sets": [
    {
      "id": "S1", "key": "53942cb55348", "rank": 0, "orientation": "reverse", "relaxation_level": 0,
      "quality": "usable", "score": 9.56,
      "primers": {
        "as_ref": {
          "role": "as_ref", "allele": "C", "three_prime_base": "G", "target_seq": "GCAGGAAAAGAAATCCTAACATCATATG",
          "dye": "FAM", "len": 28, "tm": 59.19, "gc": 35.71,
          "genomic": {
            "region": "1", "start": 11502, "end": 11529, "strand": -1, "blocks": [{ "start": 11502, "end": 11529 }]
          },
          "inserted_bases": 0,
          "neighbours": [
            {
              "key": "1:11509:T:C", "ids": ["rs5413864238"], "label": "1:11509 T/C", "start": 11509, "end": 11509,
              "alleles": "T/C", "ems": false, "distance_from_3p": 8
            },
            {
              "key": "1:11516:A:T", "ids": ["rs5413863874"], "label": "1:11516 A/T", "start": 11516, "end": 11516,
              "alleles": "A/T", "ems": false, "distance_from_3p": 15
            },
            {
              "key": "1:11523:T:A", "ids": ["rs5413863270"], "label": "1:11523 T/A", "start": 11523, "end": 11523,
              "alleles": "T/A", "ems": false, "distance_from_3p": 22
            }
          ],
          "…": "…"
        },
        "as_alt": {
          "role": "as_alt", "allele": "CGT", "three_prime_base": "A", "target_seq": "GCAGGAAAAGAAATCCTAACATCATATA",
          "dye": "HEX", "len": 28, "tm": 58.03, "gc": 32.14,
          "genomic": {
            "region": "1", "start": 11503, "end": 11529, "strand": -1, "blocks": [{ "start": 11503, "end": 11529 }]
          },
          "inserted_bases": 1, "…": "…"
        },
        "common": {
          "role": "common", "target_seq": "AGGATCTTTGCAACCCTGTGTT", "len": 22, "tm": 60.43, "gc": 45.45,
          "genomic": {
            "region": "1", "start": 11458, "end": 11479, "strand": 1, "blocks": [{ "start": 11458, "end": 11479 }]
          },
          "neighbours": [
            {
              "key": "1:11474:T:A", "ids": ["rs5413863452"], "label": "1:11474 T/A", "start": 11474, "end": 11474,
              "alleles": "T/A", "ems": false, "distance_from_3p": 6
            }
          ],
          "…": "…"
        }
      },
      "products": {
        "ref": { "size": 72, "inserted_bases": 0, "…": "…" }, "alt": { "size": 74, "inserted_bases": 2, "…": "…" }
      },
      "tm_balance": { "as_tm_diff": 1.16, "common_minus_as": 1.24 }, "neighbour_sites": 4,
      "warnings": [
        {
          "code": "AS_TM_IMBALANCE", "message": "as_ref and as_alt Tm differ by 1.16 °C (limit 1.0)",
          "details": { "diff": 1.16 }
        },
        {
          "code": "NEIGHBOUR_IN_PRIMER",
          "message": "known variants lie in the allele-specific primers: rs5413864238 (8 nt), rs5413863874 (15 nt), rs5413863270 (22 nt) from the 3′ end",
          "details": {
            "role": "as_ref", "ids": ["rs5413864238", "rs5413863874", "rs5413863270"], "distances": [8, 15, 22]
          }
        },
        {
          "code": "NEIGHBOUR_IN_PRIMER",
          "message": "known variant rs5413863452 lies in the common primer, 6 nt from its 3′ end",
          "details": { "role": "common", "ids": ["rs5413863452"], "distances": [6] }
        }
      ],
      "check": {
        "set": { "id": "S1", "ref_pair": "S1_REF", "alt_pair": "S1_ALT" },
        "pairs": [
          {
            "id": "S1_REF", "left": "AGGATCTTTGCAACCCTGTGTT", "right": "GCAGGAAAAGAAATCCTAACATCATATG",
            "expected": { "region": "1", "start": 11458, "end": 11529 }
          },
          {
            "id": "S1_ALT", "left": "AGGATCTTTGCAACCCTGTGTT", "right": "GCAGGAAAAGAAATCCTAACATCATATA",
            "expected": { "region": "1", "start": 11458, "end": 11529 }
          }
        ]
      },
      "…": "…"
    }
  ],
  "check": { "set_ids": ["S1"], "unique_primers": 3, "omitted_set_ids": [], "…": "…" },
  "warnings": [
    {
      "code": "DUPLICATE_VARIANT_IDS",
      "message": "2 Ensembl ids describe the same event 1:11502:C:CGT; they were merged",
      "details": { "key": "1:11502:C:CGT", "ids": ["tmp_1_11502_C_CGT", "rs5413863549"] }
    },
    {
      "code": "ORIENTATION_BLOCKED",
      "message": "forward: known variant rs5413863234 (G/A) lies 4 nt from the allele-specific primer's 3′ end",
      "details": { "orientation": "forward", "ids": ["rs5413863234"], "distances": [4] }
    },
    {
      "code": "DENSE_NEIGHBOURS", "message": "5 known non-EMS variants lie within 30 bp of the variant",
      "details": { "count": 5, "window": 30 }
    }
  ],
  "…": "…"
}
```

#### Example: a manual variant

rs5413864115 entered by hand in Ensembl style, as `1:11283 A/-`: a deletion inside `AAA` that can slide 2 bases.
- The variant is left-aligned to `1:11282:CA:C` and takes its id from Ensembl; `requested_id` is `null`.
- Both orientations need relaxation level 1.
- The ALT products are 1 bp shorter (`alt_offset: -1`), and the forward ALT primer spans the deletion (two `genomic.blocks`).

This is an excerpt.

<!-- example: request POST /primers/genotyping/design capture=capture-genotyping-design-manual-deletion.json#/request/body -->
```json
{
  "system_name": "sorghum_bicolor", "variant": { "region": "1", "position": 11283, "ref": "A", "alt": "-" },
  "assay": { "num_sets": 2 }
}
```

<!-- example: response POST /primers/genotyping/design 200 capture=capture-genotyping-design-manual-deletion.json#/response -->
```json
{
  "variant": {
    "key": "1:11282:CA:C", "requested_id": null, "ids": ["rs5413864115"], "label": "1:11283-11283 A/-",
    "kind": "deletion", "vcf": { "position": 11282, "ref": "CA", "alt": "C" },
    "minimal": { "start": 11283, "end": 11283, "ref": "A", "alt": "-" }, "shift": 2,
    "zone": { "start": 11282, "end": 11286 },
    "discriminating": {
      "forward": { "position": 11285, "ref_base": "A", "alt_base": "G", "alt_maps_to": 11286 },
      "reverse": { "position": 11283, "ref_base": "A", "alt_base": "C", "alt_maps_to": 11282 }
    },
    "records": [{ "id": "rs5413864115", "source": "EVA", "ems": false }],
    "submission_sequence": "AAAATATATAGAAAACAATTTTATACAGATGATTTTCCAAATGATGATTC[A/]AAGTGTGAAATTTGRAAAGWCTCTTRGASATGMTYTAAGTGGAAGGAACA",
    "…": "…"
  },
  "template": {
    "start": 10883, "end": 11685, "length": 803, "alt_length": 802,
    "features": {
      "variant": { "start": 400, "end": 401 }, "zone": { "start": 400, "end": 404 },
      "discriminating": { "forward": 403, "reverse": 401 }, "alt_offset": -1, "exempt": [[367, 38], [400, 38]]
    },
    "…": "…"
  },
  "neighbours": {
    "data": "ensembl", "window": { "start": 10883, "end": 11685 }, "variants": 53, "non_ems": 45, "ems": 8,
    "dense_non_ems": 4
  },
  "orientations": {
    "forward": {
      "status": "ok", "reason": null, "discriminating_position": 11285, "relaxation_level": 1, "sets_found": 8,
      "…": "…"
    },
    "reverse": {
      "status": "ok", "reason": null, "discriminating_position": 11283, "relaxation_level": 1, "sets_found": 8,
      "…": "…"
    }
  },
  "sets": [
    {
      "id": "S1", "key": "a5277232d8ab", "orientation": "forward", "relaxation_level": 1, "quality": "usable",
      "score": 15.48,
      "primers": {
        "as_ref": {
          "target_seq": "ACAGATGATTTTCCAAATGATGATTCAAA", "len": 29, "tm": 59.22,
          "genomic": {
            "region": "1", "start": 11257, "end": 11285, "strand": 1, "blocks": [{ "start": 11257, "end": 11285 }]
          },
          "discrimination": {
            "own_allele": { "mm_pos": [], "likelihood": "likely" },
            "other_allele": { "mm_pos": [1], "likelihood": "likely_weak" }, "terminal_mismatch_class": "weak",
            "in_shift_tract": true
          },
          "…": "…"
        },
        "as_alt": {
          "target_seq": "ACAGATGATTTTCCAAATGATGATTCAAG", "len": 29, "tm": 59.53,
          "genomic": {
            "region": "1", "start": 11257, "end": 11286, "strand": 1,
            "blocks": [{ "start": 11257, "end": 11282 }, { "start": 11284, "end": 11286 }]
          },
          "inserted_bases": 0, "…": "…"
        },
        "common": {
          "target_seq": "CCCCATGTTTTTGTTCCTTCCA", "tm": 59.3,
          "genomic": {
            "region": "1", "start": 11323, "end": 11344, "strand": -1, "blocks": [{ "start": 11323, "end": 11344 }]
          },
          "…": "…"
        }
      },
      "products": { "ref": { "size": 88, "…": "…" }, "alt": { "size": 87, "inserted_bases": 0, "…": "…" } },
      "warnings": [
        {
          "code": "NEIGHBOUR_IN_PRIMER",
          "message": "known variants lie in the common primer: tmp_1_11334_A_G (12 nt), tmp_1_11335_A_G (13 nt), tmp_1_11340_TG_T (19 nt) from the 3′ end",
          "details": {
            "role": "common", "ids": ["tmp_1_11334_A_G", "tmp_1_11335_A_G", "tmp_1_11340_TG_T"],
            "distances": [12, 13, 19]
          }
        },
        {
          "code": "SHIFT_TRACT_DISCRIMINATION",
          "message": "the as_ref primer's 3′ base lies inside the shift tract of the indel; on the other allele it may still prime through a bulge",
          "details": { "primer": "as_ref", "shift": 2 }
        }
      ],
      "…": "…"
    },
    {
      "id": "S2", "key": "cb3ef66afd37", "orientation": "reverse", "relaxation_level": 1, "quality": "usable",
      "score": 23.95,
      "primers": {
        "as_ref": {
          "target_seq": "AGAGTCTTTTCAAATTTCACACTTT", "len": 25, "tm": 56.09,
          "genomic": {
            "region": "1", "start": 11283, "end": 11307, "strand": -1, "blocks": [{ "start": 11283, "end": 11307 }]
          },
          "discrimination": {
            "own_allele": { "mm_pos": [], "likelihood": "likely" },
            "other_allele": { "mm_pos": [1], "likelihood": "likely_weak" }, "terminal_mismatch_class": "max",
            "in_shift_tract": true
          },
          "…": "…"
        },
        "as_alt": {
          "target_seq": "AGAGTCTTTTCAAATTTCACACTTG", "len": 25, "tm": 56.71,
          "genomic": {
            "region": "1", "start": 11282, "end": 11307, "strand": -1,
            "blocks": [{ "start": 11282, "end": 11282 }, { "start": 11284, "end": 11307 }]
          },
          "inserted_bases": 0, "…": "…"
        },
        "common": {
          "target_seq": "ACAAAAATAGCTCTCTAGAGTATACACA", "tm": 58.12,
          "genomic": {
            "region": "1", "start": 11179, "end": 11206, "strand": 1, "blocks": [{ "start": 11179, "end": 11206 }]
          },
          "…": "…"
        }
      },
      "products": { "ref": { "size": 129, "…": "…" }, "alt": { "size": 128, "inserted_bases": 0, "…": "…" } },
      "warnings": [
        {
          "code": "NEIGHBOUR_IN_PRIMER",
          "message": "known variants lie in the allele-specific primers: tmp_1_11298_A_G (16 nt), tmp_1_11303_A_T (21 nt) from the 3′ end",
          "details": { "role": "as_ref", "ids": ["tmp_1_11298_A_G", "tmp_1_11303_A_T"], "distances": [16, 21] }
        },
        {
          "code": "NEIGHBOUR_IN_PRIMER",
          "message": "known variant rs873026643 lies in the common primer, 25 nt from its 3′ end",
          "details": { "role": "common", "ids": ["rs873026643"], "distances": [25] }
        },
        {
          "code": "SHIFT_TRACT_DISCRIMINATION",
          "message": "the as_ref primer's 3′ base lies inside the shift tract of the indel; on the other allele it may still prime through a bulge",
          "details": { "primer": "as_ref", "shift": 2 }
        }
      ],
      "…": "…"
    }
  ],
  "check": {
    "request": {
      "system_name": "sorghum_bicolor", "mode": "region", "checks": ["specificity", "pangenome"],
      "genotyping": {
        "variant": { "region": "1", "position": 11282, "ref": "CA", "alt": "C" },
        "sets": [
          { "id": "S1", "ref_pair": "S1_REF", "alt_pair": "S1_ALT" },
          { "id": "S2", "ref_pair": "S2_REF", "alt_pair": "S2_ALT" }
        ]
      },
      "…": "…"
    },
    "set_ids": ["S1", "S2"], "unique_primers": 6, "…": "…"
  },
  "warnings": [
    {
      "code": "SHIFTABLE_INDEL",
      "message": "the deletion can slide 2 bases inside AAA; the forward and reverse primers end at different bases and the wrong-allele primer may still prime through a 1-nt bulge",
      "details": {
        "shift": 2, "forward_position": 11285, "reverse_position": 11283, "zone": { "start": 11282, "end": 11286 }
      }
    },
    {
      "code": "DENSE_NEIGHBOURS", "message": "4 known non-EMS variants lie within 30 bp of the variant",
      "details": { "count": 4, "window": 30 }
    },
    {
      "code": "RELAXED_CONSTRAINTS", "message": "forward: sets need relaxation level 1",
      "details": {
        "orientation": "forward", "level": 1,
        "changes": {
          "max_size": 32, "min_tm": 55, "max_tm": 65, "min_gc": 20, "max_gc": 80, "product_size_ranges": [[65, 150]]
        }
      }
    },
    {
      "code": "RELAXED_CONSTRAINTS", "message": "reverse: sets need relaxation level 1",
      "details": {
        "orientation": "reverse", "level": 1,
        "changes": {
          "max_size": 32, "min_tm": 55, "max_tm": 65, "min_gc": 20, "max_gc": 80, "product_size_ranges": [[65, 150]]
        }
      }
    }
  ],
  "…": "…"
}
```

#### Error example

A wrong reference allele, with `template_only` (a full design gets the same answer):

<!-- example: request POST /primers/genotyping/design capture=capture-genotyping-design-ref-mismatch.json#/request/body -->
```json
{
  "system_name": "sorghum_bicolor", "variant": { "region": "1", "position": 11109, "ref": "A", "alt": "C" },
  "template_only": true
}
```

<!-- example: response POST /primers/genotyping/design 400 capture=capture-genotyping-design-ref-mismatch.json#/response -->
```json
{
  "message": "the reference allele A does not match the genome base C at 1:11109", "code": "REF_MISMATCH",
  "details": { "region": "1", "position": 11109, "given": "A", "genome": "C" }
}
```

### `POST /primers/check`: the `genotyping` block

An optional member of the check request that says which pairs form allele-specific sets for a variant. Each set is two
ordinary `pairs[]` items that share the common primer: the REF pair (the REF-specific primer and the common primer) and
the ALT pair (the ALT-specific primer and the common primer). This subsection covers the request; the results are
[below](#results-resultsgenotyping).

| field | type | notes |
| --- | --- | --- |
| `genotyping.variant` | `{region, position, ref, alt}`, all required | `region` ≤ 255, `position` ≥ 1, `ref` and `alt` `^[ACGTacgt]{1,50}$` (VCF style, no `-`). It need not be left-aligned, and `ref` and `alt` must differ |
| `genotyping.sets` | 1–5 items | |
| `genotyping.sets[].id`, `.ref_pair`, `.alt_pair` | `^[A-Za-z0-9_.:-]{1,64}$`, all required | set ids are unique; `ref_pair` and `alt_pair` name `pairs[]` items |

**Submit-time validation.** A mislabelled or misplaced set would queue a job of up to 6,000 CPU-s and return confident,
wrong allele predictions, so it is refused at submit. The rules run in this order:

1. The block's shape (swagger, then the handler): `400 INVALID_REQUEST {field}`. Beyond the swagger definition, `ref` and
   `alt` must differ and set ids must be unique.
2. The existing pair rules (unique pair ids, at most 20 distinct primers). The 20-primer cap cannot bind a genotyping
   check: rule 4's `no_shared_common` makes a set exactly three distinct primers (two allele-specific plus one common),
   so *N* sets are 2*N* pairs and at most 3*N* distinct primers — at the 5-set ceiling, 10 pairs (the pair cap itself)
   and at most 15 primers. What can bind is the 5-set cap, the 10-pair cap and the cost guard of rule 7. A client
   mirroring these limits should name whichever one it is blocking on.
3. `mode` must be `gene` or `region`, the modes that honour `expected`: `400 INVALID_REQUEST {field: "genotyping",
   reason: "mode"}`.
4. The links between sets and pairs, without I/O: `400 GENOTYPING_SET_INVALID` with reasons `unknown_pair`, `same_pair`,
   `pair_reused`, `expected_required`, `expected_differs`, `no_shared_common`.
5. The existing rules: `expected` and `params`, the catalog, the gene, the reference assembly and its BLAST DB.
6. One FASTA read of `position ± 1,700` bp:
   - `404 UNKNOWN_REGION`, `422 NO_SEQUENCE`;
   - `400 REGION_OUT_OF_BOUNDS {region, start, end, length}`, also when left-aligning reaches the region edge;
   - `400 REF_MISMATCH {region, position, given, genome}`;
   - `400 VARIANT_TOO_REPETITIVE {region, position, shift, max}`, with `shift: null` when the event slides beyond the read;
   - then, per set, `400 GENOTYPING_SET_INVALID` with reasons `not_at_variant`, `alleles_swapped`, `too_many_edits`,
     `common_in_zone`, `expected_mismatch`.
7. The pan-genome genomes and the cost guard (`422 JOB_TOO_LARGE`).

| `reason` | the set is refused when | extra `details` |
| --- | --- | --- |
| `unknown_pair` | `ref_pair` or `alt_pair` names no pair | `pair_ids` are the missing ids |
| `same_pair` | `ref_pair` equals `alt_pair` | |
| `pair_reused` | a pair already belongs to an earlier set | |
| `expected_required` | a pair has no `expected` | `pair_ids` are the pairs without it |
| `expected_differs` | the two pairs' `expected` differ | |
| `no_shared_common` | the two pairs do not share exactly one primer, on the same side | |
| `not_at_variant` | an allele-specific primer does not end at the discriminating base on its own haplotype | `primer` (`ref_pair` or `alt_pair`), `sequence` |
| `alleles_swapped` | an allele-specific primer matches the other allele at least as well as its own | `primer`, `sequence` |
| `too_many_edits` | an allele-specific primer differs from its own allele by more than one mismatch 2 or 3 nt from its 3′ end, or by an indel | `primer`, `sequence`, `mm_pos` |
| `common_in_zone` | the common primer, placed from `expected`, does not lie beyond the zone and the base beside it | `primer: "common"`, `sequence` |
| `expected_mismatch` | `expected` is on another region, or the REF-specific primer's 5′ end is not `expected.start` (forward) or `expected.end` (reverse) | `primer: "ref_pair"`, `sequence` |

Every rejection has `details {set_id, reason, pair_ids}` plus the members listed.

- **Orientation** is derived from the shared primer and is never sent. A shared right primer makes the left primers
  allele-specific (`forward`); a shared left primer makes the right ones allele-specific (`reverse`).
- **Deliberate mismatches** need no declaration: one extra mismatch 2 or 3 nt from the 3′ end of an allele-specific
  primer is accepted as one, and its position is kept for the worker.
- **The stored request** keeps only the client's keys, with the variant uppercased and left-aligned, so `job.request`
  stays a valid check request. `POST /primers/check` never calls Ensembl.
- **Job ids.** A request with `genotyping` hashes with algorithm version `"2+g1"` (`jobId`'s `algo`) instead of `"2"`.
  Every existing job id is therefore unchanged, and a change to the allele caller (`g1`) invalidates genotyping jobs only.
  The check's `ALGORITHM_VERSION` stays `"2"`.

**Cost.** The allele caller adds `(1 + pan-genome genomes) × 0.2` CPU-s (`check.genotype_cpu_s_per_genome`) to the
[estimate](#cost-guard), and `breakdown.genotyping` appears only in such requests. The table below was computed with
`check/cost.js` and the real genome sizes (reference 708,735,318 bases; the 119 other assemblies 83,774,241,795):

| Check (region mode) | Without `genotyping` | With `genotyping` |
| --- | --- | --- |
| 1 set (3 primers), all 119 other assemblies | 1,333 | 1,357 |
| 2 sets (6 primers), all 119 | 2,666 | 2,690 |
| 4 sets (12 primers), all 119 | 5,332 | 5,356 |
| 13 distinct primers, all 119 | 5,776 | 5,800, the most accepted |
| 14 distinct primers, all 119 | 6,221 | 6,245 → `422 JOB_TOO_LARGE` |
| 2 sets (6 primers), 3 genomes | 95 | 95 |

The UI mirrors this term for display. A real job's measured CPU, and the allele caller's own cost, are under
[Measured cost](#measured-cost).

#### Example

The `check.request` of the rs871475760 design above, unchanged:

<!-- example: request POST /primers/check capture=capture-genotyping-design-rs871475760-kasp.json#/response/check/request -->
```json
{
  "system_name": "sorghum_bicolor", "mode": "region", "checks": ["specificity", "pangenome"],
  "pairs": [
    {
      "id": "S1_REF", "left": "AGCTTCTCTAAGTGGTTATCCGA", "right": "ATCTTTGACTAGCGAGAAATTCAG",
      "expected": { "region": "1", "start": 11068, "end": 11132 }
    },
    {
      "id": "S1_ALT", "left": "AGCTTCTCTAAGTGGTTATCCGA", "right": "ATCTTTGACTAGCGAGAAATTCAT",
      "expected": { "region": "1", "start": 11068, "end": 11132 }
    },
    {
      "id": "S2_REF", "left": "GGTTATCCGAATATAGTCATACTCTATTC", "right": "TCTTTGTCTACTGAGAAATCCAGA",
      "expected": { "region": "1", "start": 11081, "end": 11172 }
    },
    {
      "id": "S2_ALT", "left": "GGTTATCCGAATATAGTCATACTCTATTA", "right": "TCTTTGTCTACTGAGAAATCCAGA",
      "expected": { "region": "1", "start": 11081, "end": 11172 }
    }
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

**Response `202`.** Recorded on 2026-09-15 from the dev API and its worker on 127.0.0.1:50112 (worktree HEAD 228a76b). The
request was the one above with `genomes` added, listing the 11 research assemblies:
`sorghum_bicolorv5`, `sorghum_austrcf317961`, `sorghum_pi565121`, `sorghum_pi655972`, `sorghum_pi180348`,
`sorghum_pi276837`, `sorghum_pi656027`, `sorghum_pi510757`, `sorghum_rio`, `sorghum_pi329250`, `sorghum_s3691`.
- `progress.total` is 12 genome tasks: the reference and the 11 assemblies.
- The same pairs without `genotyping` estimate 270 CPU-s.
- The job's results are [documented below](#example-rs871475760-sets-s1-and-s2-over-11-assemblies).

<!-- example: response POST /primers/check 202 capture=capture-check-genotyping-submit.json#/response -->
```json
{
  "job_id": "8e9160d598f602137d94efa9fc409264", "status": "queued", "kind": "pangenome", "queue_position": 0,
  "progress": { "done": 0, "total": 12, "stage": "queued", "running": [] }, "estimate": { "cpu_s": 273 },
  "created_at": "2026-09-15T22:54:34.441Z", "warnings": []
}
```

**Rejections.** These were produced by the submit-time validator (`check/genotype.js`) on the real sorghum_bicolor FASTA,
without submitting a job. With `S1_ALT.left` changed to `GGTTATCCGAATATAGTCATACTCTATTC`, S1's pairs no longer share a
primer:

<!-- example: response POST /primers/check 400 capture=capture-check-genotyping-no-shared-common.json#/response -->
```json
{
  "message": "set S1: S1_REF and S1_ALT must share exactly one primer, the common primer, on the same side",
  "code": "GENOTYPING_SET_INVALID",
  "details": { "set_id": "S1", "reason": "no_shared_common", "pair_ids": ["S1_REF", "S1_ALT"] }
}
```

With `S1_REF.right` changed to `TCTTTGACTAGCGAGAAATTCAGA` (the primer Primer3 picks at 1:11108–11131 without the zone
guard):

<!-- example: response POST /primers/check 400 capture=capture-check-genotyping-not-at-variant.json#/response -->
```json
{
  "message": "set S1: the REF-specific primer TCTTTGACTAGCGAGAAATTCAGA does not end at the variant (1:11109)",
  "code": "GENOTYPING_SET_INVALID",
  "details": {
    "set_id": "S1", "reason": "not_at_variant", "pair_ids": ["S1_REF"], "primer": "ref_pair",
    "sequence": "TCTTTGACTAGCGAGAAATTCAGA"
  }
}
```

### Genotyping error codes

These are all handler errors `{message, code, details}`; no new HTTP status is used.

| Status | Code | Endpoints | When | `details` |
| --- | --- | --- | --- | --- |
| 400 | `VARIANT_WINDOW_TOO_LONG` | variants list | the window is over 50,000 bp | `{length, max}` |
| 400 | `INVALID_VARIANT` | design | a wrong combination of variant fields, or invalid or over-long alleles | `{reason: "id_or_manual" \| "alleles" \| "allele_too_long"}` |
| 400 | `REF_MISMATCH` | design, check | the reference allele is not the genome's bases | `{region, position, given, genome}` |
| 400 | `ALT_REQUIRED` | design | an id with several designable alternative alleles, and no `alt` | `{id, alts}` |
| 400 | `ALT_NOT_AT_SITE` | design | `alt` is not an allele of the id's site | `{id, alt, alleles}` |
| 400 | `UNSUPPORTED_ALLELE` | design | an allele `*`, or containing `N` | `{allele}` |
| 400 | `VARIANT_TOO_CLOSE_TO_END` | design | no orientation has room for the smallest product | `{region, position, region_length, needed}` |
| 400 | `VARIANT_TOO_REPETITIVE` | design, check | the event slides more than 1,000 bp | `{region, position, shift, max}` |
| 400 | `GENOTYPING_SET_INVALID` | check | a set's structure or position ([reasons](#post-primerscheck-the-genotyping-block)) | `{set_id, reason, pair_ids, primer?, sequence?, mm_pos?}` |
| 400 | `INVALID_REQUEST` | check | `genotyping` outside gene and region modes | `{field: "genotyping", reason: "mode"}` |
| 404 | `UNKNOWN_VARIANT` | lookup, design | Ensembl does not know the id | `{id, system_name}` |
| 422 | `NO_VARIATION_DATA` | variants list and lookup, design by id | the genome has no variation data | `{system_name}` |
| 422 | `AMBIGUOUS_VARIANT_MAPPING` | lookup, design | the id maps to several sequences of the assembly | `{id, mappings}` |
| 422 | `VARIANT_NOT_ON_ASSEMBLY` | lookup, design | the id maps to no sequence of the assembly | `{id}` |
| 500 | `THERMO_FAILED` | design | `ntthal` failed or printed nothing usable | `{binary}` |
| 503 | `THERMO_UNAVAILABLE` | design | `ntthal` is missing or could not be started | `{binary, retry_after_s}` |
| 503 | `VARIATION_SOURCE_UNAVAILABLE` | variants list and lookup, design by id | an Ensembl call failed | `{retry_after_s, reason}`, where `reason` is `timeout`, `transport`, `http_5xx`, `rate_limited`, `invalid_response`, `breaker_open` or `queue_full` |
| 503 | `FEATURE_DISABLED` | variants list and lookup, design by id | `variation.enabled` is false | `{retry_after_s: 300}` |

### Genotyping warning codes

Warnings are `{code, message, details}`.

**Variants endpoints:** `DUPLICATE_VARIANT_IDS {key, ids}`; `VARIATION_RECORDS_SKIPPED {count, reasons}` (malformed Ensembl
records dropped; an id is never a reason); `REF_MISMATCHES {count}`; `VARIANTS_TRUNCATED {returned, total, limit}`.

**Design, response level**, in the order they appear:

| Code | When | `details` |
| --- | --- | --- |
| `DUPLICATE_VARIANT_IDS` | several Ensembl ids describe the target | `{key, ids}` |
| `EMS_TARGET` | every record of the target is EMS, so natural neighbours warn instead of blocking | `{sources}` |
| `MULTIALLELIC_SITE` | the site has other alternative alleles; assemblies carrying them mismatch both AS primers | `{key, other_alts}` |
| `SHIFTABLE_INDEL` | `shift > 0`: the forward and reverse primers end at different bases | `{shift, forward_position, reverse_position, zone}` |
| `ORIENTATION_BLOCKED` | a non-EMS known variant in the last 5 nt of an AS primer | `{orientation, ids, distances}` |
| `NO_VARIATION_DATA` | the genome has no variation data, or it is switched off: neighbours were not screened | `{system_name}` |
| `NEIGHBOURS_UNAVAILABLE` | Ensembl failed during a manual design: neighbours were not screened | `{reason}` |
| `DENSE_NEIGHBOURS` | more than 2 non-EMS known variants within 30 bp of the variant | `{count, window}` |
| `SUBMISSION_NEIGHBOURS_OMITTED` | indel or multi-allelic neighbours in the submission flanks | `{ids}` |
| template and mask warnings | as in `/primers/design`: the assembly warnings of the template, and `REPEAT_MASK_FAILED`, `NO_REPEAT_MASK`, `BLAST_DEPTH_MASK`, `MOSTLY_REPEAT` with `avoid_repeats` | |
| `VARIANT_IN_REPEAT` | the repeat mask covered bases of the exempt AS window | `{orientation, masked_bases}` |
| `PRIMER3_WARNING` | Primer3 printed a warning | `{orientation, level}` |
| `ORIENTATION_SKIPPED` | too close to the region end, or an `N` in the AS window | `{orientation, reason}` |
| `ORIENTATION_NO_SETS` | the ladder was exhausted | `{orientation, levels_tried}` |
| `RELAXED_CONSTRAINTS` | the orientation's sets needed level 1 or 2 | `{orientation, level, changes}` |
| `NO_SETS` | no set in either orientation (`sets: []`; read `orientations.*.attempts[].explain`) | `{}` |
| `DESIGN_BUDGET_EXHAUSTED` | the process budget left pairs unscored | `{orientation, primer3_runs, thermo_calls, max_primer3_runs, max_thermo_calls, not_scored, sets_returned}` |

**Design, set level** (`sets[].issues[]` with a `severity`; the warn and high ones are repeated in `sets[].warnings`):

| Code | When | Severity | `details` |
| --- | --- | --- | --- |
| `ALT_PRIMER_SUBOPTIMAL` | an ordered derived oligo (the ALT primer; with a deliberate mismatch, both AS primers) has Primer3 problems but clears the floors | warn | `{oligo, problems}` |
| `AS_TM_IMBALANCE` | the AS primers' Tm differ by more than 1.0 °C | warn; high above 2.0 °C | `{diff}` |
| `COMMON_TM_OUT_OF_RANGE` | `common_minus_as` is outside −1.0…+3.0 °C | warn | `{value}` |
| `TAILED_STRUCTURE` | a tailed AS primer's hairpin or self-dimer Tm ≥ 47 °C, or a tailed cross-dimer Tm ≥ 47 °C | warn; high for a hairpin ≥ 55 °C | a primer: `{oligo, metric, value}`; a cross-dimer: `{pair, metric, value}` |
| `MISMATCH_NOT_APPLICABLE` | the base at −k differs between the AS primers, so no mismatch was applied | warn | `{position}` |
| `MISMATCH_STRUCTURE` | a deliberate-mismatch primer's hairpin or self-dimer Tm ≥ 47 °C | warn; high ≥ 55 °C | `{oligo, metric, value}` |
| `NEIGHBOUR_IN_PRIMER` | non-EMS known variants in a primer outside its last 5 nt: one issue for the AS pair (`role: "as_ref"`, both primers' variants) and one for the common primer | warn | `{role, ids, distances}` |
| `NEIGHBOUR_AT_3P` | a non-EMS known variant in the last 5 nt that was not allowed to block (EMS target, or `neighbour_policy: "ignore"`) | high | `{role, ids, distances}` |
| `WEAK_DISCRIMINATION` | an AS primer is predicted to amplify the other allele (`likely`) | high | `{primer, mm_pos}` |
| `SHIFT_TRACT_DISCRIMINATION` | an AS primer's 3′ base lies inside the indel's shift tract | warn | `{primer, shift}` |
| `WEAK_TERMINAL_CLASS` | both AS 3′ mismatches are in Little's weak class (A/G and C/T sites) | info: in `issues` only, not scored | `{}` |

In `sets[].warnings`, `TAILED_STRUCTURE` and `MISMATCH_STRUCTURE` details also carry the `severity`.

**Check, results level** (`results.warnings` of a genotyping job): `REFERENCE_CONTROL_FAILED`, `WEAK_OFF_TARGETS`,
`GENOTYPE_FALLBACK_FAILED`, `GENOTYPE_FALLBACK_BUDGET` and `GENOTYPE_FAILED`, described under
[Results warnings](#results-warnings).

### Results: `results.genotyping`

Present in the `results` of `GET /primers/check/{job_id}` only for jobs whose request carried a `genotyping` block. Other
jobs have no `genotyping` key at all (not `null`), so their results are exactly what they were before genotyping existed.
For the reference genome, which serves as a control, and for each pan-genome assembly, it reports the allele found at the
variant. For each set, it reports the predicted signal of the REF, ALT and common primers, and so the genotype an assay
would show. The definitions are `PrimerCheckGenotypingResults` and its members in `api/swagger/swagger.yaml` (spec §2.12).
Where the spec differs, this subsection describes the implementation.

**`specificity` and `pangenome` are unchanged.** A genotyping job computes them exactly as any other check does, and they
still list every off-target product, including the ones that do not change an allele prediction (see
[off-locus products](#off-locus-products-and-weak_off_targets)). On their own they cannot answer allele questions. In the
[example](#example-rs871475760-sets-s1-and-s2-over-11-assemblies), both S1 pairs amplify in all 11 assemblies
(`pangenome.pairs[].summary.amplifies: 11`), because a primer with a single 3′-terminal mismatch still gives a
`likely_weak` product. Which allele an assembly carries, and which primer would give signal there, is answered **only** by
`results.genotyping`.

| Member | Contents |
| --- | --- |
| `algorithm_version` | `"g1"`, the version of the allele caller and the prediction rules. Genotyping jobs hash with algorithm `"2+g1"`, so a new caller version changes the ids of genotyping jobs only |
| `variant` | The left-aligned variant `{key, region, position, ref, alt, shift, zone}`, plus: `flank` (K, the flank on each side of the core: `max(check.genotype_flank_min, zone length + the longer allele's length)`); `core {ref, alt}` (the reference span `[position − 1, position + len(ref) + shift]` read on each haplotype; an assembly is `ref` or `alt` only when its core is exactly one of them); and `haplotypes {ref, alt}` (the core plus K bases on each side) |
| `summary` | Allele counts over the **pan-genome genomes only** (the reference is not counted): `ref`, `alt`, `other`, `ambiguous`, `missing` and `unavailable`, which add up to `genomes_total`, the same total as `pangenome.pairs[].summary.genomes_total` |
| `genomes[]` | One entry per genome: the reference first (`is_reference: true`), then the finished pan-genome genomes in the order of the job's normalized `request.genomes`. That order is sorted by `system_name`, as in `pangenome.pairs[].genomes` |
| `sets[]` | One per request set, in request order: `id`, `ref_pair`, `alt_pair`, `orientation` (derived at submit), `deliberate_mismatch_positions` (the mismatches accepted at submit, `[]` for KASP), `specificity`, `control`, `reference` (the prediction on the reference), `summary` and `genomes[]` (one prediction per pan-genome genome, in the order of `genomes[]` without the reference) |

#### Allele calls per genome

Each `genomes[]` entry is `{system_name, display_name, is_reference, allele, observed, source, copies, orthologous_copies,
paralog_copies, reason}`.

| `allele` | Meaning | `observed` | `source` | `reason` |
| --- | --- | --- | --- | --- |
| `ref` | Every orthologous copy that covers the variant reads exactly `variant.core.ref` | the core | `amplicon` or `megablast` | `null` |
| `alt` | Every such copy reads exactly `variant.core.alt` | the core | `amplicon` or `megablast` | `null` |
| `other` | Every such copy reads a core that is neither: a third allele, such as another repeat length or a third base | the core, or `null` when the copies read different ones | `amplicon` or `megablast` | `null` |
| `ambiguous` | The copies disagree, e.g. one `ref` copy and one `alt` copy | `null` | `amplicon` or `megablast` | `null` |
| `missing` | No copy covers the variant | `null` | `null` | a missing reason (below) |
| `unavailable` | The assembly could not be searched, or the caller failed. Every prediction is `unknown` | `null` | `null` | `db_unavailable` (no usable BLAST DB), `blast_error` (its BLAST failed) or `call_failed` (the caller threw; warning `GENOTYPE_FAILED`) |

How a genome is called (spec §5.6, as implemented):

1. **Anchors.** The caller uses every product of the set's two pairs with orientation `LR` or `RL`: first the amplifying
   ones, then the `unlikely` ones, because a blocked allele-specific primer's product is still a good place to read the
   allele. Products are de-duplicated by position, at most `check.genotype_max_anchors` (50) per set.
   - `LL`/`RR` products are never anchors.
   - Approximate products (`approx: true`) are never anchors either: their ends were extrapolated.
2. **Alignment.** Each anchor's genome window is its product ± 2 × `check.genotype_amplicon_pad` (± 100 bp),
   reverse-complemented for `RL`. The set's reference segment, `expected` ± `check.genotype_amplicon_pad` (50 bp), is
   aligned to that window semi-globally:
   - the whole segment is aligned, and the genome ends are free;
   - every edit costs 1;
   - an indel stays one gap run.
3. **The call.** The two end bases of `variant.core` are mapped through the alignment. An end that falls on a gap moves
   outward to the nearest aligned base.
   - The genome bases between the two ends are the copy's `observed` core.
   - The copy's `call` is `ref` or `alt` only on an exact match, and otherwise `other`.
   - An end outside the aligned span gives `call: "missing"`: the copy does not cover the variant.
   - Edits outside the core are counted in `flank_edits` and never change the call.
4. **Copies.** Anchors on the same region and strand whose aligned variant coordinates are within the zone length of each
   other form one copy. This merges the products of both pairs, and of every set, at that locus.
   - `start`/`end` is the envelope of the products, and `anchors` is their number.
   - `identity`, `gap_compressed_identity`, `aligned_length`, `observed`, `flank_edits` and `call` come from the copy's
     alignment with the most columns. `variant_position` is the genome coordinate aligned to `variant.position`.
   - A product of the wrong size (step 5) that lies at the locus of a right-sized copy is left out: it is another product
     of one of the primers there.
5. **Orthologs and paralogs.** A copy is orthologous in either of these cases:
   - it is annotated as an ortholog (`ortholog: true`, gene mode);
   - its gap-compressed identity (each gap run counts as one edit and one column) is at least
     `check.genotype_ortholog_min_identity` (95), and every one of its products is within
     ± `check.genotype_ortholog_size_tolerance` (20 %) of its set's reference product.

   Gap compression keeps a true ortholog that carries one long indel. A 35 bp deletion gives 94.9 % raw identity but
   99.85 % gap-compressed. The pi536008 copy on scaffold2100, whose alignment runs off the scaffold end, has 88.02 % raw
   identity and 98.83 % gap-compressed.
   - `copies` lists the orthologous copies (at most `check.genotype_max_copies`, 10), and `orthologous_copies` counts them.
   - `paralog_copies` counts the other copies that still align at ≥ 80 % gap-compressed identity, such as pi180348's copy
     at 1:72.9 Mb (90 %). Less similar products are unrelated off-targets and are not counted: a real check finds 11–16
     per assembly at rs871475760, all at 66–77 %. Paralogs never vote.
6. **The genome's allele.** Copies whose call is `missing` do not vote. When the orthologous copies that cover the variant
   agree, their call is the allele (`source: "amplicon"`); when they disagree, the allele is `ambiguous`.
7. **Megablast fallback.** When no orthologous copy covers the variant, the worker runs one megablast of the reference
   `[zone.start − 200 − K, zone.end + 200 + K]` against the assembly's BLAST database.
   - Limits: single thread, no retry, a `check.genotype_megablast_timeout_ms` (20 s) timeout, and at most
     `check.genotype_max_megablast` (30) fallbacks per job.
   - An HSP is kept when its bitscore is ≥ `check.genotype_megablast_min_bitscore_frac` (0.9) × the best, its identity is
     ≥ `check.genotype_megablast_min_identity` (95 %), its query cover is ≥ `check.genotype_megablast_min_query_cover`
     (0.8), and it covers the core.
   - Kept HSPs are read with the same exact-core rule and merged into copies (`source: "megablast"`, `anchors: 0`).
     Megablast copies count as orthologous.
   - The query-cover filter misses a locus at a scaffold end. pi536008's scaffold2100 ends 176 bp past the variant, so its
     HSP covers 272 of the 431 query bases; the fallback alone would call that genome `missing`. Its amplicons call it
     `alt`.
8. **Missing reasons.** The megablast failure and budget reasons take precedence over the other two.

   | `reason` | When |
   | --- | --- |
   | `no_orthologous_copy` | No orthologous copy, and megablast kept no HSP |
   | `variant_not_covered` | An orthologous copy exists, or megablast kept an HSP, but none covers the core |
   | `fallback_failed` | Megablast failed or timed out (results warning `GENOTYPE_FALLBACK_FAILED`) |
   | `fallback_budget` | The job had already run `check.genotype_max_megablast` fallbacks (results warning `GENOTYPE_FALLBACK_BUDGET`) |

#### Amplification prediction

Each `sets[].reference` and `sets[].genomes[]` item is `{system_name, ref_primer, alt_primer, common_primer, predicted,
strength, agrees, reasons, off_locus_products}`.

**On-locus and off-locus products.** An `LR`/`RL` product of either pair is on-locus when it overlaps one of the genome's
orthologous `copies`; every other `LR`/`RL` product is off-locus. `LL`/`RR` products are not read.

**Primer calls** are `{status, likelihood, mm_pos, residual_mm_pos}`.
- `mm_pos` lists the primer's mismatch positions in the product, as distances from the 3′ end (1 is the 3′ base).
- `residual_mm_pos` is `mm_pos` without the primer's own declared deliberate mismatch, and the prediction reads it.

| `status` | Allele-specific primer (`ref_primer`, `alt_primer`): its own pair's on-locus products | Common primer (`common_primer`): its own site in the on-locus products of either pair |
| --- | --- | --- |
| `no_product` | The pair has no on-locus product; `likelihood`, `mm_pos` and `residual_mm_pos` are `null` | Neither pair has one |
| `blocked` | The product is `unlikely` | The site alone is ignored, over `max_amplifying_mismatches`, or 3′-blocked |
| `unknown` | `mm_pos` is `null` (an approximate alignment) | Same rule |
| `uncertain` | A residual mismatch at position 1, while the primer's 3′ base lies in the indel's shift tract | Not used |
| `terminal_mismatch` | A residual mismatch at position 1 | A mismatch at position 1 |
| `weak` | A residual mismatch at 2 or 3, and none at 1 | Same rule |
| `match` | Anything else: mismatches further from the 3′ end do not count | Same rule |

A primer's status is the best over its products, in the order `match`, `weak`, `terminal_mismatch`, `uncertain`,
`unknown`, `blocked`, `no_product`. Its `likelihood` and positions come from the product that gave that status; on ties
the better likelihood wins.

The common primer's status never comes from a pair's likelihood, and no deliberate mismatch is ever subtracted from it.
An AS-PCR set's deliberate mismatch therefore never makes the common primer look mismatched. In the example, s3691's S2
common primer has a mismatch at position 13 and is still `match`.

**`predicted`**: the first row that applies.

| # | Situation | `predicted` | `strength` | `reasons` |
| --- | --- | --- | --- | --- |
| 1 | The genome is `unavailable` | `unknown` | `null` | none; every primer call is `unknown` |
| 2 | Neither pair has an on-locus product | `none` | `null` | `no_orthologous_copy` when the allele is `missing` |
| 3 | The common primer is `terminal_mismatch` or `blocked`: it is in both pairs, so neither allele can give signal | `no_call` | `null` | `common_primer_3p_mismatch` |
| 4 | A primer is `unknown` or `uncertain` | `unknown` | `null` | `approx_alignment`, `shift_tract_uncertain` |
| 5 | The REF primer amplifies (`match` or `weak`) and the ALT primer does not | `ref` | `weak` when the amplifying primer or the common primer is `weak`, else `normal` | `common_primer_weak` when the common primer is `weak` |
| 6 | The ALT primer amplifies and the REF primer does not | `alt` | as in row 5 | as in row 5 |
| 7 | Both amplify | `both` | as in row 5 | as in row 5, plus `third_allele` when the allele is `other` |
| 8 | Neither amplifies | `none` | `null` | |

**`agrees`** compares the prediction with the genome's allele:
- `true` for `ref` predicted on a `ref` genome, `alt` on `alt`, and `none` on an `other` or `missing` genome;
- `null` when the allele is `ambiguous` or `unavailable`, the prediction is `unknown` or `no_call`, or a primer is
  `uncertain`;
- `false` otherwise.

#### Off-locus products and `WEAK_OFF_TARGETS`

An off-locus product can still give signal, for example a paralog elsewhere in the genome that one allele-specific primer
and the common primer amplify. It **counts** when all of these hold:
1. it amplifies (`likely` or `likely_weak`);
2. its allele-specific primer has a known alignment with no residual mismatch at position 1;
3. its common-primer site alone is `match` or `weak`;
4. **each primer has at most `check.genotype_offlocus_max_mismatches` mismatches** (default 2). The allele-specific
   primer's declared deliberate mismatch is not counted.

Each counted product adds 1 to `off_locus_products`, with the reason `ref_signal_off_locus` (a REF-pair product) or
`alt_signal_off_locus` (an ALT-pair product). It changes the prediction only by adding the signal that is missing:
- a counted REF-pair product turns `alt` into `both`, and an ALT-pair product turns `ref` into `both`;
- on a genome predicted `none` or `no_call`, the off-locus products alone decide (`ref`, `alt` or `both`), and `strength`
  stays `null`;
- a product of the pair that already amplifies (a same-allele paralog) leaves `predicted` unchanged;
- `unknown` never changes.

The reference control uses the same rule.

**`WEAK_OFF_TARGETS`.** A product that meets conditions 1–3 but has more mismatches in a primer is left out of
`off_locus_products`. It is reported instead in one job-level results warning, `WEAK_OFF_TARGETS {count, max_mismatches,
examples}`:
- `count` is the number of such products over the whole job, in every genome, set and pair.
- `max_mismatches` is the threshold.
- `examples` holds at most 5 of them, as `{system_name, set_id, pair_id, allele, region, start, end, size, orientation,
  left_mm, right_mm}`. They are ordered by genome (the reference first, then `request.genomes`), then set, REF before ALT,
  region and start. `left_mm` and `right_mm` count mismatches the way the rule does.
- The message names the genomes that have such products: the first 10, then "and N more".

These products **still appear** in `specificity.pairs[].off_targets` and in the pan-genome results. The threshold only
decides which of them change an allele prediction.

**Why 2.** The user chose the threshold on 2026-09-15 from the real rs871475760 data; see the M8 deviations at the end of
`docs/genotyping_design_spec.md`.
- **Products it leaves out.** S1_REF, and an AS-PCR set on the same common primer, each have three off-target products in
  the reference. Each product has 2–3 mismatches in each primer. The same products recur in other assemblies: S1_REF's
  9:14931544–14931612 product, for example, has a counterpart in 10 of the 11 (on chromosome 10 in pi329250). When they
  were counted, they turned all five ALT assemblies into `both`.
- **The product it keeps.** pi180348's paralog at 1:72907129–72907210 has 1 mismatch in S2's REF primer and 2 in the
  common primer. It still counts, so S2 predicts `both` on pi180348.
- **Size does not help.** Every amplifying off-target measured is 65–82 bp, as large as the on-target products (65 and
  92 bp), so a size limit cannot separate the two cases.

#### Reference control, pair verdicts and summaries

- **`sets[].reference`** is the set's prediction on the reference genome, computed like any other genome.
- **`sets[].control`** is `{status, allele, reasons}`:
  - `pass`: the reference is called `ref`, and the set predicts `ref` on it;
  - `warn`: the same, but the prediction is `weak` (reason `weak`), or it has counted off-locus products (reason
    `off_locus_products`);
  - `fail`: the reference is not `ref` (`allele_not_ref`), or the set does not predict `ref` on it (`prediction_not_ref`).
    An ALT-pair off-locus product on the reference fails the control this way, through `both`.

  `reasons` lists the control's own codes, followed by the reasons of the reference prediction. Any `fail` adds the
  results warning `REFERENCE_CONTROL_FAILED {allele, sets}`, where `sets` are the failing set ids. Treat that set's
  predictions on the other genomes as unreliable.
- **`sets[].specificity`** is `{ref_pair, alt_pair, off_target_count}`. It restates the reference verdicts of
  `results.specificity` for the set's two pairs.
  - `consistent_with_allele` is `true` when the verdict is what allele specificity leads you to expect: `specific` or
    `off_targets` for an allele-specific primer with signal on the reference, and `on_target_missing` for one that is
    `blocked` or has no product.
  - `off_target_count` is the number of distinct reference off-target products of the two pairs.
  - In the example, S1_REF's verdict is `off_targets` with 3 products, which is consistent with the allele. The control
    still passes, because none of those products counts under the off-locus rule.
- **Summaries.**
  - `genotyping.summary` counts alleles over the pan-genome genomes.
  - Each `sets[].summary` counts that set's pan-genome predictions:
    `predicted_ref + predicted_alt + both + none + no_call + unknown = genomes_total`, and
    `agree + disagree + not_comparable = genomes_total` (`agrees` true, false or null). `weak` counts the predictions
    with `strength: "weak"`.
  - The reference is never counted.
- **Specificity-only jobs** (`checks: ["specificity"]`) carry the reference entry alone, with its controls and
  predictions, and zero summaries.
- **Partial results.** `results.genotyping` first appears in the partial results written after the reference stage. These
  hold the reference entry and every set's `reference`, `control` and `specificity`.
  - Each later flush adds the pan-genome genomes finished so far, in `request.genomes` order, and recomputes the summaries.
  - `genomes` therefore grows along with `pangenome.pairs[].genomes`. The example job showed 1, 2, 9 and 10 genome entries
    while it ran, and 12 when it was done.

#### Results warnings

| Code | When | `details` |
| --- | --- | --- |
| `REFERENCE_CONTROL_FAILED` | A set's reference control is `fail` | `{allele, sets}` |
| `WEAK_OFF_TARGETS` | Off-locus products were left out of the predictions by `check.genotype_offlocus_max_mismatches` | `{count, max_mismatches, examples}` (above) |
| `GENOTYPE_FALLBACK_FAILED` | A megablast fallback failed or timed out; those genomes are `missing` (`fallback_failed`) | none; the message names the genomes |
| `GENOTYPE_FALLBACK_BUDGET` | The job's megablast budget was used up; later genomes are `missing` (`fallback_budget`) | none; the message names the genomes |
| `GENOTYPE_FAILED` | The allele caller threw for a genome; that genome is `unavailable` (`call_failed`), and the job goes on | none; the message names the genomes |

#### Caller configuration

The allele caller reads these `primers.check` keys; their defaults are in `config.js` and `config/default.yaml`. Config is
not part of the job id, so a finished job keeps the results it was computed with until it expires.

| Key | Default | What it does |
| --- | --- | --- |
| `check.genotype_cpu_s_per_genome` | 0.2 | CPU-s added to the submit-time estimate per genome task (the reference plus each pan-genome genome) |
| `check.genotype_flank_min` | 15 | The smallest K, the flank of `variant.haplotypes` |
| `check.genotype_amplicon_pad` | 50 | The reference segment aligned for a set is `expected` ± this many bases; an anchor's genome window is its product ± twice this |
| `check.genotype_ortholog_min_identity` | 95 | The lowest gap-compressed identity (%) of an orthologous copy |
| `check.genotype_ortholog_size_tolerance` | 0.2 | How far (as a fraction) an orthologous copy's products may differ in size from the set's reference product |
| `check.genotype_max_copies` | 10 | The most copies listed per genome; `orthologous_copies` still counts them all |
| `check.genotype_max_anchors` | 50 | The most products per set used as anchors in one genome |
| `check.genotype_max_megablast` | 30 | Megablast fallbacks per job; past it, a genome is `missing` with `fallback_budget` |
| `check.genotype_megablast_timeout_ms` | 20000 | Timeout of one megablast (no retry) |
| `check.genotype_megablast_min_identity` | 95 | The lowest percent identity of a kept HSP |
| `check.genotype_megablast_min_query_cover` | 0.8 | The smallest fraction of the query a kept HSP covers |
| `check.genotype_megablast_min_bitscore_frac` | 0.9 | A kept HSP's bitscore is at least this fraction of the best one |
| `check.genotype_offlocus_max_mismatches` | 2 | An off-locus product changes a prediction only when each primer has at most this many mismatches; the others go to `WEAK_OFF_TARGETS` |

These values are fixed in code: the 80 % identity floor of `paralog_copies`, the 200 bp megablast flank, the ± 2 pads of
an anchor window, and the 5 `WEAK_OFF_TARGETS` examples.

#### Measured cost

The estimate adds `check.genotype_cpu_s_per_genome` (0.2 CPU-s) per genome task
([above](#post-primerscheck-the-genotyping-block)). Measured on squam on 2026-09-15:
- **The allele caller** used 0.03–0.10 CPU-s per genome (mean 0.07) and one megablast fallback 0.09–0.13 CPU-s. This was
  for sets S1, S2 and an AS-PCR set, over the reference and the 11 research assemblies, with each call replayed in process.
  0.2 CPU-s therefore covers a call plus a fallback.
- **The whole example job** was estimated at 273 CPU-s (270 without `genotyping`). The worker and its BLAST
  processes used 176.8 CPU-s, read from `/proc`, over 36.5 s. The estimate therefore runs about **1.5×** the measured
  CPU. It is conservative, as the [cost guard](#cost-guard) intends.

#### Example: rs871475760, sets S1 and S2 over 11 assemblies

**The job.** The `check.request` of the rs871475760 KASP design ([above](#post-primerscheck-the-genotyping-block)), with
`genomes` set to the 11 research assemblies. It was submitted on 2026-09-15 to the dev API and its worker on
127.0.0.1:50112 (worktree HEAD 228a76b). The 202 response is shown with the request.

**What the excerpt leaves out** (each omission is marked `…`):
- from `results`: `engine`, `params`, `reference`, `sensitivity_note`, `primers` and `transcriptome`;
- from `request`: `params`, `pairs` and `genotyping`;
- in the lists: 10 of the 12 `genotyping.genomes` entries and most per-set genome entries.

<!-- example: request GET /primers/check/{job_id} capture=capture-check-genotyping-result.json#/request/path -->
```http
GET /sorghum_v11/primers/check/8e9160d598f602137d94efa9fc409264
```

<!-- example: response GET /primers/check/{job_id} 200 capture=capture-check-genotyping-result.json#/response -->
```json
{
  "job_id": "8e9160d598f602137d94efa9fc409264", "status": "done", "kind": "pangenome", "partial": false,
  "queue_position": null, "progress": { "done": 12, "total": 12, "stage": "done", "running": [] },
  "created_at": "2026-09-15T22:54:34.441Z", "started_at": "2026-09-15T22:54:34.640Z",
  "finished_at": "2026-09-15T22:55:11.207Z", "attempts": 1,
  "request": {
    "system_name": "sorghum_bicolor", "mode": "region", "checks": ["pangenome", "specificity"],
    "genomes": [
      "sorghum_austrcf317961", "sorghum_bicolorv5", "sorghum_pi180348", "sorghum_pi276837", "sorghum_pi329250",
      "sorghum_pi510757", "sorghum_pi565121", "sorghum_pi655972", "sorghum_pi656027", "sorghum_rio", "sorghum_s3691"
    ],
    "…": "…"
  },
  "warnings": [], "estimate": { "cpu_s": 273 },
  "results": {
    "specificity": {
      "target": "genome",
      "pairs": [
        { "id": "S1_REF", "verdict": "off_targets", "off_target_count": 3, "…": "…" },
        { "id": "S1_ALT", "verdict": "specific", "off_target_count": 0, "…": "…" },
        { "id": "S2_REF", "verdict": "specific", "off_target_count": 0, "…": "…" },
        { "id": "S2_ALT", "verdict": "specific", "off_target_count": 0, "…": "…" }
      ]
    },
    "pangenome": {
      "target": "genome",
      "pairs": [
        {
          "id": "S1_REF",
          "summary": {
            "genomes_total": 11, "single_perfect": 0, "single_mismatch": 0, "multiple": 11, "no_amplicon": 0,
            "db_unavailable": 0, "error": 0, "amplifies": 11, "truncated": 0
          },
          "…": "…"
        },
        {
          "id": "S1_ALT",
          "summary": {
            "genomes_total": 11, "single_perfect": 1, "single_mismatch": 6, "multiple": 4, "no_amplicon": 0,
            "db_unavailable": 0, "error": 0, "amplifies": 11, "truncated": 0
          },
          "…": "…"
        },
        "…"
      ]
    },
    "genotyping": {
      "algorithm_version": "g1",
      "variant": {
        "key": "1:11109:C:A", "region": "1", "position": 11109, "ref": "C", "alt": "A", "shift": 0,
        "zone": { "start": 11109, "end": 11109 }, "flank": 15, "core": { "ref": "TCT", "alt": "TAT" },
        "haplotypes": { "ref": "ATAGTCATACTCTATTCTGAATTTCTCGCTAGT", "alt": "ATAGTCATACTCTATTATGAATTTCTCGCTAGT" }
      },
      "summary": {
        "genomes_total": 11, "ref": 6, "alt": 5, "other": 0, "ambiguous": 0, "missing": 0, "unavailable": 0
      },
      "genomes": [
        {
          "system_name": "sorghum_bicolor", "display_name": "Sb bicolor BTx623 v3", "is_reference": true,
          "allele": "ref", "observed": "TCT", "source": "amplicon",
          "copies": [
            {
              "region": "1", "start": 11068, "end": 11172, "strand": 1, "variant_position": 11109, "identity": 100,
              "gap_compressed_identity": 100, "aligned_length": 192, "observed": "TCT", "flank_edits": 0,
              "call": "ref", "anchors": 2, "ortholog": null, "source": "amplicon"
            }
          ],
          "orthologous_copies": 1, "paralog_copies": 0, "reason": null
        },
        {
          "system_name": "sorghum_pi180348", "display_name": "Sb bicolor PI180348 Juar (IS 12876)",
          "is_reference": false, "allele": "alt", "observed": "TAT", "source": "amplicon",
          "copies": [
            {
              "region": "1", "start": 15028, "end": 15132, "strand": -1, "variant_position": 15091, "identity": 98.44,
              "gap_compressed_identity": 98.44, "aligned_length": 192, "observed": "TAT", "flank_edits": 2,
              "call": "alt", "anchors": 2, "ortholog": null, "source": "amplicon"
            },
            {
              "region": "1", "start": 38342, "end": 38446, "strand": 1, "variant_position": 38383, "identity": 98.45,
              "gap_compressed_identity": 98.45, "aligned_length": 193, "observed": "TAT", "flank_edits": 2,
              "call": "alt", "anchors": 2, "ortholog": null, "source": "amplicon"
            }
          ],
          "orthologous_copies": 2, "paralog_copies": 1, "reason": null
        },
        "…"
      ],
      "sets": [
        {
          "id": "S1", "ref_pair": "S1_REF", "alt_pair": "S1_ALT", "orientation": "reverse",
          "deliberate_mismatch_positions": [],
          "specificity": {
            "ref_pair": { "verdict": "off_targets", "consistent_with_allele": true },
            "alt_pair": { "verdict": "specific", "consistent_with_allele": true }, "off_target_count": 3
          },
          "control": { "status": "pass", "allele": "ref", "reasons": [] },
          "reference": {
            "system_name": "sorghum_bicolor",
            "ref_primer": { "status": "match", "likelihood": "likely", "mm_pos": [], "residual_mm_pos": [] },
            "alt_primer": {
              "status": "terminal_mismatch", "likelihood": "likely_weak", "mm_pos": [1], "residual_mm_pos": [1]
            },
            "common_primer": { "status": "match", "likelihood": "likely", "mm_pos": [], "residual_mm_pos": [] },
            "predicted": "ref", "strength": "normal", "agrees": true, "reasons": [], "off_locus_products": 0
          },
          "summary": {
            "genomes_total": 11, "predicted_ref": 6, "predicted_alt": 5, "both": 0, "none": 0, "no_call": 0,
            "unknown": 0, "weak": 0, "agree": 11, "disagree": 0, "not_comparable": 0
          },
          "genomes": [
            {
              "system_name": "sorghum_pi180348",
              "ref_primer": {
                "status": "terminal_mismatch", "likelihood": "likely_weak", "mm_pos": [1], "residual_mm_pos": [1]
              },
              "alt_primer": { "status": "match", "likelihood": "likely", "mm_pos": [], "residual_mm_pos": [] },
              "common_primer": { "status": "match", "likelihood": "likely", "mm_pos": [], "residual_mm_pos": [] },
              "predicted": "alt", "strength": "normal", "agrees": true, "reasons": [], "off_locus_products": 0
            },
            "…"
          ]
        },
        {
          "id": "S2", "ref_pair": "S2_REF", "alt_pair": "S2_ALT", "orientation": "forward",
          "deliberate_mismatch_positions": [],
          "specificity": {
            "ref_pair": { "verdict": "specific", "consistent_with_allele": true },
            "alt_pair": { "verdict": "specific", "consistent_with_allele": true }, "off_target_count": 0
          },
          "control": { "status": "pass", "allele": "ref", "reasons": [] },
          "reference": {
            "system_name": "sorghum_bicolor",
            "ref_primer": { "status": "match", "likelihood": "likely", "mm_pos": [], "residual_mm_pos": [] },
            "alt_primer": {
              "status": "terminal_mismatch", "likelihood": "likely_weak", "mm_pos": [1], "residual_mm_pos": [1]
            },
            "common_primer": { "status": "match", "likelihood": "likely", "mm_pos": [], "residual_mm_pos": [] },
            "predicted": "ref", "strength": "normal", "agrees": true, "reasons": [], "off_locus_products": 0
          },
          "summary": {
            "genomes_total": 11, "predicted_ref": 6, "predicted_alt": 4, "both": 1, "none": 0, "no_call": 0,
            "unknown": 0, "weak": 0, "agree": 10, "disagree": 1, "not_comparable": 0
          },
          "genomes": [
            {
              "system_name": "sorghum_pi180348",
              "ref_primer": {
                "status": "terminal_mismatch", "likelihood": "likely_weak", "mm_pos": [1], "residual_mm_pos": [1]
              },
              "alt_primer": { "status": "match", "likelihood": "likely", "mm_pos": [], "residual_mm_pos": [] },
              "common_primer": { "status": "match", "likelihood": "likely", "mm_pos": [], "residual_mm_pos": [] },
              "predicted": "both", "strength": "normal", "agrees": false, "reasons": ["ref_signal_off_locus"],
              "off_locus_products": 1
            },
            {
              "system_name": "sorghum_s3691",
              "ref_primer": { "status": "match", "likelihood": "likely", "mm_pos": [], "residual_mm_pos": [] },
              "alt_primer": {
                "status": "terminal_mismatch", "likelihood": "likely_weak", "mm_pos": [1], "residual_mm_pos": [1]
              },
              "common_primer": { "status": "match", "likelihood": "likely", "mm_pos": [13], "residual_mm_pos": [13] },
              "predicted": "ref", "strength": "normal", "agrees": true, "reasons": [], "off_locus_products": 0
            },
            "…"
          ]
        }
      ]
    },
    "warnings": [
      {
        "code": "WEAK_OFF_TARGETS",
        "message": "off-target products with more than 2 mismatches in a primer do not change the allele predictions; the specificity and pan-genome results still list them (sorghum_bicolor, sorghum_s3691, …, sorghum_pi510757 and 2 more)",
        "details": {
          "count": 28, "max_mismatches": 2,
          "examples": [
            {
              "system_name": "sorghum_bicolor", "set_id": "S1", "pair_id": "S1_REF", "allele": "ref", "region": "1",
              "start": 372590, "end": 372654, "size": 65, "orientation": "RL", "left_mm": 3, "right_mm": 3
            },
            {
              "system_name": "sorghum_bicolor", "set_id": "S1", "pair_id": "S1_REF", "allele": "ref", "region": "4",
              "start": 46360759, "end": 46360823, "size": 65, "orientation": "RL", "left_mm": 3, "right_mm": 3
            },
            "…"
          ]
        }
      }
    ],
    "timings_ms": {
      "reference": 10764, "sorghum_s3691": 11324, "sorghum_bicolorv5": 12905, "sorghum_rio": 12979,
      "sorghum_austrcf317961": 13010, "sorghum_pi276837": 14132, "sorghum_pi565121": 14164, "sorghum_pi656027": 14181,
      "sorghum_pi180348": 14195, "sorghum_pi510757": 12365, "sorghum_pi329250": 12366, "sorghum_pi655972": 12639,
      "total": 36542
    },
    "…": "…"
  },
  "error": null
}
```

- **Alleles.** 6 assemblies carry C (`ref`) and 5 carry A (`alt`), all called from their amplicons.
  - pi180348 has two orthologous copies, 1:15028–15132 on the minus strand and 1:38342–38446 on the plus strand (an inverted
    duplication), both `alt`. It also has one paralog, the 1:72.9 Mb copy (`paralog_copies: 1`).
  - pi276837, pi656027 and rio also have two copies each.
- **Pan-genome versus genotyping.** `pangenome` reports that both S1 pairs amplify in all 11 assemblies. `genotyping` says
  that S1's REF primer gives signal in the 6 `ref` assemblies and its ALT primer in the 5 `alt` ones.
- **S1 agrees on 11 of 11.** Its REF pair has three reference off-targets (`specificity`: `off_targets`, 3), each with
  3 mismatches in both primers.
  - They change no prediction, and the control passes.
  - `WEAK_OFF_TARGETS` reports 28 such products over the 12 genomes, all of them from S1_REF.
- **S2 agrees on 10 of 11.** On pi180348, S2's REF pair also amplifies the 1:72.9 Mb paralog, which carries C, with 1 and
  2 mismatches. The prediction there is `both`, with `reasons: ["ref_signal_off_locus"]`, `off_locus_products: 1` and
  `agrees: false`. For this panel, S1 is the better set.
- **Spec §2.13 against this run.**
  - In the spec, S1's `specificity` (`specific`, no off-targets) and pi180348's `paralog_copies` (0) were illustrative;
    the real values are above.
  - The copies, identities, observed cores, calls and primer calls of §2.13 reproduce exactly.
  - Only the S2 prediction on pi180348 differs from §2.13.

#### Two more recorded states: `ambiguous`, and a job in flight

Both come from one full-panel job on the deletion `1:11282:CA:C` (`rs5413864115`), all 119 other assemblies, run on the dev
API and its worker. The recordings are `capture-check-genotyping-ambiguous.json` and `capture-check-genotyping-running.json`.

`ambiguous` means the genome's orthologous copies disagree with each other, so no single allele can be reported. Here
`sorghum_pi544069ph352` has two copies of equal quality, 99.56 % each, reading different cores: `TCAAAG` (`ref`) and
`TCAAG` (`alt`). The genome-level `observed` is therefore `null`, and both sets predict `both`. Note that this excerpt's
`summary` counts all 119 assemblies while the recorded file keeps 17 genome entries, so the arrays deliberately do not sum
to the summary.

<!-- example: response GET /primers/check/{job_id} 200 capture=capture-check-genotyping-ambiguous.json#/response -->
```json
{
  "job_id": "afa2448f22719708b2d5dda511b3d08b", "status": "done", "kind": "pangenome", "partial": false,
  "progress": { "done": 120, "total": 120, "stage": "done", "running": [] },
  "results": {
    "genotyping": {
      "algorithm_version": "g1",
      "summary": {
        "genomes_total": 119, "ref": 111, "alt": 5, "other": 0, "ambiguous": 1, "missing": 2, "unavailable": 0
      },
      "genomes": [
        {
          "system_name": "sorghum_pi544069ph352", "display_name": "Sb bicolor PI544069 PH352", "is_reference": false,
          "allele": "ambiguous", "observed": null, "source": "amplicon",
          "copies": [
            {
              "region": "1", "start": 13134, "end": 13299, "strand": -1, "variant_position": 13196, "identity": 99.56,
              "gap_compressed_identity": 99.56, "aligned_length": 229, "observed": "TCAAAG", "flank_edits": 1,
              "call": "ref", "anchors": 2, "ortholog": null, "source": "amplicon"
            },
            "…"
          ],
          "…": "…"
        },
        "…"
      ],
      "…": "…"
    },
    "…": "…"
  },
  "…": "…"
}
```

A running job carries `partial: true` and the genotyping block it has filled so far. Immediately after the reference
stage there is one `genomes` entry, the reference, and every summary count is still zero — so a UI must not read an
empty summary as "nothing found".

<!-- example: response GET /primers/check/{job_id} 200 capture=capture-check-genotyping-running.json#/response -->
```json
{
  "job_id": "afa2448f22719708b2d5dda511b3d08b", "status": "running", "partial": true,
  "progress": { "done": 1, "total": 120, "stage": "reference", "running": [] },
  "results": {
    "genotyping": {
      "algorithm_version": "g1",
      "summary": {
        "genomes_total": 0, "ref": 0, "alt": 0, "other": 0, "ambiguous": 0, "missing": 0, "unavailable": 0
      },
      "…": "…"
    },
    "…": "…"
  },
  "…": "…"
}
```

**`other` in practice.** No genome was called `other` in either full panel run so far — not on this deletion, and not on
the five-unit TGG repeat `1:13735:TTGG:T` (`rs5413863413`), which was chosen precisely because assemblies there carry
differing repeat counts. That variant came closest: `sorghum_tx430nano` has a copy at `Scaffold_20:32561494-32561655`
reading `ATTGGTGGTGGTGGTGGT`, one TGG short of either haplotype, so that **copy** is `other`; its second copy reads `ref`,
which makes the genome `ambiguous`. So `other` is reachable per copy, and a genome-level `other` needs every orthologous
copy to miss both haplotypes.

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
- Mongo **connection** failures during a check (no genes collection, including a 15 s timeout getting it, or a genes
  query that fails with an error, in any stage) also end in exit 75: the job fails internally with `MONGO_UNAVAILABLE`,
  is requeued at the front of its queue, and the restarted worker runs it again, so no check finishes with
  `ANNOTATION_UNAVAILABLE` because mongo refused or dropped its connections (a mongod that hangs is different: see the
  next bullet). A job that is still hit on its last attempt (`check.max_attempts`, 2) ends with error
  `MONGO_UNAVAILABLE` (1 h error TTL; a re-POST runs it again) instead of being requeued forever. The driver holds
  queries while it reconnects (30 attempts, 1 s apart), so after mongod stops or restarts, the job's log lines (prefixed
  `primers worker [<site_key>] job <id>`) can first show the retry line of the next bullet, then, about 30 s after the
  loss, `primers check: annotation unavailable: failed to reconnect after 30 attempts with interval 1000 ms` (or
  another driver error, such as `Topology was destroyed`), then
  `primers worker: fatal MONGO_UNAVAILABLE; requeueing running jobs and exiting`. A mongod that is back within those
  30 s answers the held queries and the job carries on. A retry line alone therefore does not prove that a query was
  slow.
- A genes query that **times out** does not exit the worker by itself (the self-heal restart below is the only
  exception). The worker sends at most 5 genes queries to mongo at a time
  (the driver's 5 connections), and the others wait in the worker, in order. Each query has a 15 s budget counted from
  dispatch, so that wait counts. A query that runs out of time is retried once, with another 30 s: the first attempt
  keeps running and, if it had been sent, the query is sent once more; the first to settle decides (an error from
  either is a connection failure). The job log shows
  `primers check: mongo genes query timed out after 15000 ms; retrying it once with a 30000 ms budget`, or, for a query
  that had not been sent yet,
  `primers check: mongo genes query timed out after 15000 ms waiting for one of the 5 genes query slots; waiting once more with a 30000 ms budget`.
  If the retry gets an answer, nothing else changes. If it also times out, the log shows
  `primers check: annotation unavailable: mongo genes query timed out after 30000 ms`, and the job stops querying
  mongo and completes with `ANNOTATION_UNAVAILABLE`, `details: {cause: "timeout"}`.
  - **What changes in such a result.** Genome-target specificity verdicts and pan-genome statuses do not use genes and
    are unchanged. From the failure on, `genes` and `ortholog` are null, so in gene mode the pan-genome `primary`
    product and the genotyping allele calls, which prefer annotated orthologs, may differ. In transcript mode, cDNA
    products whose transcripts were not mapped yet are grouped by transcript id without the isoform suffix, and the
    cDNA-target verdicts and statuses count those groups. The reference or genomes whose annotation finished before the
    failure keep their genes.
  - **Slow or hung.** A timed-out query is not cancelled on the server. It keeps its slot until mongo answers it, so
    the queries of later jobs wait in the worker instead of piling up behind it on mongo's connections. A mongod that
    stops answering but keeps its connections open (a stalled disk, `SIGSTOP`, memory pressure) looks exactly like a
    slow one to a single query. The driver fails the queries stuck on such connections only after its 6 min socket
    timeout (errors that no check is waiting for any more) and reconnects, but the new connections hang the same way:
    the 5 slots stay held or fill up again, and every later check would wait about 45 s for mongo and complete with
    `{cause: "timeout"}`. The self-heal restart replaces the worker's connections instead of waiting on them; it
    cannot make a hung mongod answer (see the end of the next bullet).
  - **Self-heal restart.** Each time a job gives up on a genes query that timed out, the worker also looks at all of
    its genes queries (all jobs of the process) and treats mongo as lost, exactly like a connection failure, when
    either holds:
    - **hung:** all 5 query slots are held and no genes query has been answered (or has failed) for 5 min, counted
      from the later of the last answer and the moment a query was sent while no other query was outstanding, so
      time without checks does not count. Full slots alone never qualify: a large healthy pan-genome check keeps all
      5 busy for minutes, but its queries are answered every few milliseconds.
    - **repeated_timeouts:** this is the third job in a row whose genes queries timed out, with no genes query
      answered successfully in between (a late answer to an abandoned query also resets the count). The driver's
      6 min socket-timeout errors are not answers, so this still fires when those errors have restarted the 5 min
      count of the hung condition.

    With a hung mongod the worker therefore restarts at the third check in a row that times out, or earlier at the
    first check that times out while all 5 query slots are held and no genes query has been answered (or has failed)
    for 5 min; the checks that timed out before it completed with `{cause: "timeout"}`. With little traffic (small
    checks that never hold all 5 slots) the restart comes at the third such check. The worker's log then shows, after
    that job's `annotation unavailable: mongo genes query timed out` warning, these lines (the first and third at
    error level):

    ```
    primers worker [<site_key>] job <id> primers check: mongo genes queries appear hung (no answer for 312 s, 5 of 5 query slots held, 2 consecutive timed-out jobs); treating mongo as unavailable so that the check worker restarts with fresh connections
    primers worker [<site_key>] job <id> requeued (MONGO_UNAVAILABLE)
    primers worker [<site_key>] job <id>: fatal error MONGO_UNAVAILABLE (requeued): gene annotation (mongo) is unavailable (cause: hung; no genes query answered for 312 s with 5 of 5 query slots held)
    primers worker: fatal MONGO_UNAVAILABLE; requeueing running jobs and exiting
    ```

    or `mongo genes queries repeatedly timed out (...)` and
    `(cause: repeated_timeouts; genes queries timed out in 3 consecutive jobs)`, followed by `job <id> requeued (shutdown)`
    for each other running job. `no answer for` and `answered for` count the seconds since a genes query of the worker
    last succeeded (or since its first genes query); the driver's socket-timeout errors are not answers.

    The worker exits 75 and pm2 restarts it with fresh connections. Its running jobs go back
    to the front of their queues and run again after the restart, but the attempt each was on counts toward
    `check.max_attempts` (2): if that was the escalating job's second attempt, it ends with error `MONGO_UNAVAILABLE`
    (`(error)` instead of `(requeued)` in the fatal line, and the cause in the job's error message) instead of being
    requeued, and a requeued job that is hit by a fatal mongo error again ends with that error too.

    What happens after the restart depends on how mongod hangs. The restarted worker's startup probe only connects
    and gets the genes collection; it sends no genes query. If mongod no longer answers connections, the probe times
    out after 15 s (`primers worker: mongo probe failed: mongo probe timed out`, then
    `mongo is unavailable; exiting with code 75`) and pm2 keeps restarting the worker with its backoff, so no check
    runs and queued checks wait until mongod answers. If mongod still answers connections but not queries (a stalled
    disk can do that), the probe passes and the restarted worker runs checks again: they complete with
    `{cause: "timeout"}`, and at the latest every third one in a row restarts the worker again, so the self-heal lines
    recur. The two `mongosh` commands below tell this case apart (the ping answers, the genes query does not).
    Checks that completed with `{cause: "timeout"}`, before or between restarts, stay cached (below): drop them once
    mongo is fixed.
  - **Recurring timeouts.** When `{cause: "timeout"}` warnings or self-heal restarts recur, first check that mongod
    answers:

    ```bash
    timeout 10 mongosh --quiet --eval 'db.adminCommand({ping: 1})'
    timeout 20 mongosh --quiet sorghum11 --eval 'db.genes.find({_id: "SORBI_3005G072200"}, {_id: 1}).maxTimeMS(5000).toArray()'
    ```

    If it does not answer, restart mongod. If it answers, check the genes indexes (overlap queries use
    `{location.map, location.region, location.start, location.end}`) and the mongo host's disk. If all of that is
    fine and checks still time out, `pm2 restart sorghum_primers11` drops the worker's connections.
  - **Cached for 24 h.** A check that completed with `{cause: "timeout"}` is a normal `done` job: identical requests
    share it and get the degraded result for 24 h. Once mongo is healthy, drop each affected job (its id is in the log
    prefix) so that the next POST of that request runs it again; until then `GET /primers/check/<id>` answers
    `404 UNKNOWN_JOB`:

    ```bash
    K=primers:sorghum_v11:sorghum11
    redis-cli -p 6380 -n 1 DEL $K:job:<id> $K:result:<id> $K:partial:<id>
    redis-cli -p 6380 -n 1 ZREM $K:finished <id>
    ```

  - A query of the same job that fails with an error, not a timeout, is still fatal (exit 75), also after a timeout.
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
| `NTTHAL` | `ntthal` | `ntthal` binary for genotyping designs (default `/home/olson/primer3-2.6.1/bin/ntthal`); not checked at startup, a missing binary is `503 THERMO_UNAVAILABLE` |
| `PRIMERS_VARIATION_URL` | `variation.base_url` | Ensembl REST base (default `https://data.gramene.org/pansite-ensembl-115`); `https://` only, or `http://127.0.0.1[:port]` for a loopback test server; anything else is ignored with a warning |
| `PRIMERS_VARIATION_ENABLED` | `variation.enabled` | `1`/`true`/`yes`; anything else switches Ensembl off (see [Feature flags](#feature-flags)) |
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

**Genotyping keys** (defaults in `config.js` and `config/default.yaml`; see [Genotyping primers](#genotyping-primers-kasp--allele-specific-pcr)):

| Key | Default | Meaning |
| --- | --- | --- |
| `ntthal` | `/home/olson/primer3-2.6.1/bin/ntthal` | `ntthal` binary (env `NTTHAL`) |
| `variation.enabled` | `true` | Ensembl variants on or off (env `PRIMERS_VARIATION_ENABLED`) |
| `variation.base_url` | `https://data.gramene.org/pansite-ensembl-115` | Ensembl REST base (env `PRIMERS_VARIATION_URL`) |
| `variation.release` | `'115'` | reported as `source.release` and in `GET /primers/genomes` |
| `variation.species` | `{sorghum_bicolor: sorghum_bicolor}` | `system_name` → Ensembl species; no other genome is ever sent to Ensembl |
| `variation.timeout_ms` | 8000 | per outbound call |
| `variation.chunk_bp` | 10000 | overlap chunk size |
| `variation.max_chunk_bytes`, `variation.max_lookup_bytes` | 1000000, 262144 | response body caps |
| `variation.max_window` | 50000 | largest `GET /primers/variants` window |
| `variation.list_limit_default`, `variation.list_limit_max` | 2000, 5000 | `limit` default and maximum |
| `variation.max_concurrent`, `variation.queue_wait_ms` | 4, 5000 | outbound limiter per API process |
| `variation.cache_entries` | 500 | LRU size of the result cache and of the failure cache |
| `variation.cache_ttl_ms`, `variation.unknown_cache_ttl_ms`, `variation.unavailable_cache_ttl_ms` | 3600000, 300000, 30000 | how long results, unknown ids and failures are cached |
| `variation.breaker_failures`, `variation.breaker_window_ms`, `variation.breaker_open_ms` | 3, 60000, 60000 | circuit breaker |
| `variation.retry_after_s`, `variation.queue_retry_after_s` | 30, 5 | `retry_after_s` of upstream failures and of `queue_full` |
| `variation.max_allele_length` | 50 | longest allele |
| `variation.max_shift` | 1000 | farthest an indel may slide (`VARIANT_TOO_REPETITIVE`) |
| `variation.ems_source_pattern` | `^EMS_` | Ensembl sources counted as EMS |
| `genotyping.template_flank` | 400 | smallest template flank |
| `genotyping.num_return_per_run` | 20 | pairs per Primer3 design run |
| `genotyping.num_sets_default`, `genotyping.max_sets` | 6, 10 | `assay.num_sets` default and maximum |
| `genotyping.max_scored_per_orientation` | 8 | candidates scored per orientation |
| `genotyping.max_primer3_runs`, `genotyping.max_thermo_calls` | 54, 272 | process budget of one design |
| `genotyping.thermo_concurrency`, `genotyping.thermo_timeout_ms` | 4, 5000 | `ntthal` calls at a time per design, and each call's timeout |
| `genotyping.guard_gap` | 10 | `SEQUENCE_TARGET` guard length beside the zone |
| `genotyping.mask_exempt_pad` | 36 | allele-specific window kept out of the repeat mask |
| `genotyping.as_min_tm`, `genotyping.as_min_gc` | 52, 15 | hard floors |
| `genotyping.structure_warn_th`, `genotyping.structure_high_th` | 47, 55 | structure Tm thresholds (warn, high) |
| `genotyping.as_tm_diff_warn` | 1.0 | `AS_TM_IMBALANCE` threshold (high above a fixed 2.0) |
| `genotyping.common_tm_low`, `genotyping.common_tm_high` | -1.0, 3.0 | `COMMON_TM_OUT_OF_RANGE` window |
| `genotyping.neighbour_3p_window` | 5 | 3′ window for known neighbours |
| `genotyping.dense_window`, `genotyping.dense_count` | 30, 2 | `DENSE_NEIGHBOURS` |
| `genotyping.check_max_sets`, `genotyping.check_max_unique_primers` | 5, 13 | caps of the proposed check request |
| `check.genotype_cpu_s_per_genome` | 0.2 | cost of the allele caller per genome task |
| `check.genotype_flank_min`, `check.genotype_amplicon_pad` | 15, 50 | allele-caller windows, also prepared at submit |
| `check.genotype_ortholog_min_identity`, `check.genotype_ortholog_size_tolerance`, `check.genotype_max_copies`, `check.genotype_max_anchors`, `check.genotype_max_megablast`, `check.genotype_megablast_timeout_ms`, `check.genotype_megablast_min_identity`, `check.genotype_megablast_min_query_cover`, `check.genotype_megablast_min_bitscore_frac` | 95, 0.2, 10, 50, 30, 20000, 95, 0.8, 0.9 | the worker's allele caller (spec §5.6), described with `results.genotyping` |
| `check.genotype_offlocus_max_mismatches` | 2 | an off-locus product changes an allele prediction (to `both`) only when each primer has at most this many mismatches; weaker products are reported in the `WEAK_OFF_TARGETS` warning and still appear in specificity and pan-genome results |

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

The thermodynamic tables are compiled in; do not create `primer3_config/`. Genotyping designs also run `ntthal` (config
`ntthal`) for tailed structures and deliberate-mismatch duplex Tm, with the design's salts. BLAST+ 2.13.0 lives in `/home/olson/bin`.
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
- `check_genotyping.test.js` also needs `PRIMERS_IT_WORKER=1`. It submits one genotyping check job to a dev API that has
  its own worker: the `check.request` of a live rs871475760 KASP design, over the 11 research assemblies (real BLAST,
  about 40 s). It then compares `results.genotyping` with the real-data expectations.
  - The job must be new, so delete the site's Redis keys before a re-run.
  - With `PRIMERS_IT_CAPTURE_DIR` set, it also writes the two captures that the
    [results example](#example-rs871475760-sets-s1-and-s2-over-11-assemblies) cites.
- `contract.test.js` also runs every fixture through the handlers' own request rules (`check/normalize.js` up to the
  catalog lookup, `design.normalize`), and requires the check fixtures together to send every check param.
- `docs_examples.test.js` checks every example of [Genotyping primers](#genotyping-primers-kasp--allele-specific-pcr): requests validate with sway and the
  handlers' request rules, responses against their swagger definitions (allowing `null` where a definition is
  `x-nullable`), and each example against the recording it was taken from (`test/primers/fixtures/docs/`). It also
  requires the configuration section to list every environment override and genotyping config key.
- Manual checks (not automated): stop the worker with Ctrl-C during a pan-genome job and start it again — the job ends
  `done` with `attempts: 2`; start the API with `PRIMERS_REDIS_URL=redis://localhost:6399` — `POST /primers/check`
  answers `503 JOB_STORE_UNAVAILABLE` within 3 s while `/primers/design` works; start it with
  `PRIMERS_FASTA_ROOT=/nonexistent` — designs answer JSON errors (`422 NO_SEQUENCE`) and the process stays up.
