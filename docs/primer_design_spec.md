# Primer design: full v1 design spec

> Generated 2026-09-12 by a design workflow (3 designers, 3 adversarial critics, synthesis, completeness check).
> **The approved implementation plan takes precedence.** Where this spec disagrees, apply these overrides:

- Tests run with `node --test "test/primers/unit/**/*.test.js"` and `".../integration/**/*.test.js"`. Node 24 does not expand directory arguments.
- `local_max_jobs: 2`, so a specificity job isn't stuck behind a pan-genome job.
- Automatic `max_product_size` raise and `PRODUCT_TOO_LONG_TO_CHECK`.
- Rules for trusting `product_tm`, and the fixed benchmark gate.
- `Cache-Control` middleware and the bare-`/primers` 404 go **before** `register`.
- Explicit `js-yaml@3.15.0` and `sway@1.0.0` dependencies.
- `--no-track` worktree and `--ignore-scripts` npm installs.
- cDNA `-max_target_seqs` from the cDNA DB count.
- `amplifies` definition and summary invariant.
- Uppercase primer sequences in design responses.
- Tests for the dual-prefix resolver dirs.
- `soft_masked` detected by sampling for lowercase.
- `PANGENOME_TRANSCRIPT_MODELS_ONLY`.
- Restart-recovery and Redis-down become **manual** checks, not integration tests.
- The worker runs as a separate pm2 app (`worker_main.js`), not a supervisor fork. Drop `supervisor.js` and the app.js worker start; `site_key` is derived from the swagger.yaml `basePath` in both processes.
- gramene-search Parcel hack, `--ignore-scripts`, and no `^1.0.0` dependency before publish.
- Header check uses `curl -s -D - -o /dev/null`, not `-I`.
- Deploy and rollback steps as in M5.

Follow-up decisions:
- Repeat avoidance uses the real soft-mask where one exists (4 sorghum assemblies) and a BLAST copy-number mask everywhere else.
- The check worker runs as a separate pm2 app.
- Checks are anonymous, protected by host-wide caps.

## Changes after implementation (check algorithm version 2)

This spec records the v1 design, and the code has changed since. **`docs/primer_design_api.md` is authoritative for
current behaviour.** These are the deltas from the sections below:

- **Check stringency (user decision):**
  - New check param `max_amplifying_mismatches`: 0–5, default 3, and it must be below `ignore_mismatches`.
  - A product amplifies only if each primer has at most that many edits and neither primer is 3′-blocked.
  - Products with a 4- or 5-edit primer are `unlikely`.
  - `ALGORITHM_VERSION` is `'2'`.
- **Candidates:**
  - A product must be at least as long as its two primers together.
  - Products whose re-aligned primer footprints overlap are discarded.
  - When the left and right primers are identical, their products are labelled `LR`.
- **Sites:**
  - A site whose gap-aware lower bound exceeds the cap is not FASTA-re-aligned. It gets an approximate `mm` and `approx: true`.
  - The two BLAST halves of a 1-edit indel site, each with ≥ 11 aligned bases, are joined into one gapped site.
  - Re-aligned coordinates are clamped to the sequence.
- **Pan-genome:**
  - `summary.truncated` counts genomes whose search hit a cap; their statuses are lower bounds.
  - The transcript-to-gene fallback handles `_T<n>`, `.t<n>` and `.<n>.v<x>.<y>` ids.
- **Worker:**
  - `BLAST_TIMEOUT` is not retried.
  - If gene annotation fails at the connection level (no collection, or a query error) while `ctx.fatalOnMongoUnavailable` is set, the run throws `MONGO_UNAVAILABLE` (`fatal`). The worker requeues the job and exits 75.
  - A genes query that times out is retried once with twice the budget, and the first attempt stays in the race. A second timeout is not fatal: the job completes with `ANNOTATION_UNAVAILABLE` (`details.cause: "timeout"`).
  - At most 5 genes queries per mongo handle are outstanding on the driver. A timed-out query keeps its slot until mongo answers it, so abandoned queries cannot pile up on the connections across jobs.
- **Security:** BLAST+ and Primer3 run with `NCBI_DONT_USE_NCBIRC=1` in private 0700 working directories, never `/tmp`.
- **Design:**
  - `n_mask` designs recompute `product_tm` on the unmasked template (Primer3 `long_seq_tm`).
  - A multi-record FASTA paste adds the warning `MULTIPLE_RECORDS`.
  - `template_only` skips the product-range vs template-length check.
  - `ASSEMBLY_MISMATCH` ignores synthetic map bins such as UNANCHORED.
- **Cost and sensitivity:**
  - The pan-genome BLAST term is multiplied by `check.pangenome_cpu_factor`: 2.0, against 1.97 measured on the full sorghum panel.
  - `results.sensitivity_note` is request-specific: it states the actual cap, word sizes and primer lengths.
- **§10.4 under v2:**

  | Case | Result |
  |---|---|
  | P2 | Two off-targets |
  | P1 | Two off-targets, one of them `unlikely` (L4/R0) |
  | P3 | `specific` |
  | Pan-genome P3 over 353/grassl/leoti | `single_mismatch` ×3 |
  | SORBI_3001G046200 pair | `no_amplicon` in sorghum_353 and sorghum_is12661, with `nearest` |
  | Full P3 panel (119 genomes) | 66 `single_perfect`, 51 `single_mismatch`, 1 `multiple`, 1 `no_amplicon`; `truncated` 0 |

- **gramene-search host:**
  - No `persistSequence={false}`. A pasted sequence stays in the in-memory store and is stripped only from saved-view snapshots.
  - `onGeneClick` opens gene links in a new tab.
  - The demo enables primers only when `PRIMERS_API` is set.

---

# Primer design for SorghumBase v11: implementation plan

## 1. Context

**Why.** SorghumBase v11 serves 120 sorghum assemblies and 8 other genomes. Users who need PCR or qPCR primers currently leave the site for Primer3web or Primer-BLAST, and those tools do not know these assemblies. v1 adds primer design, genome specificity checks and pan-genome coverage checks to the API and to the gene pages.

**Outcome: all five v1 scopes**

| Scope | Delivered by |
|---|---|
| (a) PCR primers from a gene (with flanks), a region or a pasted sequence; target/included/excluded regions; Tm/size/GC/product-size params | `POST /primers/design`, modes `gene`, `region`, `sequence` |
| (b) qPCR on spliced cDNA with a primer spanning an exon–exon junction | mode `transcript` plus `SEQUENCE_OVERLAP_JUNCTION_LIST` |
| (c) Primer-BLAST-style genome specificity check | `POST /primers/check` → worker BLASTs the reference dna DB, plus the reference cdna DB for qPCR designs |
| (d) Keep primers out of repeats | Real soft-mask where one exists (4 sorghum assemblies); a BLAST self-depth mask everywhere else |
| (e) Pan-genome coverage | Same worker, run over every same-species assembly |

**Deliverables**
- Swagger branch `primer-design`, off `sorghum_v11`, to be merged into `release70` later.
- Primer3 2.6.1 installed at `/home/olson/bin/primer3_core`.
- New npm package `gramene-primers`.
- A Primers tab in gramene-search (GitHub master).
- `details.primers: true` in sorghum-webapp once the endpoints are live.

## 2. Final decisions

**Environment and safety**
- All swagger work happens in the worktree `/usr/local/gramene/subsites/sorghum/v11/gramene-swagger-primers` on branch `primer-design`.
- The dev server runs as `HOST=127.0.0.1 PORT=50111`.
- The live checkout, the pm2 app `sorghum_swagger11` and port `:50011` are never touched. `package-lock.json` is never committed.
- Module root is `api/helpers/primers/`, with the check algorithm under `api/helpers/primers/check/`.
- One `primers:` config block.
- Tests use `node:test` with `should`; mocha is not installed.

**Engine and coordinates**
- `primer3_core` 2.6.1 is spawned once per request: `-strict_tags`, one Boulder-IO record on stdin, no shell. It stays a separate process (GPL-2.0).
- Every template coordinate is 1-based inclusive (`PRIMER_FIRST_BASE_INDEX=1`). Intervals are `[start,length]`.
- `PRIMER_RIGHT_j=pos,len` names the rightmost template base, so the right footprint is `[pos-len+1,pos]` and `product_size = right.pos - left.pos + 1`. Verified against genome sequence.
- `PRIMER_LIBERAL_BASE=1`, and the server converts IUPAC ambiguity codes to N.
- `PRIMER_PRODUCT_MIN_TM=0` and `PRIMER_PRODUCT_MAX_TM=150` force Primer3 to print `PRODUCT_TM` (verified in `print_boulder.c:501`). `product_tm` is still nullable.
- Primer size bounds are 15–36. Primer3's `MAX_PRIMER_LENGTH` is 36, and the check needs at least 15 nt, so every designed pair is checkable unless the user allows Ns.

**Data access**
- Sequence is read in-process with `@gmod/indexedfasta@4.0.6`.
  - npm `overrides` pin `@gmod/bgzf-filehandle@4.0.0` and `generic-filehandle2@2.0.7`.
  - Reverse complement is IUPAC-aware and preserves case.
- fastaIdx is not used:
  - It picks the wrong assembly in multi-prefix dirs.
  - It cannot read dna_sm.
  - Its reverse complement only handles uppercase.
- **One** assembly resolver (`assemblies.js`) serves design, genomes and the worker. The worker never re-resolves; it reads `job.resolved`.
- Exon coordinates are gene-relative, 1-based, in transcription order.
  - Junctions are always derived from the exon segments.
  - `exon_junctions` is absent for non-coding transcripts, because `dump_genes.js` only builds it per translation.

**Repeat avoidance**
- Genome is `soft_masked` (dna_sm size differs from dna): take the mask from dna_sm.
- Otherwise use a megablast self-depth mask (depth ≥ 3, HSP ≥ 50 bp, ≥ 85% identity).
- Default mode `n_mask`: masked bases become N and `PRIMER_MAX_NS_ACCEPTED=0`.
- Mode `three_prime`: `PRIMER_LOWERCASE_MASKING=1`.
- Verified today: of 84 sorghum `dna_sm.toplevel.fa.gz` files, only 4 are really masked (rio, tx2783pac, tx430nano, tx436pac). The other 80 are same-size copies, and 37 of the 121 sorghum dirs have no dna_sm at all.

**Check contract**
- One frozen `PrimerCheckRequest` (§A.2.3):
  - Tuning goes in a `params{}` object.
  - `checks` ⊆ {specificity, pangenome}, default `[specificity]`; specificity always runs.
  - ≤ 10 pairs and ≤ 20 unique primers; primers match `^[ACGTacgt]{15,36}$`.
  - A genome from another species → 400.
  - No `max_mismatches` alias and no `transcriptome` value; the cDNA target follows from `mode: transcript`.
- `expected` is honoured only in `gene` and `region` modes. In transcript mode the genome target never returns `on_target_missing`.
- No DELETE and no `cancelled` status in v1, because job ids are shared across identical requests.

**Check algorithm**
- Primer-BLAST-style scoring: reward 1, penalty −1, ungapped, e-value 30000, fixed `-searchsp 1.5e10`.
- **Word size 5 for the reference** (genome and cDNA) and **6 for the pan-genome**. Measured on sorghum:
  - ws7 misses a 20-mer with 2 mismatches and a 22-mer with 3 (not found at word size 7).
  - ws6 finds the 20-mer; ws5 finds both.
  - Single-thread CPU per primer·Gb: ws5 5.2 s, ws6 2.2 s, ws7 1.2 s.
- A gap-aware parse filter runs first. Primer-BLAST amplification rules are applied only after DP re-alignment, where a gap counts as a mismatch.
- cDNA amplicons are grouped by gene (isoforms collapsed) before counting or classifying.
- Every result carries a computed per-primer sensitivity guarantee.

**Jobs**
- Job store is Redis 6380 db 1 through `ioredis@^5.10.1`. Job ids are deterministic, results are gzipped, every key has a TTL, and each site keeps at most 500 finished jobs.
- Each site has two queues, `spec` and `pan`.
- Host-wide caps:
  - ≤ 2 running jobs, of which ≤ 1 pan-genome job.
  - ≤ 4 BLAST processes/threads per specificity job and ≤ 8 per pan-genome job, so ≤ 12 BLAST CPUs host-wide.
  - BLAST runs at nice 10.
- If Redis is down, checks fail closed.
- app.js forks the worker under a supervisor. On restart, running jobs are requeued (max 2 attempts). Partial results are visible while a job runs.

**API hygiene**
- Error and warning codes are UPPER_SNAKE. Handler errors are `{message, code, details}`. Swagger validator errors keep their own shape, `{message:'Validation errors', errors:[{code}]}`; that is where `INVALID_CONTENT_TYPE` appears.
- Every `/primers` response sends `Cache-Control: no-store`.
- 503 bodies carry `details.retry_after_s`, and `Access-Control-Expose-Headers: Retry-After` is set.
- Unknown paths and wrong methods under `/primers` return JSON 404/405.
- Every async entry point is promise-wrapped, because Node 24 throws on unhandled rejections.

**Front end**
- `gramene-primers`:
  - TypeScript, built as ESM + CJS + d.ts, with React 18 as a peer and zero runtime dependencies.
  - `gpr-`-scoped CSS injected at runtime, and also shipped as `dist/gramene-primers.css`.
  - A `mount()` helper; no IIFE bundle in v1.
- gramene-search work starts from GitHub master (v2.23.1, React 18), not the stale React-17 submodule. The tab offers all four modes. Dev linking uses a tarball.

## 3. Architecture overview

```
 www.sorghumbase.org (sorghum-webapp search_app)
   gramene-search GeneList ─ Primers tab ─ <PrimerDesigner> (npm gramene-primers)
        │ fetch(cache:'no-store'), CORS *
        ▼
 Apache gorgonzola (TLS, mod_cache) ──► squam:50011  pm2 sorghum_swagger11 (app.js, swagger-express-mw, /sorghum_v11)
   POST /primers/design ─► controllers/primers.js ─► helpers/primers/design.js
        semaphore(4) ─► template.js ─► sequence.js ─► @gmod/indexedfasta ─► /scratch/olson/fasta/<sys>/dna/<Prefix>.dna[_sm].toplevel.fa.gz
                     └► repeat_mask.js ─► blastn -task megablast (only if no real soft-mask)
        primer3.js ─► spawn /home/olson/bin/primer3_core  (Boulder-IO stdin → stdout)
   GET  /primers/genomes ─► genomes.js (mongo maps+taxonomy) + assemblies.js (single resolver)
   POST /primers/check ─► jobs/index.js ─► check/index.js normalize() ─► Redis 6380 db1 (queues, job/result keys)
   GET  /primers/check/{id} ┘                                                ▲  host-wide slot ZSETs
        │ child_process.fork (jobs/supervisor.js)                             │
        ▼                                                                     │
   jobs/worker_main.js ── claim (Lua) ───────────────────────────────────────┘
        check/index.js run(request, ctx)
          blast.js (blastn -num_threads 1 × ≤8)  ─► <Prefix>.dna.toplevel / <Prefix>.cdna.all BLAST DBs
          sites.js ─► amplicons.js ─► realign.js (indexed FASTA) ─► classify.js
          specificity.js / pangenome.js ─► annotate.js (mongo genes: overlap, transcript→gene, orthologs)
```

## 4. Step 0: Primer3 2.6.1 and dev environment

### 0.1 Build and install Primer3 (as olson; same layout as BLAST)
```bash
cd /home/olson/src
curl -fL -o primer3-2.6.1.tar.gz https://github.com/primer3-org/primer3/archive/refs/tags/v2.6.1.tar.gz
tar xzf primer3-2.6.1.tar.gz && cd primer3-2.6.1/src
make -j4                                          # gcc/g++ 8.5, links -lm only
make install PREFIX=/home/olson/primer3-2.6.1
ln -s ../primer3-2.6.1/bin/primer3_core /home/olson/bin/primer3_core
/home/olson/bin/primer3_core -about               # must print 2.6.1
```
- Do not create `primer3_config/` and do not set `PRIMER_THERMODYNAMIC_PARAMETERS_PATH`; the thermodynamic tables are compiled in.
- **Runtime gate** before fixing `max_template_length`: a 50 kb template with 20 pairs must return 20 `PRODUCT_TM` lines in under 10 s. Otherwise set `design.max_template_length: 20000`.
```bash
SEQ=$(samtools faidx /scratch/olson/fasta/sorghum_bicolor/dna/Sorghum_bicolor.Sorghum_bicolor_NCBIv3.dna.toplevel.fa.gz 4:7400001-7450000 | grep -v '>' | tr -d '\n')
time printf 'SEQUENCE_ID=bench\nSEQUENCE_TEMPLATE=%s\nPRIMER_TASK=generic\nPRIMER_FIRST_BASE_INDEX=1\nPRIMER_EXPLAIN_FLAG=1\nPRIMER_LIBERAL_BASE=1\nPRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT=1\nPRIMER_NUM_RETURN=20\nPRIMER_PRODUCT_SIZE_RANGE=100-1000\nPRIMER_PRODUCT_MIN_TM=0\nPRIMER_PRODUCT_MAX_TM=150\n=\n' "$SEQ" \
  | /home/olson/bin/primer3_core -strict_tags | grep -c '^PRIMER_PAIR_[0-9]*_PRODUCT_TM='
```

### 0.2 Worktree for branch `primer-design`
```bash
cd /usr/local/gramene/subsites/sorghum/v11/gramene-swagger
git status --porcelain                            # must print only: ?? package-lock.json
git fetch origin
git worktree add -b primer-design ../gramene-swagger-primers origin/sorghum_v11
cd ../gramene-swagger-primers
cp ../gramene-swagger/package-lock.json .         # never git add it; add "package-lock.json" to .gitignore on the branch
cp -a ../gramene-swagger/node_modules .           # node_modules/gramene-mongodb-config -> ../../gramene-mongodb-config resolves identically
source ~/.nvm/nvm.sh && nvm use 24.14.1
npm install --no-audit --no-fund @gmod/indexedfasta@4.0.6 ioredis@^5.10.1 config@1.31.0
# then add to package.json: "overrides": {"@gmod/bgzf-filehandle": "4.0.0", "generic-filehandle2": "2.0.7"} and run npm install again
```
Never run npm in `../gramene-swagger` or `../gramene-mongodb-config`; the live instance uses both directories.

### 0.3 Dev server (after the app.js change in §A.9)
```bash
cd /usr/local/gramene/subsites/sorghum/v11/gramene-swagger-primers
HOST=127.0.0.1 PORT=50111 SWAGGER_HOST=localhost:50111 SWAGGER_SCHEMES=http \
PRIMERS_SITE_KEY=sorghum_v11_dev PRIMERS_GLOBAL_MAX_JOBS=1 PRIMERS_MAX_QUEUED=5 \
node app.js                                        # plain node, not pm2
# from the workstation: ssh -L 50111:localhost:50111 -L 5174:localhost:5174 -L 1234:localhost:1234 squam.cshl.edu
# docs: http://localhost:50111/sorghum_v11/docs/?url=/sorghum_v11/swagger
```
- Dev jobs use the same Redis db 1 and the same host-wide slot keys as production; the site key `sorghum_v11_dev` separates the namespace. Dev jobs therefore count against the host cap.
- **Guard after every session:**
  - `git -C ../gramene-swagger status --porcelain` still prints only `?? package-lock.json`.
  - `pm2 describe sorghum_swagger11` shows the same pid and uptime.
  - `ss -ltnp | grep 50011` shows the same pid.
- Tests never `require('app.js')`, because it listens on `PORT||50011`. Integration tests call the dev URL.

---

## 5. Part A: gramene-swagger (worktree `gramene-swagger-primers`)

### A.1 Files

**Modified**

| Path | Change |
|---|---|
| `api/swagger/swagger.yaml` | Add 4 paths before `/{collection}` (tag `Primer design`, controller `primers`) and the definitions from §A.2 |
| `app.js` | Optional SWAGGER_HOST override, `HOST` bind, primers CORS expose, supervisor start, JSON 404/405 for `/primers` (§A.9) |
| `config/default.yaml` | `primers:` block (§A.7) |
| `package.json` | Dependencies `@gmod/indexedfasta 4.0.6`, `ioredis ^5.10.1`, `config 1.31.0`; `overrides`; scripts `test:primers` (`node --test test/primers/unit/`), `test:primers:it` (`node --test test/primers/integration/`) |
| `.gitignore` | Add `package-lock.json` |

**Created: API side (`api/controllers`, `api/helpers/primers`)**

| Path | Purpose |
|---|---|
| `api/controllers/primers.js` | `designPrimers`, `primerGenomes`, `submitPrimerCheck`, `getPrimerCheck`. Each is wrapped as `(req,res)=>Promise.resolve().then(()=>h(req,res)).catch(e=>sendError(res,e))` and sends `Cache-Control: no-store` |
| `api/helpers/primers/node_compat.js` | Idempotent `util.is*` polyfill, copied from app.js lines 3–25 |
| `api/helpers/primers/config.js` | Requires node_compat. Sets `process.env.NODE_CONFIG_DIR ||= path.resolve(__dirname,'../../../config')`. Merges DEFAULTS with `config.has('primers') ? config.util.toObject(config.get('primers')) : {}` and env vars (§A.7). Exposes `setBasePath()`, `siteKey()`, `_setForTests()` |
| `api/helpers/primers/errors.js` | `PrimerHttpError(status, code, message, details)`, `sendError(res, err)`; 503 errors put `retry_after_s` in `details` and in the `Retry-After` header |
| `api/helpers/primers/semaphore.js` | Counting semaphore with a bounded waiter queue and a wait timeout |
| `api/helpers/primers/genomes.js` | Catalog: maps (type genome) + taxonomy, species = first `rank:'species'` ancestor; `sameSpecies(system_name)` |
| `api/helpers/primers/assemblies.js` | The single resolver (§A.4.1) |
| `api/helpers/primers/sequence.js` | BgzipIndexedFasta LRU (32 handles), `regionLength`, `fetch`, IUPAC `revcomp` |
| `api/helpers/primers/coords.js` | Footprints, genomic and spliced mappers, junction overlap, mask run merging |
| `api/helpers/primers/template.js` | `buildTemplate(req, deps)` for the four modes |
| `api/helpers/primers/repeat_mask.js` | Soft-mask application and megablast depth mask with an LRU cache (§B.14) |
| `api/helpers/primers/boulder.js` | `serialize` (injection-guarded), `parse`, `parseExplain`, `extractPairs` |
| `api/helpers/primers/primer3.js` | Spawn with timeout and output caps; `version()` |
| `api/helpers/primers/design.js` | Normalize and validate, presets, build tags, orchestrate under semaphore and deadline, build the response |

**Created: job infrastructure (`api/helpers/primers/jobs`)**

| Path | Purpose |
|---|---|
| `jobs/index.js` | `submit(body)` (async normalize, resolve, id, cost guard, store.submit) and `status(id)` |
| `jobs/redis_store.js` | ioredis store |
| `jobs/memory_store.js` | Same interface in memory, for tests |
| `jobs/lua.js` | Submit, claim and lock-refresh scripts |
| `jobs/worker.js` | Lock, requeue, claim loop, `runJob`, `ctx` |
| `jobs/worker_main.js` | Forked child entry; signal handling |
| `jobs/supervisor.js` | Fork, restart backoff, in-process mode |

**Created: check algorithm (`api/helpers/primers/check`)**

| Path | Purpose |
|---|---|
| `check/index.js` | Plug-in interface (§A.8.5) |
| `check/blast.js` | BLAST argument arrays and streaming line parser |
| `check/sites.js` | Hit → site conversion and filters |
| `check/amplicons.js` | Pair forward/reverse-facing sites into candidates |
| `check/realign.js` | FASTA fetch and DP re-alignment |
| `check/classify.js` | Primer-BLAST rules |
| `check/annotate.js` | Mongo annotation: gene overlap, transcript→gene, orthologs |
| `check/specificity.js` | Reference verdicts |
| `check/pangenome.js` | Per-genome status |
| `check/cost.js` | Cost estimate |

**Created: tests and docs**

| Path | Purpose |
|---|---|
| `test/primers/unit/*.test.js`, `test/primers/integration/*.test.js`, `test/primers/fixtures/**` | See §A.12 |
| `docs/primer_design_api.md` | See §A.11 |

### A.2 Endpoints (all under basePath `/sorghum_v11`)

#### A.2.1 `POST /primers/design` (synchronous)

**Request fields.** Swagger `PrimerDesignRequest` has `additionalProperties:false` and `consumes: application/json`.

| Field | Type / limits | Notes |
|---|---|---|
| `mode` | enum gene, transcript, region, sequence (required) | |
| `gene_id`, `transcript_id` | string ≤ 255 | Equality lookups only |
| `system_name` | `^[a-z0-9_]+$` ≤ 128 | |
| `region` | `{region ≤255, start ≥1, end ≥1, strand 1\|-1}` | |
| `sequence` | string ≤ 60000 | FASTA headers allowed |
| `flank_up`, `flank_down` | int 0–10000 | gene mode |
| `target`, `included` | `[start,length]` | 1-based template coordinates |
| `excluded` | ≤ 50 × `[start,length]` | |
| `avoid_repeats` | bool | |
| `repeat_mask_mode` | enum `n_mask` (default), `three_prime` | |
| `junction_spanning` | bool | transcript mode, default true |
| `template_only` | bool | Return the template (with mask) and no pairs; skips Primer3 |
| `params` | PrimerDesignParams | Closed set; see the next table |

**PrimerDesignParams**

| Field | Limits |
|---|---|
| `opt_size`, `min_size`, `max_size` | int **15–36** |
| `opt_tm`, `min_tm`, `max_tm` | 30–90 |
| `opt_gc`, `min_gc`, `max_gc` | 0–100 |
| `max_tm_diff` | 0–30 |
| `max_poly_x` | 0–10 |
| `gc_clamp` | 0–5 |
| `max_end_stability` | 0–100 |
| `max_ns` | 0–5 |
| `salt_monovalent` | 0–1000 |
| `salt_divalent` | 0–100 |
| `dntp_conc` | 0–100 |
| `dna_conc` | 0–10000 |
| `num_return` | 1–20 |
| `min_3_prime_overlap_of_junction`, `min_5_prime_overlap_of_junction` | 1–20 |
| `product_size_ranges` | 1–10 × `[a,b]`, each 20–50000 |

**Mode requirements**

| mode | required | optional | template |
|---|---|---|---|
| gene | `gene_id` | `system_name` (must match the gene), `transcript_id` (feature overlay), flanks | Genomic gene span ± flanks, in transcription orientation |
| transcript | `gene_id` | `transcript_id` (default canonical), `junction_spanning` | Spliced cDNA of the transcript |
| region | `system_name`, `region` | | Genomic region on the given strand |
| sequence | `sequence` | `system_name` (validated; used later for checks) | Cleaned pasted sequence |

**Examples**
```jsonc
{"mode":"gene","gene_id":"SORBI_3001G000200","flank_up":200,"flank_down":100,"params":{"product_size_ranges":[[300,800]]}}
{"mode":"transcript","gene_id":"SORBI_3001G000200","junction_spanning":true}
{"mode":"region","system_name":"sorghum_bicolor","region":{"region":"1","start":11080,"end":15099,"strand":-1},"target":[499,50],"excluded":[[1000,40]]}
{"mode":"sequence","sequence":">amp\nATGGCCRYT...","system_name":"sorghum_bicolor","avoid_repeats":true}
```

**Response (200).** The coordinates and sequences in this example were verified against the genome; Tm, penalty and similar numbers are illustrative.
```jsonc
{
  "template": {"mode":"transcript","system_name":"sorghum_bicolor","gene_id":"SORBI_3001G000200","transcript_id":"SORBI_3001G000200.1",
    "region":"1","start":11180,"end":14899,"strand":-1,"length":1982,"seq":"…",
    "masked":false,"mask_source":null,"mask":[],"masked_fraction":0,
    "features":{"gene":null,
      "exons":[{"id":"EER93047-1","start":1,"end":397,"genomic":{"start":14503,"end":14899}},"…"],
      "cds":{"start":299,"end":1651},"junctions":[397,493,597,869,992,1081,1187,1285,1546,1630]}},
  "pairs":[{"rank":0,"penalty":0.41,"product_size":124,"product_tm":81.2,"compl_any_th":0.0,"compl_end_th":0.0,
    "left":{"seq":"ATTACATCAAATAGGCCTTG","start":856,"end":875,"len":20,"tm":59.1,"gc":35.0,"self_any_th":0,"self_end_th":0,"hairpin_th":0,"end_stability":3.1,"penalty":0.2,
            "junction":{"position":869,"overlap_5p":14,"overlap_3p":6},
            "genomic":{"region":"1","start":13395,"end":13650,"strand":-1,"blocks":[{"start":13395,"end":13400},{"start":13637,"end":13650}]}},
    "right":{"seq":"AACTTCTTTGTCGATCCATG","start":960,"end":979,"len":20,"tm":59.8,"gc":40.0,"…":"…","junction":null,
            "genomic":{"region":"1","start":13291,"end":13310,"strand":1,"blocks":[{"start":13291,"end":13310}]}},
    "product":{"start":856,"end":979,"genomic":{"region":"1","start":13291,"end":13650,"strand":-1},"genomic_size":360}}],
  "explain":{"left":{"raw":"considered 9120, low tm 3280, ok 1022","considered":9120,"low tm":3280,"ok":1022},"right":{"raw":"…"},"pair":{"raw":"considered 91, ok 5"}},
  "settings":{"preset":"qpcr","junction_spanning":true,"avoid_repeats":false,"repeat_mask_mode":null,"params":{"opt_size":20,"…":"…"}},
  "engine":{"primer3":"2.6.1"},
  "warnings":[]
}
```
- In sequence mode every `genomic` field is null.
- `genomic_size` is present only for transcript templates.
- `product_tm` may be null.

**Errors**

| Status | code | When |
|---|---|---|
| 400 | `{message:'Validation errors', errors:[{code:'INVALID_CONTENT_TYPE'\|'OBJECT_ADDITIONAL_PROPERTIES'\|'PATTERN'\|…}]}` | Swagger layer |
| 400 | `INVALID_REQUEST`, `INVALID_PARAMS`, `INVALID_SEQUENCE`, `TEMPLATE_TOO_LONG`, `REGION_OUT_OF_BOUNDS`, `INTERVAL_OUT_OF_BOUNDS`, `SYSTEM_NAME_MISMATCH`, `PRIMER3_INPUT_ERROR` (`details.primer3_error`) | Handler or Primer3 |
| 404 | `UNKNOWN_GENE`, `UNKNOWN_TRANSCRIPT`, `UNKNOWN_GENOME`, `UNKNOWN_REGION` | |
| 413 | body-parser | Body > 100 kb |
| 422 | `NO_SEQUENCE`, `AMBIGUOUS_ASSEMBLY` | |
| 500 | `PRIMER3_FAILED`, `GENE_STRUCTURE_MISMATCH`, `INTERNAL` | |
| 503 | `BUSY` (`retry_after_s: 5`), `PRIMER3_UNAVAILABLE`, `FEATURE_DISABLED`, `MONGO_UNAVAILABLE` | |
| 504 | `DEADLINE_EXCEEDED` | Total 45 s or Primer3 30 s |

#### A.2.2 `GET /primers/genomes?system_name=X`

`system_name` has the same pattern as above. Unknown → 404 `UNKNOWN_GENOME`.
```jsonc
{"system_name":"sorghum_bicolor","species":{"taxon_id":4558,"name":"Sorghum bicolor"},
 "counts":{"total":120,"with_blastdb":120,"with_cdna_blastdb":120},
 "genomes":[
  {"system_name":"sorghum_bicolor","display_name":"Sb bicolor BTx623 v3","taxon_id":4558006,"map_id":"GCA_000003195.3","is_query":true,
   "has_sequence":true,"has_blastdb":true,"has_cdna_blastdb":true,"repeat_masking":"unmasked_copy","total_bases":708735318,"warnings":[]},
  {"system_name":"sorghum_rio","…":"…","repeat_masking":"soft_masked"}]}
```
- Species comes from the first ancestor with rank `species`. Never compute it as `taxon_id/1000`.
- The query genome is listed first, then the rest by `display_name`.
- Filesystem paths are never returned.
- Assemblies are resolved with concurrency 8 and cached for 10 minutes.

#### A.2.3 `POST /primers/check` (frozen contract)
```yaml
PrimerCheckRequest:
  type: object
  required: [system_name, pairs]
  additionalProperties: false
  properties:
    system_name: { type: string, pattern: "^[a-z0-9_]+$", maxLength: 128 }
    mode: { type: string, enum: [gene, transcript, region, sequence] }        # default region
    gene_id: { type: string, maxLength: 255 }
    transcript_id: { type: string, maxLength: 255 }
    checks: { type: array, minItems: 1, maxItems: 2, items: { type: string, enum: [specificity, pangenome] } }   # default [specificity]
    genomes: { type: array, maxItems: 150, items: { type: string, pattern: "^[a-z0-9_]+$", maxLength: 128 } }
    params:
      type: object
      additionalProperties: false
      properties:
        max_product_size: { type: integer, minimum: 50, maximum: 10000 }     # 4000
        ignore_mismatches: { type: integer, minimum: 3, maximum: 6 }         # 6
        min_total_mismatches: { type: integer, minimum: 0, maximum: 6 }      # 2
        min_3p_mismatches: { type: integer, minimum: 1, maximum: 5 }         # 2
        three_prime_window: { type: integer, minimum: 3, maximum: 10 }       # 5
        include_unlikely: { type: boolean }                                  # false
        repeat_site_threshold: { type: integer, minimum: 1, maximum: 100 }   # 5
    pairs:
      type: array
      minItems: 1
      maxItems: 10
      items:
        type: object
        required: [id, left, right]
        additionalProperties: false
        properties:
          id: { type: string, pattern: "^[A-Za-z0-9_.:-]+$", maxLength: 64 }
          left: { type: string, pattern: "^[ACGTacgt]{15,36}$" }
          right: { type: string, pattern: "^[ACGTacgt]{15,36}$" }
          expected:
            type: object
            required: [region, start, end]
            additionalProperties: false
            properties: { region: { type: string, maxLength: 255 }, start: { type: integer, minimum: 1 }, end: { type: integer, minimum: 1 } }
```

**Handler rules** (async `check.normalize`)
- Pair ids must be unique, else `DUPLICATE_PAIR_ID`. At most 20 unique primers, else `TOO_MANY_PRIMERS`.
- `mode: transcript` requires `gene_id`. `expected` is stripped, with warning `EXPECTED_IGNORED`.
- Gene and transcript modes look up the gene doc:
  - `gene.system_name` must equal `system_name`, else `SYSTEM_NAME_MISMATCH`.
  - `transcript_id` defaults to the canonical transcript.
  - `resolved.gene` stores the location, transcript ids and ortholog ids grouped by system_name (one `genes.find({_id:{$in:orthologIds}},{fields:{system_name:1}})`).
- The reference genome needs a dna BLAST DB (else 422 `NO_BLASTDB`), plus a cdna DB in transcript mode.
- `genomes` applies only when checks include pangenome.
  - Omitted → every same-species genome with `has_blastdb` (or `has_cdna_blastdb` in transcript mode), minus `system_name`.
  - Not in the catalog → 404 `UNKNOWN_GENOME`.
  - Other species, or missing the required DB → 400 `GENOME_NOT_CHECKABLE` with `details.genomes`.
- Cost guard (§B.13): over the limit → 422 `JOB_TOO_LARGE {estimate_cpu_s, limit}`.
- `kind` is `pangenome` if checks include pangenome, else `specificity`.

**Deterministic job id**
```js
job_id = sha256(canonicalJSON({ v: 1, algo: check.ALGORITHM_VERSION, request: normalized,
                                dbs: { [sys]: assemblies[sys].fingerprint } })).hex.slice(0, 32)
```
`normalized` is built as follows:
- Primers uppercased; pair order kept.
- `checks` sorted and unique, with `specificity` added.
- `genomes` expanded and sorted.
- All params defaults filled.
- `mode` defaulted; transcript defaulted.

**Responses**
```jsonc
// 202 newly queued, 200 already existed (re-POST of an errored job re-queues it)
{"job_id":"9f2c0a4be1d34c7a8e5f00112233aabb","status":"queued","kind":"pangenome","queue_position":0,
 "progress":{"done":0,"total":4,"stage":"queued","running":[]},"estimate":{"cpu_s":16},"created_at":"2026-09-12T20:01:02.000Z","warnings":[]}
```
Errors:
- 400: validation, `INVALID_REQUEST`, `DUPLICATE_PAIR_ID`, `TOO_MANY_PRIMERS`, `GENOME_NOT_CHECKABLE`, `SYSTEM_NAME_MISMATCH`.
- 404: `UNKNOWN_GENOME`, `UNKNOWN_GENE`.
- 422: `JOB_TOO_LARGE`, `NO_BLASTDB`, `AMBIGUOUS_ASSEMBLY`.
- 503: `QUEUE_FULL` (`retry_after_s: 60`), `JOB_STORE_UNAVAILABLE`, `FEATURE_DISABLED`, `MONGO_UNAVAILABLE`.

#### A.2.4 `GET /primers/check/{job_id}`

`job_id` matches `^[0-9a-f]{32}$`.
```jsonc
{"job_id":"…","status":"running","kind":"pangenome","partial":true,"progress":{"done":37,"total":121,"stage":"pangenome","running":["sorghum_is929","sorghum_rio"]},
 "created_at":"…","started_at":"…","request":{"…normalized…":""},"warnings":[],"results":{"…partial, same schema as §B.12…":""}}
{"job_id":"…","status":"done","partial":false,"progress":{"done":121,"total":121,"stage":"done"},"finished_at":"…","results":{"…§B.12…":""}}
{"job_id":"…","status":"error","error":{"code":"JOB_TIMEOUT","message":"check exceeded 30 min"}}
```
- `status` is one of `queued`, `running`, `done`, `error`.
- Unknown or expired job → 404 `UNKNOWN_JOB`.

### A.3 Coordinate conventions (API contract)
- **Template:** 1-based inclusive, in template orientation. Intervals are `[start,length]`.
- **Footprints:** left `[pos, pos+len-1]`, right `[pos-len+1, pos]`.
- **Junction value `p`:** the boundary between 1-based cDNA bases `p` and `p+1`. It is passed verbatim with `PRIMER_FIRST_BASE_INDEX=1`.
  - A left primer `[a,b]` spans `j` when `j-a+1 ≥ min5` and `b-j ≥ min3`.
  - A right primer spans `j` when `j-a+1 ≥ min3` and `b-j ≥ min5`.
- **Genomic strand of a primer:** the genomic strand whose 5′→3′ sequence equals the primer. Left = `template.strand`; right = `-template.strand`.
- **`genomic.blocks`:** ascending `{start,end}` pieces, more than one when a cDNA primer spans a junction. `product.genomic` is the envelope of both primers' blocks.
- **Gene-relative position `p` to genomic:** `strand==1 ? loc.start+p-1 : loc.end-p+1`.

### A.4 Template construction (`template.js`)

#### A.4.1 Assembly resolution: `assemblies.resolve(system_name)`, cached 10 min
1. `system_name` must match the regex **and** be in the maps catalog, else 404.
   - `dir = path.join(fasta_root, system_name)`; assert `path.dirname(dir) === fasta_root`.
   - File names come only from `readdir`.
2. **Candidates** are the union of:
   - prefixes from `dir/dna/*.dna.toplevel.fa.gz` that also have `.fai` and `.gzi`;
   - prefixes from `dir/*.dna.toplevel.nal|.nin`, excluding volume files (`/\.\d\d\.nin$/`).
3. **Pick**, in order:
   1. `cfg.assembly_overrides[system_name]`.
   2. The unique candidate ending in `'.' + map._id`.
   3. The unique best score, where score = number of `maps.regions` entries (`{names[], lengths[]}`) found in the candidate's `.fai` with the same length.
   4. The only candidate.
   5. A tie where the tied candidates have identical `.fai` (name,length) sets and identical `.nsq` sizes: newest BLAST DB mtime, with warning `AMBIGUOUS_ASSEMBLY`. This covers vitis_vinifera 12X vs IGGP_12x: both .nsq 121,706,180 bytes, same regions.
   6. Otherwise → 422 `AMBIGUOUS_ASSEMBLY` until an override is configured.
4. **Always validate:** fewer than 80% of `maps.regions` matched by name and length → warning `ASSEMBLY_MISMATCH`. For reference, sorghum_bicolor matches 10 of 11; its UNPLACED bin does not match.
5. **Result object:**
   - `{system_name, taxon_id, display_name, map_id, prefix, dir}`
   - `fasta:{dna, dna_sm}`, each requiring `.fai` and `.gzi`
   - `blastdb:{dna: <dir>/<prefix>.dna.toplevel (.nal preferred, else .nin), cdna: <dir>/<prefix>.cdna.all (.nal|.nin) | null}`
   - `repeat_masking`: `soft_masked` if dna_sm exists and its size differs from dna; `unmasked_copy` if same size; else `absent`. `cfg.repeat_masking_overrides` wins.
   - `total_bases` = sum of `.fai` lengths; `num_sequences` = number of `.fai` lines (the API process never spawns blastdbcmd).
   - `fingerprint` = sha1 of the prefix plus `size:mtimeMs` of the dna `.fa.gz`, `.fai`, the dna `.nal` or `.nin`+`.nsq`, and the cdna `.nal` or `.nin`+`.nsq`.
   - `warnings`

#### A.4.2 Sequence access
- `sequence.fetch(fastaPath, region, start1, end1, strand)` calls `getSequence(region, start1-1, end1)`, preserves case, and reverse-complements (IUPAC-aware) on strand −1.
- `regionLength` undefined → 404 `UNKNOWN_REGION`. A short read → 400 `REGION_OUT_OF_BOUNDS`. Internal fetch cap is 2 Mb.

#### A.4.3 Mode `gene`
1. `genes.findOne({_id: gene_id}, {fields:{_id,name,system_name,taxon_id,location,gene_structure,homology}})`.
2. Resolve the assembly; no `fasta.dna` → 422 `NO_SEQUENCE`.
3. Extent, with `up = flank_up||0`, `down = flank_down||0`:
   - strand +1: `gStart=max(1,start-up)`, `gEnd=min(regionLen,end+down)`, `effUp=start-gStart`.
   - strand −1: `gStart=max(1,start-down)`, `gEnd=min(regionLen,end+up)`, `effUp=gEnd-end`.
   - Length > `max_template_length` → 400 `TEMPLATE_TOO_LONG`, with a hint to use transcript or region mode.
4. Features, using `transcript_id` or else the canonical transcript:
   - `gene {effUp+1, effUp+geneLen}`.
   - `exons` in rank order, each `exon.start/end + effUp`.
   - `cds`: cDNA → gene-relative through the exon segments, then `+effUp`.
   - Mapper: template `t` → genomic `strand==1 ? gStart+t-1 : gEnd-t+1`.
   - Verified: SORBI_3001G000200 with flanks 200/100 gives `1:11080-15099(-1)`, length 4020, exon1 at 201–597, CDS start t=499 at genomic 14601.

#### A.4.4 Mode `transcript`
1. Fetch the gene span on its strand, then slice each exon in `transcript.exons` order.
   - Segment: `{id, t_start: cum+1, t_end: cum+n}`.
   - Genomic range: strand +1 `[loc.start+e.start-1, loc.start+e.end-1]`; strand −1 `[loc.end-e.end+1, loc.end-e.start+1]`.
2. Assert `cdna.length === transcript.length`, else 500 `GENE_STRUCTURE_MISMATCH` (logged).
3. Junctions = cumulative segment ends except the last. Compare with `exon_junctions` only when that field exists.
4. `transcript.exons.length === 1` with junction spanning requested → design without the constraint and warn `SINGLE_EXON_TRANSCRIPT`.
5. More than 200 junctions (`PR_MAX_INTERVAL_ARRAY`): keep the junctions inside `included`, else the first 200, with warning `JUNCTIONS_TRUNCATED`.
6. Features: `exons[{id,start,end,genomic}]`, `cds` (null for non-coding), `junctions`.
7. Mapper: clip the span to each segment and map each piece to genomic blocks.
   - strand +1: `[g_start+lo-t_start, g_start+hi-t_start]`.
   - strand −1: `[g_end-(hi-t_start), g_end-(lo-t_start)]`.
   - Sort blocks ascending.
   - Verified: cDNA 856–875 → blocks `[13395–13400, 13637–13650]` on −1; cDNA 960–979 → `[13291–13310]`, right primer on +1.

#### A.4.5 Mode `region`
- `1 ≤ start ≤ end ≤ regionLength`, and length ≤ max.
- Same mapper as gene mode with `effUp=0`; `features={}`.

#### A.4.6 Mode `sequence`
1. Drop lines starting with `>`, then remove whitespace and digits.
2. What remains must match `/^[ACGTRYKMSWBDHVN]+$/i`, else 400 `INVALID_SEQUENCE`. Length must be 20–50000.
3. Non-ACGTN characters become N, with warning `IUPAC_CONVERTED`.
4. Case is kept only when `avoid_repeats` is set (lowercase = the user's own mask); otherwise uppercase.

#### A.4.7 `avoid_repeats`
See §B.14. The template records `masked`, `mask_source` (`softmask`, `blast_depth`, `user_lowercase` or null), `mask` (merged `[start,len]` runs) and `masked_fraction`. `template.seq` is always returned uppercase.

### A.5 Primer3 invocation

**Always sent (server-owned):**
```
SEQUENCE_ID=<sanitised [^A-Za-z0-9_.:-]→_>   SEQUENCE_TEMPLATE=<one line>   PRIMER_TASK=generic
PRIMER_PICK_LEFT_PRIMER=1   PRIMER_PICK_RIGHT_PRIMER=1   PRIMER_PICK_INTERNAL_OLIGO=0
PRIMER_FIRST_BASE_INDEX=1   PRIMER_EXPLAIN_FLAG=1   PRIMER_LIBERAL_BASE=1
PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT=1   PRIMER_THERMODYNAMIC_TEMPLATE_ALIGNMENT=0
PRIMER_PRODUCT_MIN_TM=0   PRIMER_PRODUCT_MAX_TM=150      (forces PRIMER_PAIR_i_PRODUCT_TM; non-constraining)
P3_FILE_FLAG=0   PRIMER_NUM_RETURN=<n>   PRIMER_PRODUCT_SIZE_RANGE=<a-b c-d>
=
```

**Conditional tags**

| API input | Boulder tag | Format |
|---|---|---|
| `target` | `SEQUENCE_TARGET` | `s,l` |
| `included` | `SEQUENCE_INCLUDED_REGION` | `s,l` |
| `excluded` | `SEQUENCE_EXCLUDED_REGION` | `s,l s,l …` |
| transcript + `junction_spanning` + junctions | `SEQUENCE_OVERLAP_JUNCTION_LIST` | derived junctions, space-separated |
| same | `PRIMER_MIN_3_PRIME_OVERLAP_OF_JUNCTION` (default 4), `PRIMER_MIN_5_PRIME_OVERLAP_OF_JUNCTION` (default 7) | int |
| same | `PRIMER_INTERNAL_MIN_3_PRIME_OVERLAP_OF_JUNCTION`, `PRIMER_INTERNAL_MIN_5_PRIME_OVERLAP_OF_JUNCTION` | Set to the primer values, to pass Primer3's internal-oligo check |
| `avoid_repeats`, `n_mask` | `PRIMER_MAX_NS_ACCEPTED=0` (forced) | Masked bases replaced with N |
| `avoid_repeats`, `three_prime` | `PRIMER_LOWERCASE_MASKING=1` | Masked bases lowercase |
| size | `PRIMER_OPT_SIZE`, `PRIMER_MIN_SIZE`, `PRIMER_MAX_SIZE` | int |
| Tm | `PRIMER_OPT_TM`, `PRIMER_MIN_TM`, `PRIMER_MAX_TM` | float |
| GC | `PRIMER_OPT_GC_PERCENT`, `PRIMER_MIN_GC`, `PRIMER_MAX_GC` | float |
| `max_tm_diff` | `PRIMER_PAIR_MAX_DIFF_TM` | float |
| `max_poly_x` | `PRIMER_MAX_POLY_X` | int |
| `gc_clamp` | `PRIMER_GC_CLAMP` | int |
| `max_end_stability` | `PRIMER_MAX_END_STABILITY` | float |
| `max_ns` | `PRIMER_MAX_NS_ACCEPTED` | int |
| salt / dNTP / DNA | `PRIMER_SALT_MONOVALENT`, `PRIMER_SALT_DIVALENT`, `PRIMER_DNTP_CONC`, `PRIMER_DNA_CONC` | float |

**Presets** (filled before the client's `params` are merged; effective values echoed in `settings`)

| Preset | Modes | Size opt/min/max | Tm opt/min/max | GC | Other |
|---|---|---|---|---|---|
| `pcr` | gene, region, sequence | 20/18/25 | 60/57/63 | 30–70 | max ΔTm 3, poly-X 4, product `[[100,1000]]`, 5 pairs |
| `qpcr` | transcript | 20/18/24 | 60/58/62 | 35–65 | max ΔTm 2, poly-X 4, product `[[70,150]]`, 5 pairs, `junction_spanning` true |

**Cross-field checks** → 400 `INVALID_PARAMS`
- min ≤ opt ≤ max for size, Tm and GC.
- Every range has `a < b`, and at least one range has `a ≤ template length`.
- Target, included and excluded intervals lie inside the template.
- When a junction list is sent, `min_5` and `min_3` must each be ≤ `floor(max_size/2)`. Primer3 otherwise aborts with a global error (`libprimer3.cc:7494-7533`).

**Serializer**
- Keys must match `/^[A-Z0-9_]+$/`.
- Values are finite numbers or pre-cleaned strings. A value containing `\n`, `\r` or a leading `=` throws.
- Output is `k=v\n` lines followed by `=\n`.

**Parser**
- Read lines up to the first `=`, splitting each on its first `=`.
- `parseExplain` turns `considered 9120, low tm 3280, ok 1022` into `{raw, considered, 'low tm', ok}`.
- For `i < PRIMER_PAIR_NUM_RETURNED`:
  - `[lp,ll]=PRIMER_LEFT_i`, `[rp,rl]=PRIMER_RIGHT_i`.
  - Assert `PRIMER_PAIR_i_PRODUCT_SIZE == rp-lp+1`.
  - Read `_TM`, `_GC_PERCENT`, `_SELF_ANY_TH`, `_SELF_END_TH`, `_HAIRPIN_TH`, `_END_STABILITY` and `_PENALTY` for each primer, and `_PENALTY`, `_COMPL_ANY_TH`, `_COMPL_END_TH`, `_PRODUCT_SIZE` and `_PRODUCT_TM` (nullable) for the pair.

**Process (`primer3.js`)**
- `spawn(cfg.primer3_core, ['-strict_tags'], {cwd: os.tmpdir(), env:{PATH:'/usr/bin:/bin'}})`; `stdin.end(record)`; EPIPE on stdin is swallowed.
- **Order:** take the design semaphore (4 running, 16 waiters, 10 s wait, else 503 `BUSY`) **before** building the template. One 45 s deadline covers template, mask and Primer3. Primer3 gets `min(30 s, remaining)`, then SIGKILL → 504 `DEADLINE_EXCEEDED`.
- **Caps:** stdout 2 MB → SIGKILL + 500; stderr 64 KB, logged only.
- **Exit handling:**
  - `ENOENT`/`EACCES` → 503 `PRIMER3_UNAVAILABLE`.
  - Non-zero exit → 500 `PRIMER3_FAILED`; the log gets the exit code and 500 characters of stderr.
  - `PRIMER_ERROR` → 400 `PRIMER3_INPUT_ERROR`.
  - `PRIMER_WARNING` → warning `PRIMER3_WARNING`.
  - Zero pairs → 200 with `pairs:[]` and warning `NO_PAIRS`.
- `version()` runs `-about` lazily (2 s timeout) and caches the result.

### A.6 Validation, limits, security

| Limit | Value |
|---|---|
| Template length | ≤ 50,000 (drops to 20,000 if the Step 0 gate fails) |
| Flank | ≤ 10,000 |
| Internal genomic fetch | ≤ 2 Mb |
| num_return | ≤ 20 |
| Product size ranges | ≤ 10 |
| Excluded intervals | ≤ 50 |
| Body size | ≤ 100 kb |
| Design concurrency per instance | 4 running + 16 waiting |
| Design deadline | 45 s |
| Check: pairs / unique primers / genomes | 10 / 20 / 150 |
| Queued jobs per site (both queues) | ≤ 50 |
| Running jobs host-wide | 2, of which ≤ 1 pan-genome |
| BLAST per specificity job / pan-genome job | 4 / 8 (≤ 12 host-wide) |
| Estimated job CPU | ≤ 6,000 CPU-s |
| Job timeout | 30 min |
| Stored result size | ≤ 5 MB uncompressed |

- **Two validation layers.** Swagger enforces types, `additionalProperties:false`, patterns and maxima. The handler adds mode-conditional required fields, cross-field ranges, intervals-in-template, and catalog/filesystem membership.
- **Filesystem.** Regex check plus catalog membership before any `path.join`; the dirname is asserted; file names only from `readdir`.
- **Mongo.** `typeof === 'string'`, equality lookups only.
- **Processes.** Argument arrays, minimal env, no shell. The worker's `ctx.spawn` allow-list is `cfg.blastn`; the API side may spawn `cfg.primer3_core` and `cfg.blastn` (megablast mask).
- **Auth.** None (same as `/gene_lists/validate`). CORS stays `*`; `Retry-After` is exposed.

### A.7 Config

`config/default.yaml`, appended as a top-level block:
```yaml
primers:
  enabled: true
  site_key: null                    # default "<basePath w/o slash>:<getMongoConfig().db>" e.g. sorghum_v11:sorghum11
  primer3_core: /home/olson/bin/primer3_core
  blastn: /home/olson/bin/blastn
  fasta_root: /scratch/olson/fasta
  tmp_dir: /tmp
  assembly_overrides: {}            # { system_name: prefix }
  repeat_masking_overrides: {}      # { system_name: soft_masked|unmasked_copy|absent }
  design: { deadline_ms: 45000, primer3_timeout_ms: 30000, max_concurrent: 4, max_waiting: 16, wait_timeout_ms: 10000,
            max_template_length: 50000, max_flank: 10000, max_fetch_length: 2000000, max_stdout_bytes: 2000000 }
  repeat_mask: { enabled: true, min_hsp_len: 50, perc_identity: 85, evalue: 1.0e-10, min_depth: 3, threads: 2,
                 chunk: 20000, exon_pad: 100, timeout_ms: 20000, cache_entries: 200 }
  check:
    store: redis                    # redis | memory
    redis_url: redis://localhost:6380/1
    in_process: false
    global_max_jobs: 2
    pangenome_max_jobs: 1
    local_max_jobs: 1
    spec_job_procs: 4
    pangenome_job_procs: 8
    nice: 10
    max_queued: 50
    max_pairs: 10
    max_unique_primers: 20
    max_genomes: 150
    max_job_cpu_s: 6000
    cpu_s_per_primer_gb: { ws5: 5.2, ws6: 2.2, ws7: 1.2 }
    reference_word_size: 5
    pangenome_word_size: 6
    blast_timeout_ms: 600000
    job_timeout_ms: 1800000
    heartbeat_ms: 10000
    stale_ms: 60000
    queue_sweep_ms: 600000
    progress_min_interval_ms: 5000
    max_attempts: 2
    ttl_active_s: 21600
    ttl_done_s: 86400
    ttl_error_s: 3600
    max_finished_jobs: 500
    max_result_bytes: 5000000
    max_candidates_per_pair: 5000
    max_realign_sites_per_genome: 20000
    realign_concurrency: 32
    max_offtargets_listed: 100
    max_annotated_amplicons: 200
    defaults: { max_product_size: 4000, ignore_mismatches: 6, min_total_mismatches: 2, min_3p_mismatches: 2,
                three_prime_window: 5, include_unlikely: false, repeat_site_threshold: 5 }
```
- **Environment overrides**, read in `config.js` rather than by the config library:
  - `PRIMERS_ENABLED`, `PRIMERS_IN_PROCESS`: booleans, `/^(1|true|yes)$/i`.
  - `PRIMERS_SITE_KEY`, `PRIMER3_CORE`, `BLASTN`, `PRIMERS_FASTA_ROOT`, `PRIMERS_JOB_STORE`, `PRIMERS_REDIS_URL`.
  - `PRIMERS_GLOBAL_MAX_JOBS`, `PRIMERS_MAX_QUEUED`: integers.
- **No `custom-environment-variables.yaml`.** A malformed `__format: json` value would crash startup of the live runner.

### A.8 Job infrastructure

#### A.8.1 Redis layout

Prefix `P = primers:<site_key>:`.

| Key | Type | Content | TTL |
|---|---|---|---|
| `P job:<id>` | string JSON | `{job_id, kind, status, request, resolved, created_at, started_at, finished_at, progress, attempts, worker, warnings, error}`; `resolved` is never returned | queued/running 6 h, refreshed; done 24 h; error 1 h |
| `P partial:<id>` | gzip buffer | Partial results, ≤ 1 write per 5 s | 6 h, deleted on completion |
| `P result:<id>` | gzip buffer | Final results | 24 h |
| `P queue:spec`, `P queue:pan` | list | Job ids | none (ids only) |
| `P running` | zset | id → heartbeat ms | none |
| `P finished` | zset | id → finish ms; trimmed to 500 by deleting the oldest job/result keys | none |
| `P worker` | string | Worker lock, `SET NX PX 30000`, refreshed every 10 s | 30 s |
| `primers:global:slots` | zset | `<site_key>:<id>` → heartbeat | stale entries pruned at 60 s |
| `primers:global:pan_slots` | zset | Same, for pan-genome jobs | pruned at 60 s |

**Lua scripts** (`ioredis.defineCommand`; they derive job keys inside the script, which is fine on standalone Redis 5)
```lua
-- primersSubmit KEYS: job, targetQueue, otherQueue   ARGV: jobJson, id, ttlActiveS, maxQueued
local cur = redis.call('GET', KEYS[1])
if cur and cjson.decode(cur).status ~= 'error' then return {'EXISTS', cur} end
if redis.call('LLEN', KEYS[2]) + redis.call('LLEN', KEYS[3]) >= tonumber(ARGV[4]) then return {'FULL', ''} end
redis.call('SET', KEYS[1], ARGV[1], 'EX', tonumber(ARGV[3])); redis.call('RPUSH', KEYS[2], ARGV[2])
return {'QUEUED', ARGV[1]}

-- primersClaim KEYS: queueSpec, queuePan, running, globalSlots, panSlots
--              ARGV: nowMs, staleMs, globalMax, localMax, panMax, siteKey, jobKeyPrefix
local now = tonumber(ARGV[1])
redis.call('ZREMRANGEBYSCORE', KEYS[4], '-inf', now - tonumber(ARGV[2]))
redis.call('ZREMRANGEBYSCORE', KEYS[5], '-inf', now - tonumber(ARGV[2]))
if redis.call('ZCARD', KEYS[4]) >= tonumber(ARGV[3]) or redis.call('ZCARD', KEYS[3]) >= tonumber(ARGV[4]) then return false end
local function pop(q)
  while true do
    local id = redis.call('LPOP', q); if not id then return nil end
    local j = redis.call('GET', ARGV[7] .. id)
    if j and cjson.decode(j).status == 'queued' then return id end      -- drop expired/duplicate ids
  end
end
local id, pan = pop(KEYS[1]), false
if not id and redis.call('ZCARD', KEYS[5]) < tonumber(ARGV[5]) then id = pop(KEYS[2]); pan = id ~= nil end
if not id then return false end
redis.call('ZADD', KEYS[3], now, id); redis.call('ZADD', KEYS[4], now, ARGV[6] .. ':' .. id)
if pan then redis.call('ZADD', KEYS[5], now, ARGV[6] .. ':' .. id) end
return {id, pan and 'pan' or 'spec'}
```
- `complete(id)`: one MULTI doing SET result (gzip) EX ttl_done, SET job EX ttl_done, DEL partial, ZREM running/slots/pan_slots, ZADD finished; then trim `finished`.
- `fail(id)`: the same shape with ttl_error.
- `queuePosition(id)`: `LRANGE <queue> 0 49` and `indexOf` (Redis 5 has no LPOS).

#### A.8.2 Redis clients
- **API process:** `new Redis(url, {lazyConnect:true, enableOfflineQueue:false, maxRetriesPerRequest:1, connectTimeout:2000})`.
  - Add an `'error'` listener.
  - `await connect()` once; reset on failure.
  - Every store call has a 3 s deadline; failure → 503 `JOB_STORE_UNAVAILABLE`.
  - Design and genomes endpoints do not depend on Redis.
- **Worker:** default reconnect. Any loop error → sleep 5 s. BLAST never runs without a claimed slot, so the worker fails closed.

#### A.8.3 Supervisor and worker
- **app.js start:** `Promise.resolve().then(() => require('./api/helpers/primers/jobs/supervisor').start({basePath})).catch(e => console.error('primers worker not started', e))`.
- **Supervisor:**
  - `!enabled` → no-op.
  - `store: memory` or `in_process` → run the worker in-process.
  - Otherwise `fork('jobs/worker_main.js', [], {env:{...process.env, PRIMERS_SITE_KEY: siteKey}, execArgv:['--max-old-space-size=2048'], stdio:'inherit'})`, with `'error'` and `'exit'` listeners.
  - Restart with 5→60 s backoff; give up (logged) after 10 crashes in 10 min.
  - Exit code 75 (mongo unavailable) is restarted with backoff.
  - Parent `exit` → `child.kill('SIGTERM')`.
- **Worker loop:**
  1. `acquireWorkerLock`; if another live worker holds it, retry every 5 s.
  2. On first acquisition, `requeueRunning`: for each running id, if `attempts ≥ max_attempts` → `fail(WORKER_LOST)`; else status `queued`, LPUSH it to its queue, ZREM its slots.
  3. Every 10 min, EXPIRE every queued job key to `ttl_active`, so queued jobs don't expire while waiting.
  4. Every 1 s, `claim`, then `runJob`.
- **runJob:**
  1. Set `status: running`, `started_at`, `attempts+1`.
  2. Start an `AbortController` with the 30 min timeout and a 10 s heartbeat (ZADD running/slots/pan_slots, EXPIRE job).
  3. `results = await check.run(request, ctx)`.
  4. Serialized size > 5 MB → fail `RESULT_TOO_LARGE`; otherwise `complete`.
  5. On error → `fail({code: e.code || (aborted ? 'JOB_TIMEOUT' : 'CHECK_FAILED'), message})`.
  6. `finally`: clear timers, remove the tmpdir.
- **Shutdown** (SIGTERM, SIGINT, parent disconnect): abort the running job, set it back to `queued`, LPUSH it to the front of its queue, ZREM slots, exit within 1 s. After a hard kill, the next lock holder requeues within about 30 s. Clients see running → queued → running → done without re-POSTing.
- **Mongo failure:** if `mongoCollection()` resolves undefined in the worker, exit 75 so the supervisor restarts it. gramene-mongodb-config caches a failed connection forever.

#### A.8.4 site_key
- Default `"${basePath.slice(1)}:${require('gramene-mongodb-config').getMongoConfig().db}"`: `sorghum_v11:sorghum11` here, `sorghum_v11:search70` in release70.
- Logged at worker start. Dev uses `sorghum_v11_dev`.

#### A.8.5 Plug-in interface (`check/index.js`)
```js
module.exports = {
  ALGORITHM_VERSION: '1',                                               // bump => new job ids
  async normalize(body, { catalog, assemblies, mongo, cfg }) -> { request, resolved, kind, warnings, estimate: { cpu_s, total } },
  async run(request, ctx) -> results                                    // §B.12 schema
};
ctx = {
  jobId, siteKey, signal /* AbortSignal */, config /* frozen primers cfg */, resolved /* job.resolved */,
  procs /* 4 specificity job | 8 pangenome job */,
  progress({ done, total, stage, running }),     // ≤ 1 write / 5 s (+ final)
  partial(results),                              // ≤ 1 write / 5 s to P partial:<id>
  tmpdir(), log: { info, warn, error },
  spawnLines(cmd, args, { stdin, timeoutMs, onLine }) -> Promise<{ code, signal, stderrTail }>,  // allow-list, niced, SIGTERM→SIGKILL 3 s on abort
  mongo: require('gramene-mongodb-config')
};
```

### A.9 app.js changes (all guarded; unchanged behaviour when env is unset)
1. Before `SwaggerExpress.create`: if `process.env.SWAGGER_HOST` is set, load `api/swagger/swagger.yaml` with js-yaml `safeLoad`, set `host` and `schemes = (SWAGGER_SCHEMES||'http').split(',')`, and pass the object as `config.swagger`.
2. In the create callback, before `register`: `app.use(basePath + '/primers', cors({ exposedHeaders: ['Retry-After'] }))`.
3. After `register`:
   - `require('./api/helpers/primers/config').setBasePath(basePath)`, then the supervisor start from §A.8.3.
   - `app.use(basePath + '/primers', (req, res) => res.status(404).set('Cache-Control','no-store').json({message:'unknown primers endpoint', code:'NOT_FOUND'}))`.
   - `app.use(basePath + '/primers', (err, req, res, next) => res.headersSent ? next(err) : res.status(err.statusCode||err.status||500).set('Cache-Control','no-store').json({message: err.message, code: err.code || (err.statusCode===405 ? 'METHOD_NOT_ALLOWED' : 'INTERNAL'), errors: err.errors}))`.
4. `app.listen(port, process.env.HOST)`.

### A.10 Swagger additions
- **Paths:**
  - `/primers/design` (post `designPrimers`)
  - `/primers/genomes` (get `primerGenomes`)
  - `/primers/check` (post `submitPrimerCheck`; responses 200/202/400/404/422/503)
  - `/primers/check/{job_id}` (get `getPrimerCheck`)
  - All use `x-swagger-router-controller: primers` and `consumes: [application/json]`, and are inserted before `/{collection}`.
- **Request definitions** are strict: `PrimerDesignRequest`, `PrimerDesignParams`, `PrimerInterval`, `PrimerRegion`, `PrimerCheckRequest`.
- **Response definitions** are documentation only (`validateResponse:false`): `PrimerError`, `PrimerDesignResponse`, `PrimerGenomesResponse`, `PrimerCheckJob`, `PrimerCheckResults` (§B.12).

### A.11 `docs/primer_design_api.md` (style of `docs/gene_lists_api.md`)
1. Base path and engine; no auth; `Content-Type: application/json` required; 100 kb body limit.
2. Error shapes: the validator's `errors[].code` including `INVALID_CONTENT_TYPE`, and handler codes.
3. Coordinate conventions (§A.3).
4. Design: modes, fields, params↔tags table, presets, examples, response, errors, limits.
5. Repeat avoidance, with an honest note: `soft_masked` only for rio, tx2783pac, tx430nano and tx436pac; `blast_depth` elsewhere also masks multi-copy gene families.
6. Genomes endpoint.
7. Check endpoints:
   - lifecycle, deterministic ids, TTLs, restart behaviour, queue limits, cost guard;
   - verdict table (§B.8), results schema;
   - sensitivity guarantee: word-size runs; mid-primer indels not detected.
8. Operations:
   - config and env vars;
   - Redis keys and `redis-cli -p 6380 -n 1 --scan --pattern 'primers:<site_key>:*'` inspection;
   - Primer3 install, dev server, tests.

### A.12 Tests (`node:test` + `should`)

**Unit tests** (offline; `npm run test:primers`)

| File | Covers |
|---|---|
| `boulder.test.js` | Serializer rejects newlines, lowercase keys and leading `=`. Parser handles `=` inside values. `extractPairs` on the captured `primer3_output_transcript.txt` (right footprint, `product_size === rpos-lpos+1`, `product_tm` present). `parseExplain` |
| `coords.test.js` | SORBI_3001G000200 (−): cDNA 856–875 → blocks `[13395–13400, 13637–13650]`; 960–979 → `[13291–13310]` on +1; gene template t=499 → 14601. SORBI_3001G000700 (+) equivalents. Junction overlap formulas |
| `template.test.js` | Stub sequence source: cDNA length; derived junctions equal `exon_junctions` when present; non-coding multi-exon fixture (AT3G05735.1-like, 2 exons, no `exon_junctions`) still gets junctions; flank clamping; `n_mask`; `IUPAC_CONVERTED`; `TEMPLATE_TOO_LONG` |
| `design_validate.test.js` | Mode-required fields; cross-field params; junction overlap > `max_size/2` → `INVALID_PARAMS`; intervals outside the template; FASTA cleaning |
| `assemblies.test.js` | Fixture dirs in `os.tmpdir()`: map._id suffix; region scoring; vitis-style tie with identical .fai/.nsq → newest + warning; non-identical tie → `AMBIGUOUS_ASSEMBLY`; `.nal` preferred over `.nin`; `unmasked_copy` detection; `../x` rejected |
| `genomes_species.test.js` | 39946/1736656/1736659 → 4530; 112509001 → 4513; 4558006 → 4558 |
| `jobs_memory.test.js` | Idempotent submit; FULL; claim limits (global, local, pan); stale pruning; claim drops expired ids; `requeueRunning` + `max_attempts`; queue TTL sweep; TTL expiry → 404; progress and partial visible via `status()` |
| `job_id.test.js` | Key order and default filling do not change the id; a fingerprint change does |
| `sites.test.js` | Plus/minus conversion; tail bounds; `del18` (q1–17, t3=4) and `ins12` (q1–12, t3=11) survive the parse filter |
| `realign.test.js` | P3_L in 353 footprint → positions `[1]`; P2_L on chr5 `GGACAGCTCCACGACCTATCAG` → `[10,7]`; `del18` → mm 1 with one gap near the 3′ end |
| `amplicons_classify.test.js` | LR/RL/LL/RR pairing; size bounds; candidate cap; rule matrix including terminal mismatch → `likely_weak` |
| `verdicts.test.js` | Mode × target table (§B.8): transcript genome target never `on_target_missing`; cDNA grouping by gene; pan-genome status matrix; primary ordering; `mismatch_in_3p_window` |
| `sensitivity.test.js` | `guaranteed_max_mismatches(L=20,W=5)=3`, `(20,6)=2`, `(22,6)=2`, `(25,5)=4` |
| `controller_wrap.test.js` | A rejected handler promise becomes a JSON error, never an unhandled rejection |
| `contract.test.js` | Loads `swagger.yaml` with sway and validates every JSON in `test/primers/fixtures/contract/requests/*.json` (generated by gramene-primers `npm run fixtures`); fails on any Validation errors. Also asserts §A.2.3 example bodies validate, and legacy bodies (top-level `max_mismatches`, `checks:['transcriptome']`) are rejected |

**Integration tests** (`PRIMERS_IT_BASE=http://localhost:50111/sorghum_v11 npm run test:primers:it`)

| File | Covers |
|---|---|
| `design_endpoints.test.js` | All modes. For every returned primer, fetch its `genomic.blocks` from fastaIdx `:8888` and assert they equal the primer (reverse complement by strand). Transcript mode: ≥ 1 primer per pair satisfies the junction overlap rule. Error cases from §10 |
| `genomes_endpoint.test.js` | sorghum_bicolor: 120 genomes, query first, sorghum_rio `soft_masked` |
| `check_jobs.test.js` | Fixtures from §10.4: idempotent 202 then 200; results; restart recovery; Redis-down behaviour |

---

## 6. Part B: check worker (`api/helpers/primers/check/`)

### B.1 Inputs
`job.resolved` holds:
- `assemblies[sys]` from §A.4.1 for the reference and every pan-genome genome;
- `gene` (gene/transcript modes): location, transcript ids, and ortholog ids by system_name from `homology.homologous_genes.ortholog_one2one`, `ortholog_one2many`, `ortholog_many2many` (the real keys);
- species info.

### B.2 Tasks and scheduling
```
tasks = [ {kind:'reference', target:'genome'} ]
      + (mode=='transcript' ? [ {kind:'reference', target:'cdna'} ] : [])
      + (checks has 'pangenome' ? genomes.sortBy(total_bases, system_name).map(g => ({kind:'pangenome', genome:g, target: mode=='transcript' ? 'cdna' : 'genome'})) : [])
progress = {done, total: tasks.length, stage: 'reference'|'transcriptome'|'pangenome'|'annotate', running:[system_name]}
```
- Reference tasks run first, as one blastn with `-num_threads ctx.procs`. Results go to `ctx.partial` as soon as the specificity block exists.
- Pan-genome tasks then run as ≤ `ctx.procs` concurrent single-thread blastn processes, in a stable order so repeat jobs reuse the page cache.
- A BLAST failure is retried once after 2 s. A second failure marks that genome `error`, but a reference failure fails the job.
- `reference_size` for pan-genome = size of the on-target (or inferred) reference amplicon; null if none.

### B.3 BLAST
```js
const args = ['-task','blastn-short','-reward','1','-penalty','-1','-word_size', String(ws /* 5 reference, 6 pangenome */),'-ungapped',
  '-evalue','30000','-searchsp','15000000000','-dust','no','-soft_masking','false',
  '-max_target_seqs', String(Math.max(5000, asm.num_sequences)), '-max_hsps','100000',
  '-num_threads', String(threads), '-db', target === 'cdna' ? asm.blastdb.cdna : asm.blastdb.dna, '-query','-',
  '-outfmt', target === 'cdna' ? '6 qseqid sseqid qlen qstart qend sstart send sstrand mismatch qseq sseq'
                               : '6 qseqid sseqid qlen qstart qend sstart send sstrand mismatch'];
```
- stdin is `>q<i>\n<SEQ>\n`, one entry per unique uppercase primer.
- stdout is parsed line by line and never buffered; only the last 2 KB of stderr is kept. Timeout 600 s.
- Multi-volume DBs: pass the `.nal` prefix.
- `sseqid` equals the `.fai` name for genome DBs and the transcript id for cDNA DBs (no gene field in the title, so gene mapping goes through mongo; §B.9).

**Measured** (read-only, squam, BLAST+ 2.13.0; 15 primers including planted mutants vs sorghum_bicolor at 8 threads)

| word size | wall | hits | `m20_7_14` (2 mm, runs 6/6/6) | `m22_6_12_18` (3 mm) | `del18` / `ins12` |
|---|---|---|---|---|---|
| 7 | 5.9 s | 232,600 | missed | missed | q1–17 / q1–12 only |
| 6 | 8.5 s | 278,461 | found (mm 2) | missed | same |
| 5 | 15.3 s | 304,888 | found | found (mm 3) | same |

Single thread, 10 primers:
- ws5 on sorghum_bicolor: 37.2 CPU-s, i.e. **5.2 CPU-s per primer·Gb**.
- ws6 on sorghum_353: 14.6 CPU-s → **2.2**.
- ws7 on sorghum_353: 7.9 CPU-s → **1.2**.

**Sensitivity guarantee** (reported per primer)
- A substitution-only site is found if it contains an exact run of length ≥ W and its ungapped score is ≥ about 11 (score 12 was reported, 10 was not).
- `guaranteed_max_mismatches(L, W)` = the largest k with `ceil((L-k)/(k+1)) ≥ W` and `L - 2k ≥ 12`. Examples: 20-mer W5 → 3, W6 → 2; 22-mer W6 → 2; 25-mer W5 → 4.
- Indels are found when one side of the indel still scores ≥ 11. Mid-primer indels in short primers can be missed; the `sensitivity_note` says so.

### B.4 Hits to sites (`sites.js`)
```
t5 = qstart-1 ; t3 = qlen-qend
plus : face 'F'; p5 = sstart - t5; p3 = send + t3          // 3' end points to higher coordinates
minus: face 'R'; p5 = sstart + t5; p3 = send - t3
lbU  = mismatch + ceil(t5/2) + ceil(t3/2)                   // substitution-only lower bound (ungapped X-drop)
lbI  = mismatch + (t5>0) + (t3>0)                           // gap-aware bound: each unaligned tail costs ≥ 1 edit
parse filter (param-independent):  keep if lbU <= 5  ||  (lbI <= 5 && t5 + t3 <= 8)
job filter (before pairing):       lb = (t5+t3 <= 8) ? min(lbU, lbI) : lbU ;  drop if lb >= ignore_mismatches
                                   // the 'blocked' 3'-rule is NOT applied before re-alignment
near_perfect(primer) = count of hits with t5==0 && t3==0 && mismatch <= 1   (reference genome)
```
- Sites are stored as per-primer columnar `Int32Array`s. The parse filter kept about 62% of hits, i.e. about 13k sites per primer per sorghum genome.
- A primer over `max(50000, 60 × total_bases/1e6)` sites is marked `truncated: true, repetitive: true`.

### B.5 Amplicon calling (`amplicons.js`)
```
role(seq) = seq == pair.left ? 'L' : 'R'                    // identical L/R → 'L'
maxSize   = task.kind=='pangenome' && reference_size ? max(max_product_size, ceil(1.5*reference_size)+500) : max_product_size
for pair P:
  S = sites(P.left) ∪ sites(P.right), grouped by subject; F[] (face F), R[] (face R), both sorted by p5
  for subject, for f in F:
     for r in R starting at lowerBound(r.p5 >= f.p5 + max(len(f),len(r)) - 1) while r.p5 - f.p5 + 1 <= maxSize:
        cands.push({subject, fwd:f, rev:r, orientation: role(f)+role(r)})       // LR | RL | LL | RR
        if cands.length > max_candidates_per_pair (5000): P.truncated = true; break all
sitesToRealign = unique(c.fwd, c.rev)   (cap 20000/genome → truncated)
realign(sitesToRealign)                 // §B.6, realign_concurrency 32
for c: c.start = fwd.p5', c.end = rev.p5', c.size = end-start+1; c.class = classify(c); drop 'ignored'
target cdna → group by gene (§B.9)
```

### B.6 Re-alignment (`realign.js`)
```
genome target: window = FASTA(asm.fasta.dna)[subject][min(p5,p3) - 3 .. max(p5,p3) + 3]; face R → IUPAC revcomp (case-insensitive compare)
DP primer p[1..n] vs window w[1..m]: match 0, mismatch 1, gap 1; free leading/trailing window bases;
   primer fully consumed; last primer base aligned to a window base (no terminal gap); ties → end nearest m-3
traceback → edits with distance from primer 3' end (1 = terminal):
   mm = #edits, mm_3p = #edits with dist <= three_prime_window, terminal_mm = (dist 1 is a mismatch), gaps, p5' (genomic)
cdna target (cdna .fa.gz not indexed): exact mismatch positions over the aligned part from qseq/sseq,
   plus lbU for the tails; approx = (t5+t3 > 0)
```
The fetch uses the same `sequence.js` handle LRU; the samtools fallback is dropped.

### B.7 Classification (`classify.js`, Primer-BLAST rules)
```
blocked(s)  = s.mm >= min_total_mismatches && s.mm_3p >= min_3p_mismatches
ignored     : fwd.mm >= ignore_mismatches || rev.mm >= ignore_mismatches      (dropped)
unlikely    : blocked(fwd) || blocked(rev)                                    (counted; listed only if include_unlikely)
likely_weak : fwd.terminal_mm || rev.terminal_mm
likely      : otherwise
repetitive(primer) = near_perfect_sites > repeat_site_threshold
```
`likely` and `likely_weak` both count as amplicons.

### B.8 Verdicts by mode × target (`specificity.js`)

| mode | target | on-target | verdicts |
|---|---|---|---|
| gene, region (with `expected`) | genome | Candidate with `subject==E.region`, `\|start-E.start\|≤2`, `\|end-E.end\|≤2`, orientation LR/RL. `start`/`end` = 5′ of the forward-facing and reverse-facing primers, i.e. `product.genomic` | `specific`, `off_targets`, `on_target_missing`, `truncated`, `error` |
| gene/region without `expected`, sequence | genome | Exactly one candidate with both mm 0 and LR/RL → `on_target_inferred: true` | `specific`, `off_targets`, `unverified_target`, `truncated`, `error` |
| transcript | genome | None. Likely amplicons with both sites inside the gene span → `gdna_products[]` | `specific` (no likely amplicons outside the gene), `off_targets`, `truncated`, `error`; never `on_target_missing` |
| transcript | cdna (`results.transcriptome`) | Any likely amplicon on a transcript of `gene_id`, grouped | `specific`, `off_targets`, `on_target_missing`, `truncated`, `error` |

- `off_targets` also covers LL and RR single-primer products.
- `off_target_count` is always exact; the listed off-targets are capped at 100.

### B.9 cDNA grouping (`annotate.js`)
- Subject → gene mapping: `genes.find({taxon_id: asm.taxon_id, 'gene_structure.transcripts.id': {$in: subjects}}, {fields:{_id:1,'gene_structure.transcripts.id':1}})`.
  - Batched by 500; cached per genome for the job.
  - An unmapped subject becomes its own gene (id with `/\.\d+$/` stripped), with warning `TRANSCRIPT_GENE_UNMAPPED`.
- Group likely amplicons by `(gene_id, orientation)` into `{gene_id, isoforms:[{transcript_id,start,end,size}], size_min, size_max, likelihood (best), left_mm, right_mm, left_3p_mm, right_3p_mm, approx}`.
- Counts, verdicts and pan-genome statuses all use groups. Measured need: sorghum_bicolor has 6,852 multi-isoform genes; sorghum_353 and rio have none.

### B.10 Pan-genome (`pangenome.js`)

| status | Rule (genome target: distinct amplicons; cDNA target: distinct gene groups) |
|---|---|
| `db_unavailable` | No dna DB (cdna DB in transcript mode), or `AMBIGUOUS_ASSEMBLY` |
| `error` | BLAST failed twice or timed out |
| `no_amplicon` | 0 likely; `nearest` = best unlikely, if any |
| `single_perfect` | Exactly 1 likely; both mm == 0 |
| `single_mismatch` | Exactly 1 likely; any mm > 0; flags `mismatch_in_3p_window` (`left_3p_mm>0 \|\| right_3p_mm>0`) and `terminal_mismatch` |
| `multiple` | ≥ 2 likely; `primary` + `others[≤10]` + `other_amplicons` |

- **Primary order:** `ortholog` desc, then `left_mm+right_mm` asc, then `|size - reference_size|` asc. `size_delta = size - reference_size`.
- **`ortholog`:** true when an overlapping gene (or the cDNA group's gene) is in the design gene's ortholog set; null for region and sequence modes.
- **`ortholog_annotated`:** whether the design gene lists any ortholog in that genome. With `no_amplicon`, this separates "no ortholog annotated" from "primers fail".
- Coordinates are never compared across assemblies.
- **Per-pair summary:** `{genomes_total, single_perfect, single_mismatch, multiple, no_amplicon, db_unavailable, error, amplifies}`.

### B.11 Annotation
- Gene overlap for genome amplicons: `genes.find({'location.map': asm.map_id, 'location.region': subject, 'location.start': {$lte: end, $gte: start - 1000000}, 'location.end': {$gte: start}}, {fields:{_id,name,biotype,location}}).limit(20)`.
  - Uses the index `{location.map, location.region, location.start}`; measured 52 keys examined, 36 ms.
  - At most 200 amplicons per task.
- Mongo unavailable: `genes: null`, `ortholog: null`, warning `ANNOTATION_UNAVAILABLE`; classification still completes.

### B.12 Results schema (`PrimerCheckResults`)
```jsonc
{
  "engine": {"algorithm_version":"1","blast":"2.13.0","reference":"blastn-short r1 p-1 ws5 ungapped e30000 searchsp1.5e10","pangenome":"… ws6 …"},
  "params": {"max_product_size":4000,"ignore_mismatches":6,"min_total_mismatches":2,"min_3p_mismatches":2,"three_prime_window":5,"include_unlikely":false,"repeat_site_threshold":5},
  "reference": {"system_name":"sorghum_bicolor","map_id":"GCA_000003195.3","total_bases":708735318},
  "sensitivity_note": "Sites are detected only if they contain an exact match of at least the word size and an ungapped score of about 11 or more; mid-primer indels in short primers may be missed. See primers[].sensitivity.",
  "primers": {"GGACAGCTCCACAACATATCAG": {"len":22,"near_perfect_sites":2,"repetitive":false,"truncated":false,
               "sensitivity":{"reference":{"word_size":5,"guaranteed_max_mismatches":3},"pangenome":{"word_size":6,"guaranteed_max_mismatches":2}}}},
  "specificity": {"target":"genome","pairs":[{
     "id":"P2","verdict":"off_targets","on_target_inferred":false,"truncated":false,
     "on_target":{"region":"4","start":7423537,"end":7423746,"size":210,"orientation":"LR","likelihood":"likely","left_mm":0,"right_mm":0,"left_3p_mm":0,"right_3p_mm":0,"left_mm_pos":[],"right_mm_pos":[],"genes":[{"id":"SORBI_3004G087700","strand":1}]},
     "off_target_count":2,"unlikely_count":0,
     "off_targets":[
       {"region":"4","start":7437317,"end":7437526,"size":210,"orientation":"LR","likelihood":"likely","left_mm":0,"right_mm":0,"left_3p_mm":0,"right_3p_mm":0,"left_mm_pos":[],"right_mm_pos":[],"approx":false,"genes":[{"id":"SORBI_3004G087800","strand":1}]},
       {"region":"5","start":66890615,"end":66890824,"size":210,"orientation":"RL","likelihood":"likely","left_mm":2,"right_mm":2,"left_3p_mm":0,"right_3p_mm":0,"left_mm_pos":[10,7],"right_mm_pos":[15,10],"approx":false,"genes":[{"id":"SORBI_3005G183900","strand":-1}]}],
     "gdna_products":[]}]},
  "transcriptome": null,   // mode transcript: {"target":"cdna","pairs":[{"id","verdict","on_target":{"gene_id","isoforms":[…]},"off_target_count","off_targets":[{"gene_id","isoforms":[…],"size_min","size_max","likelihood","left_mm",…,"approx"}]}]}
  "pangenome": {"target":"genome","pairs":[{
     "id":"P3","reference_size":580,
     "summary":{"genomes_total":3,"single_perfect":0,"single_mismatch":0,"multiple":0,"no_amplicon":0,"db_unavailable":0,"error":0,"amplifies":3},
     "genomes":[{"system_name":"sorghum_353","display_name":"Sb verticilliflorum 353","status":"…","ortholog_annotated":true,
        "primary":{"region":"4","start":7499931,"end":7500510,"strand":1,"size":580,"size_delta":0,"orientation":"LR","likelihood":"likely_weak",
                   "left_mm":1,"right_mm":0,"left_3p_mm":1,"right_3p_mm":0,"left_mm_pos":[1],"right_mm_pos":[],"mismatch_in_3p_window":true,"terminal_mismatch":true,
                   "genes":[{"id":"353.004G093500"}],"ortholog":true},
        "other_amplicons":0,"others":[],"nearest":null}]}]},
  "warnings": [{"code":"…","message":"…"}],
  "timings_ms": {"reference": 11000, "sorghum_353": 1500}
}
```
- In the pan-genome example, the status and summary counts depend on the paralog site in 353 (illustrative). The primary coordinates and mismatches are verified.
- All warning codes are UPPER_SNAKE: `NO_FASTA_FOR_REALIGN`, `AMBIGUOUS_ASSEMBLY`, `ASSEMBLY_MISMATCH`, `TRANSCRIPT_GENE_UNMAPPED`, `ANNOTATION_UNAVAILABLE`, `EXPECTED_IGNORED`.

### B.13 Performance, concurrency, cost

**Concurrency**
- ≤ 2 running jobs host-wide across all sites and dev, of which ≤ 1 pan-genome job.
- A specificity job uses 4 BLAST threads/processes; a pan-genome job uses 8. Maximum 12 BLAST CPUs, at nice 10.
- Design megablast masking in the API processes adds at most 4 × 2 threads per instance.
- gramene-blast (`MAX_RUNNING_JOBS 3`) is not counted.

**Cost estimate** (`check/cost.js`; server is authoritative, and the client mirrors it for display)
```
est_cpu_s = uniq_primers × [ ref_Gb × c(ws_ref) + (transcript ? ref_cdna_Gb × c(ws_ref) : 0)
                             + Σ_pan (genome_Gb or cdna_Gb) × c(ws_pan) ]      c = {5: 5.2, 6: 2.2, 7: 1.2} CPU-s per primer·Gb
reject > 6000 → 422 JOB_TOO_LARGE
```
`genome_Gb` is `total_bases / 1e9`. The cDNA size comes from the cdna `.nin` header via a one-time `blastdbcmd -info` in the worker; use 0.15 Gb in the API estimate.

| Workload (sorghum) | Estimate |
|---|---|
| Specificity only, 5 pairs / 10 primers, ws5, 4 threads | ≈ 37 CPU-s → ≈ 15–20 s wall, plus ≈ 2 s re-align and annotate |
| + cDNA (transcript mode) | + ≈ 5 CPU-s |
| Pan-genome, all 119 other assemblies, 10 primers, ws6, 8 procs | ≈ 1,800 CPU-s → ≈ 4 min |
| Pan-genome, 20 primers | ≈ 3,600 CPU-s → ≈ 7.5 min (under the cap) |
| Pan-genome, transcript mode (cDNA ≈ 0.1 Gb each), 10 primers | ≈ 260 CPU-s → ≈ 35 s |
| Hits per primer per sorghum genome | 18–21k; parse filter keeps ≈ 62%; about 1–3k candidate sites to re-align per genome |

- **No site cache in v1.** Identical requests reuse results for 24 h through the deterministic job id.
- **Worker memory:** ≤ 8 genome tasks in flight × about 30 MB of columnar sites, well inside `--max-old-space-size=2048`.

### B.14 Repeat avoidance at design time (`repeat_mask.js`)
```
if !avoid_repeats → no mask
elif asm.repeat_masking == 'soft_masked':       fetch same coords from asm.fasta.dna_sm → lowercase runs → mask_source 'softmask'
                                                (transcript mode: splice the dna_sm gene sequence → mask carries through; dna/dna_sm .fai identical in 310/310 sets)
elif cfg.repeat_mask.enabled && asm.blastdb.dna:
   queries = genomic/region templates: 20 kb chunks of the template;   transcript mode: each exon ± 100 bp (genomic), projected onto cDNA
   blastn -task megablast -db <asm.blastdb.dna> -query - -dust no -soft_masking false -evalue 1e-10 -perc_identity 85
          -max_target_seqs 500 -max_hsps 200 -num_threads 2 -outfmt "6 qseqid qstart qend length"
   depth[i] = #HSPs with length >= 50 covering i (self-hit counts); masked where depth >= 3 → mask_source 'blast_depth'
   runs inside the design semaphore and deadline (timeout 20 s → warning REPEAT_MASK_FAILED, design proceeds unmasked)
   LRU cache key: fingerprint|region|start|end|strand|mode
else: warning NO_REPEAT_MASK
apply: n_mask (default) → masked bases → 'N' + PRIMER_MAX_NS_ACCEPTED=0 ;  three_prime → lowercase + PRIMER_LOWERCASE_MASKING=1
warnings: MOSTLY_REPEAT if masked_fraction > 0.8 ; BLAST_DEPTH_MASK ("multi-copy gene families are masked as repeats") when source blast_depth
```
- **Validation** (sorghum_tx436pac 4:7547610–7564601, 16,992 bp, against RepeatMasker): depth ≥ 3 gives precision 0.91, recall 0.60, 0.63 s at 2 threads.
- `SEQUENCE_EXCLUDED_REGION` is not used for masks (200-interval cap). `PRIMER_MASK_TEMPLATE` is not used (no plant k-mer lists).
- The check's `primers[].repetitive` flag is shown next to the pair in the UI.

### B.15 Error handling

| Situation | Behaviour |
|---|---|
| Reference DB missing at submit | 422 `NO_BLASTDB` |
| Pan-genome genome has no DB or an ambiguous assembly | `db_unavailable` |
| FASTA missing for re-alignment | Lower bounds with `approx: true`; warning `NO_FASTA_FOR_REALIGN` |
| blastn exit ≠ 0 | One retry, then genome `error` (reference: job error `REFERENCE_BLAST_FAILED`); stderr tail with paths reduced to basenames |
| Mongo unavailable in the worker | Exit 75 → supervisor restart → job requeued |
| Candidate or site caps hit | `truncated` |
| Abort (timeout or shutdown) | Children SIGTERM → SIGKILL after 3 s; job error `JOB_TIMEOUT`, or requeued on shutdown |

---

## 7. Part C: `gramene-primers` React package (`/home/olson/src/warelab/gramene-primers`, later `github.com/warelab/gramene-primers`, MIT)

### C.1 Packaging and build
- **Language and deps:** TypeScript strict (tbrowse tsconfig). Zero runtime dependencies. Peers `react`/`react-dom` `^18.2.0`. State with `useReducer` + context. Plain SVG and tables (no d3, ag-grid, react-bootstrap or zustand).
- **Build outputs** (`vite build`, library mode, `target: es2020`):
  - `dist/gramene-primers.js` (ESM), `dist/gramene-primers.cjs` (CJS)
  - `dist/index.d.ts`, plus `dist/index.d.cts` copied by `scripts/copy-dts.mjs`
  - `dist/gramene-primers.css`, emitted by a 10-line plugin from `src/styles/primers.css`, which is imported `?inline` for runtime injection
  - Externals: `react`, `react-dom`, `react-dom/client`, `react/jsx-runtime`
- **package.json:**
```json
{ "name":"gramene-primers","version":"0.1.0","type":"module",
  "main":"./dist/gramene-primers.cjs","module":"./dist/gramene-primers.js","types":"./dist/index.d.ts",
  "exports":{".":{"import":{"types":"./dist/index.d.ts","default":"./dist/gramene-primers.js"},
                  "require":{"types":"./dist/index.d.cts","default":"./dist/gramene-primers.cjs"}},
             "./style.css":"./dist/gramene-primers.css","./dist/*":"./dist/*","./package.json":"./package.json"},
  "sideEffects":["*.css"],"files":["dist","CHANGELOG.md"],
  "scripts":{"dev":"vite","build":"rm -rf dist && vite build && node scripts/copy-dts.mjs","typecheck":"tsc --noEmit",
             "test":"vitest run","test:it":"vitest run --config vitest.it.config.ts","fixtures":"vitest run test/fixtures.gen.test.ts",
             "lint:pkg":"publint && attw --pack .","pack:local":"npm run build && npm pack",
             "prepublishOnly":"npm run build && npm run typecheck && npm test && npm run lint:pkg"},
  "peerDependencies":{"react":"^18.2.0","react-dom":"^18.2.0"},"dependencies":{} }
```
  - No `source` field. Parcel follows `source` for symlinked packages.
- **devDependencies:** `vite ^7.3.6`, `@vitejs/plugin-react ^5.2.0`, `vite-plugin-dts ^4.5.4`, `typescript ~5.9.3`, `vitest ^4.1.11`, `jsdom ^27.4.0`, `@testing-library/react ^16.3.3`, `@testing-library/dom ^10.4.1`, `@testing-library/user-event ^14.6.1`, `@testing-library/jest-dom ^6.9.1`, `axe-core ^4.13.0`, `react`/`react-dom ^18.3.1`, `@types/react ^18.3`, `publint`, `@arethetypeswrong/cli`.
  - Do not take latest: npm on 2026-09-12 has vite 8.3, plugin-react 6.1 (vite 8 only), vitest 5, TypeScript 7.0.2.
  - jsdom 30 needs Node ≥ 24.15; the host runs 24.14.1.
- **Playground mode** (`vite`, `examples/playground`): port 5174; alias `gramene-primers → src/index.ts`; proxy `/sorghum_v11 → PRIMERS_PROXY_TARGET || http://localhost:50111`.
- **Releases and CI:**
  - `0.x` published with `--tag next`; `1.0.0` once `/primers` is live on data.sorghumbase.org.
  - `.github/workflows/ci.yml` on Node 20/24 runs typecheck, test, build and lint:pkg.

### C.2 Public API
```ts
// components
export { PrimerDesigner, PairsTable, TemplateMap, SpecificityResults, PangenomeMatrix, GenomePicker };
// headless
export { createPrimersClient, PrimersApiError, isAbortError };
export { buildDesignRequest, buildCheckRequest, isCheckablePrimer, CHECK_LIMITS, PRESETS, EXPLAIN_HINTS, revcomp,
         pairsToTSV, primersToFasta, ampliconsToFasta, offTargetsToTSV, pangenomeToTSV, estimateCheckCpu,
         initialDesignerState, normalizeDesignerState };
export { mount, ensureStylesInjected, VERSION };
export type { DesignMode, DesignRequest, DesignResponse, PrimerTemplate, PrimerPair, PrimerOligo, PrimerGenomic, PrimerWarning,
  DesignParams, CheckRequest, CheckParams, CheckJob, CheckResults, SpecificityPairResult, OffTarget, TranscriptomePairResult,
  PangenomePairResult, PangenomeGenomeResult, GenomesResponse, GenomeEntry, GrameneGene, PrimersClient, PrimersClientOptions,
  PollOptions, PrimerDesignerProps, PrimerDesignerState };
```
Types mirror §A.2 and §B.12. Readers are tolerant: unknown fields are ignored and optional fields are guarded.

**Props**
```ts
interface PrimerDesignerProps {
  apiBase: string;                                   // e.g. https://data.sorghumbase.org/sorghum_v11
  client?: PrimersClient;
  gene?: GrameneGene; geneId?: string; systemName?: string;
  region?: { region: string; start: number; end: number; strand?: 1 | -1 }; sequence?: string;
  modes?: DesignMode[]; defaultMode?: DesignMode; defaultParams?: Partial<DesignParams>;
  state?: PrimerDesignerState; onStateChange?: (s: PrimerDesignerState) => void; persistSequence?: boolean;  // default true
  onDesign?: (res: DesignResponse, req: DesignRequest) => void; onCheckUpdate?: (job: CheckJob) => void; onError?: (e: PrimersApiError) => void;
  features?: { check?: boolean; pangenome?: boolean; export?: boolean; map?: boolean };                   // default true
  geneHref?: (geneId: string, systemName?: string) => string; onGeneClick?: (geneId: string, systemName?: string) => void;
  geneLabel?: string; theme?: 'light' | 'dark' | 'auto'; className?: string; style?: React.CSSProperties;
  injectStyles?: boolean;                                                                                   // default true
  poll?: Partial<Pick<PollOptions, 'initialDelayMs' | 'maxDelayMs' | 'factor'>>;
}
```

**Serializable state** (`v: 1`, kept in the host store; responses and results are never stored)
```ts
interface PrimerDesignerState {
  v: 1; mode: DesignMode; transcriptId?: string; flankUp?: number; flankDown?: number;
  region?: { region: string; start: number; end: number; strand: 1 | -1 }; sequence?: string; systemName?: string;
  target?: [number, number]; included?: [number, number]; excluded?: [number, number][];
  junctionSpanning?: boolean; avoidRepeats?: boolean; repeatMaskMode?: 'n_mask' | 'three_prime';
  preset?: 'pcr' | 'qpcr'; params?: Partial<DesignParams>; designed?: boolean;
  selectedRank?: number; checkedRanks?: number[];
  check?: { checks: ('specificity' | 'pangenome')[]; genomes?: string[]; params?: Partial<CheckParams>;
            jobId?: string; submitted?: { id: string; left: string; right: string }[] };
  view?: { resultsTab: 'pairs' | 'specificity' | 'transcriptome' | 'pangenome'; explainOpen?: boolean };
}
```

**Client**
```ts
createPrimersClient({ apiBase, fetch?, headers?, timeouts?: { design?: 60000; other?: 15000 } }): PrimersClient
interface PrimersClient {
  design(req: DesignRequest, o?: { signal?: AbortSignal }): Promise<DesignResponse>;       // POST /primers/design
  submitCheck(req: CheckRequest, o?): Promise<CheckJob & { created: boolean }>;            // POST /primers/check (202 ⇒ created)
  getCheck(jobId: string, o?): Promise<CheckJob>;                                          // GET /primers/check/{job_id}
  pollCheck(jobId: string, o?: PollOptions): Promise<CheckJob>;                            // resolves at done | error
  runCheck(req: CheckRequest, o?: PollOptions): Promise<CheckJob>;
  listGenomes(systemName: string, o?): Promise<GenomesResponse>;                           // memoized, evicted on error
  getGene(geneId: string, o?): Promise<GrameneGene | null>;                                // GET /genes?idList=
}
interface PollOptions { signal?: AbortSignal; onUpdate?(j: CheckJob): void; initialDelayMs?: 1000; maxDelayMs?: 10000; factor?: 1.5;
                        maxConsecutiveErrors?: 5; resubmit?: CheckRequest /* only for jobs started in this mount */; pauseWhenHidden?: true }
mount(el: Element | string, props: PrimerDesignerProps): { update(p: Partial<PrimerDesignerProps>): void; unmount(): void }  // createRoot; unmount aborts requests
```
- **Transport:** `fetch` with `Accept: application/json` (plus `Content-Type` when there is a body), `cache:'no-store'`, `credentials:'omit'`, `mode:'cors'`. Path and query segments are `encodeURIComponent`-ed. No extra query params are ever added (mongo routes treat unknown params as filters).
- **`PrimersApiError {status, code, message, details, errors, retryAfterMs}`:**
  - `{message, code, details}` maps directly; the validator shape `{message, errors}` → `code: 'VALIDATION'`.
  - A non-JSON body → `HTTP_<status>`; a TypeError → `NETWORK`; a timeout → `TIMEOUT`.
  - `retryAfterMs` comes from `details.retry_after_s`, then the `Retry-After` header, then null.
  - `AbortError` is rethrown untouched.
- **Poll loop:**
  - Stop when the job is `done` or `error`.
  - On 404: resubmit once if `resubmit` is set; otherwise throw `UNKNOWN_JOB`.
  - On 503, 5xx, `NETWORK` or `TIMEOUT`: sleep `retryAfterMs ?? delay`; give up after 5 consecutive errors.
  - Delay starts at 1 s, ×1.5 per poll, capped at 10 s. It resets when `progress.done` changes and is ≥ 2 s while queued.
  - Polling pauses while the document is hidden.

### C.3 UI structure
- **Header:** ModeTabs (WAI-ARIA tabs), limited to the modes offered.
- **Left column (inputs):**
  - **gene:** transcript overlay select (canonical badge); flank up/down slider plus number input (step 50, capped so the template stays ≤ 50 kb); template-size readout; genes > 50 kb disable the mode with the hint "use Transcript or Region".
  - **transcript:** transcript select (exon count, cDNA length, CDS); junction toggle enabled when `exons.length > 1`; "Restrict to CDS" / "Clear".
  - **region:** fields prefilled from `gene.location`.
  - **sequence:** FASTA textarea with live cleaned length, invalid characters and the 20–50,000 limit; genome select (needed for checks).
  - `RepeatOptions`: avoid-repeats toggle and mode; shows the genome's `repeat_masking` and notes that BLAST-depth masking also masks gene families.
  - `IntervalsEditor`: target, included, excluded as start/length rows.
  - `ParamsPanel`:
    - PCR/qPCR presets mirroring §A.5; placeholders come from `settings.params` after a design.
    - Fields: size, Tm, GC, ΔTm, product-size chips, `num_return`; advanced fields collapsed.
    - Inline validation mirrors the server: min ≤ opt ≤ max; `a < b`; junction overlaps ≤ `floor(max_size/2)`.
    - Only changed params are sent.
  - Buttons: **Design primers** (becomes **Cancel**, which aborts the request, while running) and **Preview template** (`template_only`).
- **Right column (results):** Warnings → ErrorBanner → TemplateMap → result tabs **Pairs | Specificity | Transcriptome** (if present) **| Pan-genome** (if requested) → CheckPanel → ExportMenu.
- **TemplateMap (SVG):**
  - Tracks: ruler (tooltip shows template and genomic positions); gene-mode exons/introns/CDS or transcript-mode exon blocks with junction ticks and CDS band; repeat-mask hatch; target/included/excluded; pair lanes (greedy packing, arrows, junction notch).
  - Click or Enter selects a pair; Fit and Zoom-to-pair buttons. No drag-select in v1.
- **PairsTable** (native `<table>`):
  - Columns: checkbox (checkable only if both primers match `^[ACGTacgt]{15,36}$`, within 10 pairs / 20 unique primers); #; Left and Right 5′→3′ with copy; Tm; GC; product bp (plus a "gDNA n" line for transcripts); penalty; compl; hairpin; junction badge; after a check, a verdict chip and "amplifies n/total"; expand.
  - `PairDetail` shows stats, genomic blocks and `AmpliconSequence` (60 nt lines, footprints highlighted, junction bars, copy as FASTA).
- **CheckPanel:**
  - Genome specificity is always on.
  - The Pan-genome toggle is hidden when `listGenomes().counts.total ≤ 1`.
  - `GenomePicker`: filter, All/None; disables `!has_blastdb`, or `!has_cdna_blastdb` in transcript mode; `genomes` is omitted when everything is selected.
  - CPU estimate via `estimateCheckCpu` from `total_bases`; Submit is disabled above 6000 CPU-s.
  - Advanced check params, collapsed.
  - `buildCheckRequest`:
    ```
    { system_name, mode, gene_id?, transcript_id?, checks, genomes?, params?,
      pairs: checked.map(p => ({ id: `P${p.rank+1}`, left: p.left.seq, right: p.right.seq,
               expected: (mode === 'gene' || mode === 'region') && p.product.genomic
                         ? { region, start, end } : undefined })) }
    ```
  - Progress bar `role="progressbar"` with `aria-valuetext` (e.g. "37 of 121 · pangenome · running: sorghum_rio"); shows `queue_position` while queued.
  - **Stop watching** only detaches polling.
- **Result views:**
  - `SpecificityResults`: verdict chips (icon, colour and text); on-target line; sortable off-target table with `MismatchGlyph`; "showing 100 of N"; `gdna_products` note ("DNase-treat RNA"); sensitivity note; repetitive-primer flags.
  - `SpecificityResults` for the transcriptome block: one row per gene with an isoform list.
  - `PangenomeMatrix` (`role="grid"`, Okabe-Ito palette, glyph plus colour plus text):

    | status | glyph | colour | cell text |
    |---|---|---|---|
    | single_perfect | ✓ | #009E73 | size_delta |
    | single_mismatch | ≈ | #56B4E9 | "3′" corner flag when `mismatch_in_3p_window` or `terminal_mismatch` |
    | multiple | ×N | #E69F00 | |
    | no_amplicon | ∅ | #D55E00 | "no ortholog annotated" when `ortholog_annotated` is false |
    | db_unavailable | ? | grey hatch | |
    | error | ! | #555 | |
    | pending | … | dotted border | |

    Sort by name or worst status; "issues only" filter; roving tabindex; Enter opens CellDetail; Legend always visible.
- **ExplainPanel:** opens automatically on `NO_PAIRS`. `EXPLAIN_HINTS` is keyed by the exact 2.6.1 labels, including `too many Ns`, `no overlap of required point`, `low faction bound`, `unacceptable product size`.
- **ErrorBanner:** by code:
  - `BUSY`: auto-retry up to 3× with a countdown.
  - `QUEUE_FULL`: retry after the delay.
  - `JOB_TOO_LARGE`: shows the estimate and the limit.
  - `PRIMER3_UNAVAILABLE`, `FEATURE_DISABLED`: disable the form.
  - `NETWORK`, `HTTP_5xx`: Retry button.
  - `VALIDATION`: lists `errors[]`.
- **ExportMenu:** pairs TSV, primers FASTA (`>{label}_P{n}_F|R tm= gc= {genomic}`), amplicons FASTA, off-targets TSV, pan-genome TSV. TSV strips tabs and newlines and prefixes `'` to non-numeric cells starting with `= + - @`. Downloads use Blob with `<a download>`; copy uses `navigator.clipboard` with a textarea fallback.
- **Accessibility:** labels plus `aria-describedby`; polite and assertive live regions; table `<caption>`/`scope`; SVG `role="group"` with `<title>`; colour never the only channel; axe-core in tests.

### C.4 State and lifecycle
- **State mode:** controlled when `state` is passed; uncontrolled otherwise. Every change emits a JSON-serializable state. `persistSequence:false` drops `sequence`.
- **Identity:** `gene._id ?? geneId ?? systemName+region ?? hash(sequence)`. A change resets state and aborts both design and polling. Hosts also pass `key`.
- **useDesign:** one AbortController per request (the previous one is aborted); a monotonically increasing request id drops stale responses; `BUSY` auto-retry.
- **Restore:**
  - `designed: true` re-runs the design once.
  - With `check.jobId`, call `getCheck` once: running or queued → poll; done → render; 404 → show "Results expired — Re-run check" (no automatic resubmit).
  - Automatic resubmit on 404 applies only to jobs started in the current mount.
- **Result matching:** check results map to pairs by primer sequence (left+right), not by rank. Pairs that no longer match show "not checked".
- **Partial results:** rendered while `status: running` and `partial: true`.

### C.5 Styling
- **Delivery:** one `primers.css`, injected once as `<style id="gramene-primers-styles">` from `useInsertionEffect` in each top-level component. `injectStyles={false}` plus `gramene-primers/style.css` serves CSP-strict hosts.
- **Naming:** every class starts with `gpr-` and every variable with `--gpr-*`. No Bootstrap class names; state goes in `data-state`/`aria-*`.
- **Scoping:**
  - Everything sits under `.gpr-root`.
  - A defensive reset (`.gpr-root button|input|select|textarea|table|label|code|pre|h3|h4`, specificity 0,1,1) beats Bootstrap 4 reboot, e.g. `code{color:#e83e8c}`.
  - No `!important`.
- **Themes and layout:** `.gpr-theme-light`, `.gpr-theme-dark`, and `auto` via `prefers-color-scheme`. `container-type: inline-size` with a two-column layout at ≥ 960px. `prefers-reduced-motion` respected.
- **Colours:** feature colours copied from gramene-search `sequences.css` (utr5 #aaccaf, cds #a7b4d3, utr3 #c5a3bf, intron #aaaaaa), under `gpr-` names.

### C.6 Tests (vitest + jsdom + Testing Library)
- **Unit:**
  - `client.test.ts`: URL, headers, `no-store`, no extra params, 202 vs 200, error mapping incl. validator shape, HTML 502, network, timeout, `retry_after_s` body and header, header-less 503.
  - `poll.test.ts` (fake timers): backoff 1000 → 1500 → 2250 … capped at 10000; reset on progress; ≥ 2000 while queued; terminal states; `maxConsecutiveErrors`; 404 resubmit only with `resubmit`; visibility pause; abort.
  - `request.test.ts`: per-mode requests; changed params only; `expected` only for gene/region; genomes omitted when all selected; `isCheckablePrimer`.
  - `state.test.ts`, `exporters.test.ts`, `explain.test.ts`, `pangenome.test.ts`, `coords.test.ts` (verified pair §A.4.4; `MismatchGlyph` index = len − pos), `cost.test.ts` (matches the server coefficients).
- **Components** (fake `PrimersClient`):
  - Gene: SORBI_3001G000200; the gene doc is not mutated.
  - Transcript: SORBI_3001G046200, canonical `.2` single-exon → junction toggle disabled; `.1` enables it.
  - No pairs → ExplainPanel.
  - Errors: VALIDATION, BUSY countdown, FEATURE_DISABLED.
  - Reset on gene change aborts requests.
  - Check flow queued → running (partial) → done: verdict chips, matrix keyboard navigation, CellDetail.
  - Pan-genome toggle hidden when `counts.total` is 1.
  - Controlled restore: expired job shows "Re-run".
  - axe-core on the main states; `mount` / `update` / `unmount`.
- **Contract fixtures:**
  - `npm run fixtures` writes `buildDesignRequest` / `buildCheckRequest` output for every mode to `test/fixtures/contract/requests/*.json`; copy them into the swagger `test/primers/fixtures/contract/requests/`.
  - `scripts/capture-fixtures.mjs` (run against `:50111`) stores real design, genomes and finished-check responses in `test/fixtures/api/`, used by the component tests.
- **Integration** (opt-in, `PRIMERS_IT_BASE`): real client against the dev server; guards pass; `product_size == right.end - left.start + 1`.
- **Package checks:** `npm run lint:pkg` is clean; `grep -c 'react.production' dist/gramene-primers.js` returns 0.

### C.7 Playground (local only in v1)
- Pages:
  - Gene: SORBI_3001G000200 (−), SORBI_3001G000700 (+)
  - Transcript: SORBI_3004G087700.3, SORBI_3001G046200
  - Region, Sequence
  - Check: P1/P2/P3/P5_L; J_L+P1_R
  - Theme toggle and a controlled-state JSON inspector
- `?api=mock` replays fixtures; `?api=live` uses `/sorghum_v11` through the proxy: `PRIMERS_PROXY_TARGET=http://localhost:50111 npx vite --port 5174`.

---

## 8. Part D: Primers tab in gramene-search

**Base.** Fresh clone of GitHub master, not the stale submodule at `/usr/local/gramene/release69/gramene-meta/modules/gramene-search` (255 commits behind, React 17):
```bash
git clone git@github.com:warelab/gramene-search.git /home/olson/src/warelab/gramene-search
cd /home/olson/src/warelab/gramene-search && git checkout -b primers-tab && npm install
```

**Changes** (line numbers refer to master 5d135c3)

| File | Change |
|---|---|
| `package.json` | `"gramene-primers": "^1.0.0"` in dependencies; during development install the tarball with `--no-save` |
| `src/components/results/details/Primers.js` | New file (below) |
| `src/components/results/GeneList.js` | After line 12: `import Primers from "./details/Primers"`. In `inventory` (line 21): `primers: Primers,`. In `allDetails`, after the `sequences` entry (~line 146): `{ id: 'primers', label: 'Primers', popup: 'Design PCR/qPCR primers; check specificity and pan-genome coverage', available: true }`. The existing filter at line 196 already treats `true` or an object as enabled |
| `src/bundles/uiViewState.js` | Header comment: `primers: PrimerDesignerState`. Reducer case `'UI_PRIMERS_SET': return setGene(state, payload.geneId, { primers: payload.state })`, next to `UI_EXPRESSION_SET`. Action `doSetPrimersState: ({ geneId, state }) => ({ dispatch }) => dispatch({ type: 'UI_PRIMERS_SET', payload: { geneId, state } })`, next to `doSetExpressionState` |
| `src/bundles/viewSnapshot.js` | Mirror `sequences` at lines 34 (schema comment), 92 (`if (entry.primers && Object.keys(entry.primers).length) return true;`), 152 (`primers: e.primers ? {...e.primers} : {}`), 235 (`primers: e.primers && Object.keys(e.primers).length ? {...e.primers} : undefined`) |
| `src/demo.js` | Sorghum `details`: `primers: process.env.PRIMERS_API ? { apiBase: process.env.PRIMERS_API } : true` |
| `README.md` | Document `config.details.primers` (`true \| {apiBase, pangenome}`) and tarball linking |

```jsx
// src/components/results/details/Primers.js
import React, { useCallback } from 'react'
import { connect } from 'redux-bundler-react'
import { PrimerDesigner } from 'gramene-primers'

const Primers = props => {
  const geneId = props.searchResult.id
  const gene = props.geneDocs[geneId]                       // Gene.ensureGene() guarantees it is loaded
  const conf = props.config.details.primers                 // true | { apiBase?, pangenome? }
  const apiBase = (conf && typeof conf === 'object' && conf.apiBase) || props.grameneAPI
  const byGene = props.uiViewState && props.uiViewState.byGene && props.uiViewState.byGene[geneId]
  const { doSetPrimersState } = props
  const onStateChange = useCallback(state => doSetPrimersState({ geneId, state }), [geneId, doSetPrimersState])
  return (
    <PrimerDesigner key={geneId} apiBase={apiBase} gene={gene} systemName={gene.system_name}
      modes={['gene', 'transcript', 'region', 'sequence']} defaultMode="gene"
      state={byGene && byGene.primers} onStateChange={onStateChange} persistSequence={false}
      features={{ pangenome: !(conf && conf.pangenome === false) }}
      geneHref={id => `?idList=${encodeURIComponent(id)}`} geneLabel={gene.name || gene._id} />
  )
}
export default connect('selectGrameneAPI', 'selectUiViewState', 'doSetPrimersState', Primers)
```
- No CSS import (styles are injected at runtime) and a static import (no lazy chunk, which avoids `--public-url` issues).

**Dev linking** (tarball, not `npm link` or `file:`, which would create a second React and trigger `source` resolution)
```bash
cd /home/olson/src/warelab/gramene-primers && npm run pack:local            # gramene-primers-0.1.0.tgz
cd /home/olson/src/warelab/gramene-search && npm install --no-save ../gramene-primers/gramene-primers-0.1.0.tgz && rm -rf .parcel-cache
PRIMERS_API=http://localhost:50111/sorghum_v11 SUBSITE=sorghum npx parcel src/sorghum.html --port 1234
```

**Release order**
1. Merge the swagger branch and deploy the API (§11).
2. `npm publish gramene-primers@1.0.0`.
3. Publish gramene-search 2.24.0.
4. In `warelab/sorghum-webapp` `search_app`: bump `gramene-search`, add `primers: true` to `config.json` `details`, run `npm run build`.
5. gramene-sites and other subsites stay off until their swagger has `/primers`.

---

## 9. Reuse

| Existing code | Path | Reused for |
|---|---|---|
| `util.is*` polyfill | `gramene-swagger/app.js` lines 3–25 | `node_compat.js` (forked worker, tests) |
| Mongo collections, `getMongoConfig()` | `../gramene-mongodb-config/collections.js` | genes, maps, taxonomy access; `site_key` |
| async controller style | `api/controllers/gxa.js` | `primers.js` |
| MAX_IDS guard style | `api/controllers/solr.js` | Limit checks |
| Validator and `json_error_handler` bagpipes | `config/default.yaml`, `node_modules/swagger-node-runner/fittings/` | Error shape; `sway` in `contract.test.js` |
| Docs style | `docs/gene_lists_api.md` | `docs/primer_design_api.md` |
| `BgzipIndexedFasta` usage | `/home/olson/src/warelab/fastaIdx/api.js` | `sequence.js` (IUPAC revcomp replaces its uppercase-only one) |
| spawn argument arrays, job hashing | `/usr/local/gramene/gramene-blast/api.js` | `blast.js`, `primer3.js` |
| Gene-structure semantics | `../gramene-mongodb/search/dump_genes.js` (`add_exons_to_transcripts`, `add_translations_to_transcripts`) | Exon, CDS and junction conventions |
| Solr capability list (not changed; `available:true` used) | `../gramene-solr/genes/mongo2solr.join.js` | Confirms no reindex is needed |
| BLAST DBs, `cdna.all` DBs, bgzip FASTA with `.fai`/`.gzi` | `/scratch/olson/fasta/<system_name>/` | Resolver inputs |
| Detail-tab registry, `ensureGene`, `FullscreenContainer` | gramene-search `GeneList.js`, `details/generic.js` | Tab wiring |
| `connect('selectGrameneAPI')` | gramene-search `details/Homology.js` (line 632) | Tab API base |
| `sequences` view-state slice | gramene-search `bundles/uiViewState.js`, `bundles/viewSnapshot.js` | `primers` slice |
| `decorateSeq` colours | gramene-search `details/Sequences.js`, `sequences.css` | Palette only (copied, not imported) |
| Packaging and runtime style injection | `github.com/warelab/tbrowse` (`package.json`, `vite.config.ts`) | gramene-primers build |
| `should`, `supertest` | swagger devDependencies | Tests |

---

## 10. Verification

### 10.0 Environment guards (before and after each step)
```bash
git -C /usr/local/gramene/subsites/sorghum/v11/gramene-swagger status --porcelain   # only: ?? package-lock.json
pm2 describe sorghum_swagger11 | grep -E 'pid|uptime'                               # unchanged
```

### 10.1 Primer3
- `/home/olson/bin/primer3_core -about` prints 2.6.1.
- The 50 kb benchmark (§0.1) prints `20` in under 10 s.

### 10.2 Swagger unit and contract tests
```bash
cd /usr/local/gramene/subsites/sorghum/v11/gramene-swagger-primers && npm run test:primers     # all pass, incl. contract.test.js
```

### 10.3 Design and genomes against the dev server
```bash
BASE=http://localhost:50111/sorghum_v11
post() { curl -s -w '\n%{http_code}\n' -H 'Content-Type: application/json' -d "$2" "$BASE$1"; }
```

| # | Request | Expected |
|---|---|---|
| 1 | `post /primers/design '{"mode":"transcript","gene_id":"SORBI_3001G000200"}'` | 200; `template.length` 1982; `features.junctions[:3]` `[397,493,597]`; ≤ 5 pairs; every pair has a primer with non-null `junction`; `product_size` 70–150; `product_tm` numeric |
| 2 | `post /primers/design '{"mode":"gene","gene_id":"SORBI_3001G000200","flank_up":200,"flank_down":100}'` | `template` `1:11080-15099`, strand −1, length 4020; `features.exons[0]` 201–597; `features.cds.start` 499 |
| 3 | Coordinate cross-check through fastaIdx | `for r in 1:13637..13650:-1 1:13395..13400:-1; do curl -s localhost:8888/sequence/region/sorghum_bicolor/$r \| python3 -c 'import json,sys;print(json.load(sys.stdin)["seq"],end="")'; done` → `ATTACATCAAATAGGCCTTG`; `…/1:13291..13310:1` → `AACTTCTTTGTCGATCCATG` (verified). The integration test does this for every returned primer |
| 4 | `post /primers/design '{"mode":"gene","gene_id":"SORBI_3001G000700"}'` (+ strand) | 200; genomic blocks verify through fastaIdx |
| 5 | `post /primers/design '{"mode":"sequence","sequence":">x\nACGTRYACGT…(≥20 nt)"}'` | 200 with warning `IUPAC_CONVERTED` (not `PRIMER3_INPUT_ERROR`) |
| 6 | `curl -s -d '{}' $BASE/primers/design` (no Content-Type) | 400; `message` "Validation errors"; `errors[].code` includes `INVALID_CONTENT_TYPE` |
| 7 | `post /primers/design '{"mode":"gene","gene_id":"X","bogus":1}'` | 400 Validation errors |
| 8 | `post /primers/design '{"mode":"gene","gene_id":"NOPE"}'` | 404 `UNKNOWN_GENE` |
| 9 | `post /primers/design '{"mode":"region","system_name":"../etc","region":{"region":"1","start":1,"end":100}}'` | 400 Validation errors (pattern) |
| 10 | `post /primers/design '{"mode":"transcript","gene_id":"SORBI_3001G000200","params":{"max_size":24,"min_5_prime_overlap_of_junction":13}}'` | 400 `INVALID_PARAMS` |
| 11 | `post /primers/design '{"mode":"gene","gene_id":"SORBI_3004G087700","avoid_repeats":true,"template_only":true}'` | `mask_source` `blast_depth`; warning `BLAST_DEPTH_MASK`; `pairs: []` |
| 12 | `post /primers/design '{"mode":"region","system_name":"sorghum_tx436pac","region":{"region":"4","start":7547610,"end":7564601},"avoid_repeats":true}'` | `mask_source` `softmask`; `masked_fraction > 0`; no primer overlaps `mask` |
| 13 | `curl -s "$BASE/primers/genomes?system_name=sorghum_bicolor"` | `counts.total` 120; `genomes[0].system_name` sorghum_bicolor; sorghum_rio `repeat_masking` `soft_masked`; no filesystem paths |
| 14 | `curl -s -X PUT $BASE/primers/check/0123456789abcdef0123456789abcdef` | JSON 405 `METHOD_NOT_ALLOWED`, not an HTML stack trace |
| 15 | `curl -sI -X POST -H 'Content-Type: application/json' -d '{"mode":"transcript","gene_id":"SORBI_3001G000200"}' $BASE/primers/design` | `Cache-Control: no-store` |

### 10.4 Checks (fixtures from the specificity analysis, sorghum_bicolor reference)

| Case | Body (abridged) | Expected |
|---|---|---|
| P2 | `mode gene, gene_id SORBI_3004G087700, pairs [{id P2, left GGACAGCTCCACAACATATCAG, right GGACATTTGAAGCCCATGGCC, expected 4:7423537-7423746}]` | First POST 202, second POST 200 with the same `job_id`. `verdict off_targets`; 4:7437317–7437526 LR 0/0 (SORBI_3004G087800); 5:66890615–66890824 RL, L mm `[10,7]`, R mm `[15,10]`, `likely` |
| P1 | L `GATCGACAATCCGACGATAGAAG`, R `GTGAACATCATGCTGCCCGATG`, expected 4:7422190–7422843 | off-targets 4:7435938–7436628 (691 bp) and 5:66891530–66892235 (706 bp); P1_L site at 4:7432183 ignored (7 mm) |
| P3 | L `GATATCAGTGGAATCATAAGACCG`, R `CATCGATATCAGGATCTGGCTT`, expected 4:7422482–7423061 | off-target 4:7436221–7436845 `likely_weak` |
| P5_L repetitive | pair P5_L `CGACGAGGAGGCGCTGCGATG` + P2_R | `primers["CGACGAGGAGGCGCTGCGATG"].repetitive` true |
| Sensitivity | pair `m20_7_14` `GGACAGATCCACATCATATC` + P2_R, no expected | Site at 4:7423537 with left_mm 2 is present (found at ws5; would be missed at ws7) |
| Indel | pair `del18` `GGACAGCTCCACAACATTCAG` + P2_R | Amplicon at 4:7423537 with left_mm 1 (gap), `likely`, not dropped |
| qPCR | `mode transcript, gene_id SORBI_3004G087700`, pair J_L `CCAACAAAGTCATGGATGCACT` + P1_R | genome specificity verdict ≠ `on_target_missing`; `transcriptome.on_target.isoforms` = .1/.2/.3 (278 bp); `transcriptome.off_targets` = 2 gene groups, SORBI_3005G183900 (isoforms .1, .2) and SORBI_3004G087800 |
| Pan-genome | P3, `checks [specificity,pangenome]`, genomes `[sorghum_353,sorghum_grassl,sorghum_leoti]` | 353 and grassl primary 580 bp, `terminal_mismatch` true, `likely_weak`; leoti primary 581 bp, right_mm 1 |
| Presence/absence | L `TTACCTCTCCTCCTTCCCGATC`, R `GCATAGTCTCATTAGCTAGCAC`, gene SORBI_3001G046200, pangenome over `[sorghum_353,sorghum_is12661,sorghum_ji2731,sorghum_leoti]` | 353 and is12661 `no_amplicon` with `ortholog_annotated` false; ji2731 `size_delta` −177; leoti 0 mm |
| Size variation | L `CCCTCCAACGTTGAAGGGATTG`, R `CTCATCTCATTACTTTCCACCC`, gene SORBI_3001G088450 (− strand) | reference 1:6877536–6879899 (2364 bp, RL); 353 `size_delta` +22; s3691 −90 |
| Contract | the §A.2.3 example verbatim | 202/200; a legacy body with top-level `max_mismatches` → 400 Validation errors |
| Cost guard | 10 pairs × 20 primers, pangenome over all genomes plus a synthetic large list | `estimate.cpu_s` reported; > 6000 → 422 `JOB_TOO_LARGE` |

**Operational checks during the first full-panel run**
- `vmstat 5` for I/O wait.
- `redis-cli -p 6380 INFO memory` (used memory delta).
- `redis-cli -p 6380 -n 1 --scan --pattern 'primers:sorghum_v11_dev:*'`, then `ttl` on each key > 0.
- Wall time within §B.13 estimates ±50%; `timings_ms` recorded in the docs.

**Restart recovery:** Ctrl-C the dev server during a pan-genome job, restart it; the job goes back to `queued` and finishes `done` with `attempts: 2`.

**Redis down:** restart dev with `PRIMERS_REDIS_URL=redis://localhost:6399`.
- `curl $BASE/genes?idList=SORBI_3001G000200` still answers.
- `POST /primers/check` → 503 `JOB_STORE_UNAVAILABLE` within 3 s.
- `/primers/design` still works.

**Unhandled rejection:** with `PRIMERS_FASTA_ROOT=/nonexistent`, design returns JSON errors and the process stays up.

### 10.5 gramene-primers
- `npm run typecheck && npm test && npm run build && npm run lint:pkg`: all green.
- `grep -c 'react.production' dist/gramene-primers.js` returns 0.
- `npm run test:it` with `PRIMERS_IT_BASE` passes.
- Playground `?api=live` designs SORBI_3001G000200 in all four modes.

### 10.6 gramene-search / sorghum-webapp manual checklist (Parcel dev on :1234 via SSH tunnel)
1. Tab "Primers" appears after "Sequences" only when `details.primers` is set.
2. Gene, transcript (junction badge), region and pasted-FASTA sequence designs work.
3. Map selection syncs with the table.
4. The specificity check shows verdict chips.
5. Pan-genome matrix renders with keyboard navigation.
6. Saved-view restore: design re-runs; an expired job shows "Re-run check"; there is no automatic resubmit.
7. No console errors; Sequences tab styling unaffected; primers UI unaffected by Bootstrap 4 reboot in a local sorghum-webapp `search_app` build using the tarball.

---

## 11. Risks and follow-ups

**Risks**
- Primer3 runtime on 20–50 kb templates is only gated by the Step 0 benchmark; `max_template_length` may drop to 20 kb.
- **BLAST sensitivity** is bounded and documented:
  - ws5 guarantees ≤ 3 substitutions for 20–22-mers; ws6 guarantees ≤ 2.
  - Mid-primer indels in short primers can be missed.
  - Pan-genome costs about 2× ws7. Coefficients (5.2 / 2.2 / 1.2 CPU-s per primer·Gb) were measured on sorghum and must be re-measured for release70 genomes such as barley.
- The BLAST-depth mask has recall 0.60 and masks multi-copy gene families; `n_mask` can leave no primers there (MOSTLY_REPEAT / explain hints).
- The Primer-BLAST rules are heuristics. Paralog-rich genes will often show `multiple`.
- Redis 6380 is shared with gramene-blast (3.3 GB, noeviction). It is bounded by gzip, TTLs, the 500-finished cap and 5 s progress throttling.
- The host CPU cap ignores gramene-blast and solr/mongo I/O. One pan-genome job at a time plus nice 10 mitigates this; watch `vmstat`.
- The forked worker lives under the pm2 app. A hard kill leaves at most one stale-slot window (60 s).
- gramene-mongodb-config never retries a failed connection. Mitigated by worker exit 75; the API returns 503.
- **Merging into release70:**
  - Its `app.js` has an uncommitted PORT edit (`50011 → 50070`) plus untracked `app.js.bak.port50011` and `package-lock.json`, so git will refuse the merge until the owner moves PORT into the pm2 env and reverts or commits.
  - Its basePath is still `/sorghum_v11`; the site key includes the mongo db (`search70`).
- A caching Apache front end sits in front of the API; `no-store` must survive the proxy (verify with `curl -D -` on data.sorghumbase.org after deploy).
- Coordinated releases: gramene-primers → gramene-search → sorghum-webapp. The tab stays hidden until the API is live.
- A saved view re-runs Primer3; if presets change, pairs may differ from stored check pairs (handled by sequence matching).
- No auth or per-IP limits. Global caps protect the host but not fairness between clients.

**Deploy (performed by the user, later)**
1. Merge `primer-design` into `sorghum_v11`.
2. Build Primer3.
3. `npm install` in the live checkout (lockfile moved aside first).
4. `pm2 restart sorghum_swagger11`.
5. Verify §10.3 against `https://data.sorghumbase.org/sorghum_v11`.
6. Then release70.

**Deferred from v1**
- IIFE/UMD standalone bundle and GitHub Pages playground.
- Map drag-select and zoom popover.
- CSV, order-sheet and JSON-bundle exports.
- Disk site cache for BLAST sites; ops CLI.
- `DELETE /primers/check/{id}` and `cancelled` status.
- `expected_transcript` coordinates.
- Genome-DB pass for transcript-mode pan-genome.
- gramene-mongodb-config reconnect fix, as a separate reviewed change.
- `PRIMER_MASK_TEMPLATE` k-mer lists.
- Real soft-masking of the 80 unmasked sorghum dna_sm copies.
- bgzip and index of `cdna.all` FASTA (exact cDNA re-alignment) and of the 13 masked-but-unindexed dna_sm files.
- React 19 peer range.
- Enabling other subsites (release70 barley BLAST DB naming mismatch: only 2 of 77 have local DBs).
- `trust proxy` and rate limiting.
