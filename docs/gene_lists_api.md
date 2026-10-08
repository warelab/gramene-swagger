# Gene Lists API

User "gene lists" (aka saved searches) let a signed-in user save, share, and re-open a set of genes.
This documents the current `/gene_lists*` endpoints served by this gramene-swagger instance.

- **Base path:** `https://data.sorghumbase.org/sorghum_v11` (local dev: `http://localhost:50011/sorghum_v11`)
- **Controllers:** `api/controllers/mongo.js` (list CRUD + save) and `api/controllers/solr.js` (`validate`)
- **Storage:** the `genelists` mongo collection lives in the shared **`userData1`** database (not the
  release-specific db), so a user's lists persist across site releases. The `site` field scopes a list to a
  site (e.g. `sorghum`). Configured in `gramene-mongodb-config/collections.js`.

> **Clients written against the pre-2026-08 API must change.** `validate` no longer returns a hash and no
> longer saves anything, and saving moved from query params to a JSON body. See
> [Migrating an existing client](#migrating-an-existing-client) — one of the changes fails **silently**.

## Authentication

Firebase **ID token** passed as `Authorization: Bearer <token>` (verified server-side via `firebase-admin`).

| Endpoint | Auth |
| --- | --- |
| `POST /gene_lists/validate` | none |
| `GET /gene_lists` | optional (anonymous sees only public **active** lists; `uid` defaults to 0) |
| `POST /gene_lists` (save) | **required** |
| `PATCH /gene_lists` (update) | required |
| `DELETE /gene_lists` (delete) | required |
| `POST /gene_lists/restore` | required |

> **All POST requests must send `Content-Type: application/json`**, otherwise swagger rejects them with
> `400 INVALID_CONTENT_TYPE` before your handler is reached.

## Data model — a `genelists` document

```jsonc
{
  "_id":       "1435321525 MHudLfNAnWTQakT4ZSGG07CELuq1",  // composite key: "<hash> <uid>"
  "hash":      1435321525,        // murmur3-x86-32, computed BY THE SERVER (see POST /gene_lists)
  "uid":       "<firebase uid>",
  "owner":     "Sunita Kumari",   // decoded token name / email / uid
  "site":      "sorghum",
  "label":     "Coexpressed genes in shoot",
  "n_genes":   32,                // derived from the posted ids, never taken from the client
  "isPublic":  true,
  "createdAt": "2026-06-05T21:05:41.724Z",
  "deletedAt": "2026-07-06T17:00:00.000Z"  // present ONLY when soft-deleted; absent = active
}
```

The member **gene IDs are not stored** on the document. The `hash` is the link to the genes Solr core: a
list's genes are the docs matching `saved_search:<hash>` in `sorghum_genes11`.

---

## The two-step save flow

Validation and saving are now **separate calls with no shared state**. Nothing is persisted until step 2,
so a user who abandons the dialog leaves no trace.

```
1. POST /gene_lists/validate   (no auth)  -> resolved / ambiguous / unknown
2.  ...client resolves ambiguity, user confirms...
3. POST /gene_lists            (auth)     -> server computes the hash, saves, tags the core
```

---

## `POST /gene_lists/validate`

Resolve submitted identifiers against the genes core. **Read-only** — no hash, no writes, no side effects
of any kind. Safe to call on every keystroke or paste.

- **Auth:** none. **Body:** a JSON array of identifier strings.

```
POST /gene_lists/validate
Content-Type: application/json

["SORBI_3001G000200", "sobic.001g000200", "LOC_Os09g01930", "NOT_A_GENE"]
```

```jsonc
{
  "resolved":  [ { "input": "SORBI_3001G000200", "id": "SORBI_3001G000200" },
                 { "input": "sobic.001g000200",  "id": "SORBI_3001G000200" } ],
  "ambiguous": [ { "input": "LOC_Os09g01930", "matches": ["Os09g0106200", "Os09g0354900"] } ],
  "unknown":   [ "NOT_A_GENE" ]
}
```

**Matching.** Each input is checked against both the stable `id` field and the `alt_id` field, so
`Sobic.*`, `LOC_Os*`, `GRMZM*` and similar alternate identifiers now resolve (previously they all came back
as missing). Every input you send appears in exactly one of the three arrays; order within each array
follows the submitted order.

**Ambiguity is the client's job.** Any input matching more than one gene is reported as `ambiguous`,
*including* when one of the matches is an exact stable id — the server deliberately does not pick a winner.
This is common: `loc_os07g23485` matches 15 genes, `grmzm2g034428` matches 10. Present the `matches` to the
user and send back only the chosen ids in step 2.

**Case sensitivity is asymmetric, by index design.** `id` is case-sensitive; `alt_id` is indexed
lowercased. One real consequence: `Sobic.001G000200` comes back **ambiguous** — a gene in
`sorghum_bicolort2tcas` literally has that stable id, while `sorghum_bicolor`'s `SORBI_3001G000200` carries
it as an alt_id — but the lowercase `sobic.001g000200` resolves uniquely to `SORBI_3001G000200`, because
lowercase misses the case-sensitive `id` match. Not a bug; don't "fix" it by lowercasing user input.

**Whitespace and delimiters are the client's responsibility.** The server only drops empty strings and
duplicates. Split and trim before posting.

**Limits** (also declared in swagger, so violations are rejected before the handler runs):

| rule | response |
| --- | --- |
| more than **6000** ids | `400` |
| any id longer than **255** chars | `400` |
| a non-string array entry | `400` |
| body is not a JSON array | `400` |
| solr unreachable | `502` |

Injection is not a concern — the lookup uses Solr's `{!terms}` parser, which treats values literally.
`"SORBI_3001G000200) OR id:(*"` comes back as `unknown`.

---

## `POST /gene_lists` — save a list

Takes the **resolved** identifiers in a JSON body. The server computes the hash, derives `n_genes`, upserts
the document, and tags the member genes in the Solr core.

- **Auth:** required.

```
POST /gene_lists
Content-Type: application/json
Authorization: Bearer <firebase id token>

{
  "label":    "my list",
  "site":     "sorghum",
  "isPublic": false,
  "ids":      ["SORBI_3001G000200", "SORBI_3001G000300"]
}
```

| field | type | required | notes |
| --- | --- | --- | --- |
| `label` | string | yes | non-empty after trim; max 255 |
| `site` | string | yes | e.g. `sorghum`; max 64 |
| `ids` | string[] | yes | non-empty; max 6000; each max 255 chars |
| `isPublic` | boolean | no | defaults to `false` (anything other than literal `true` is false) |

**Do not send `hash` or `n_genes`.** Both are server-derived. A body containing `hash` is rejected with
`400` rather than silently ignored, so an un-migrated client fails loudly.

```jsonc
{ "message": "list saved", "hash": 1435321525,
  "_id": "1435321525 MHudLfNAnWTQakT4ZSGG07CELuq1",
  "n_genes": 32, "upserted": true }
```

- `upserted` is `true` on first save, `false` when it updated an existing doc.
- **Idempotent.** `_id` is `"<hash> <uid>"` and the hash is content-derived, so re-saving the same id set as
  the same user updates in place rather than creating a duplicate. Solr tagging uses `add-distinct`.
- **Order-independent:** ids are deduplicated and sorted before hashing, so `[A,B]` and `[B,A]` produce the
  same list.
- **Side effect:** adds `hash` to `saved_search` on every id, awaited and committed before the response —
  so the list is queryable the moment you get a `200`. (It used to be fired off after responding, which
  meant an immediate read could miss it.)

| failure | response |
| --- | --- |
| no / malformed `Authorization` header | `401 Authorization header missing or malformed` |
| token rejected | `401 Authorization failed` |
| body contains `hash` | `400 hash is computed by the server; remove it and post 'ids' instead` |
| missing `label` / `site` / empty `ids` | `400` with the field named |
| >6000 ids, or an id >255 chars | `400` |
| mongo or solr write failed | `500 Failed to save gene list` |

---

## `GET /gene_lists`

List saved lists.

| param | in | type | required | notes |
| --- | --- | --- | --- | --- |
| `site` | query | string | yes | e.g. `sorghum` |
| `isPublic` | query | boolean | yes | `true` = **public** lists for the site; `false` = the **caller's own** lists (needs token) |
| `rows` | query | integer | no | default 20; `-1` also yields the default 20 |
| `includeDeleted` | query | string | no | `active` (default) or `trash` |

- `trash` returns **the caller's own soft-deleted lists**; `401` if anonymous.
- **Response:** a streamed JSON array of `genelists` documents. Note it is streamed via `JSONStream`, so
  the raw text has newlines and leading commas between elements — fine for `JSON.parse`/`response.json()`,
  but don't hand-parse it line by line.

---

## `PATCH /gene_lists` — update a list

- **Auth:** required. Query `listId` (the `_id`, required) + JSON body `{ "label"?: string, "isPublic"?: boolean }`.
- Only those two fields can change; a list's membership is immutable (save a new list instead — a different
  id set produces a different hash and therefore a different `_id`).
- **Response:** `{ "message": "list updated", "updated": {...} }`; `400` if no valid fields; `404` if not owned/found.

---

## `DELETE /gene_lists` — delete a list

| param | in | type | required | notes |
| --- | --- | --- | --- | --- |
| `listId` | query | string | yes | the `_id` |
| `force` | query | boolean | no | `true` = permanent now; omit/false = 30-day soft delete |

- **Default (soft):** sets `deletedAt`; hidden from the active view but **restorable for 30 days**, then
  purged by a daily cron. → `{ "message": "list marked for deletion (restorable for 30 days)" }`
- **`force=true`:** permanent. → `{ "message": "list permanently deleted" }`
- `404` if not found / not owned (soft delete also `404`s if already soft-deleted).
- **Never writes Solr `saved_search`.** A deleted list's genes keep their tag; cross-release propagation and
  cleanup are handled separately by `gramene-solr/scripts/sync_saved_search.js`.

---

## `POST /gene_lists/restore` — restore a soft-deleted list

- **Auth:** required. `Content-Type: application/json`. Query `listId` (required).
- Clears `deletedAt` (only within the 30-day window before purge).
- **Response:** `{ "message": "list restored" }`; `404` if nothing to restore.

---

## Opening a list — fetching its genes

A list's genes are the ones tagged with its hash. **Query `/search`, not `/genes`:**

```
GET /search?q=saved_search:1435321525&rows=100&fl=id,name
```

`/search` is the Solr-backed endpoint and understands `saved_search`. `/genes` is a **different, Mongo-backed
endpoint** that takes `q` / `idList` and knows nothing about `saved_search` — pointing it at
`saved_search:<hash>` returns an empty array with a `200`, no error. (Earlier revisions of this document said
either would work; that was wrong.)

To hydrate full gene documents, take the ids from `/search` and pass them to Mongo:

```
GET /genes?idList=Sobic.001G066900,Sobic.001G189500
```

---

## Migrating an existing client

### 1. `validate` no longer returns `hash`, and no longer saves

| before | now |
| --- | --- |
| `{ ids: [...], missing: [...], hash: <n> }` | `{ resolved: [{input,id}], ambiguous: [{input,matches}], unknown: [...] }` |
| tagged `saved_search` on the core as a side effect | writes nothing |
| stable ids only | stable ids **and** alternate ids |
| no way to express a multi-match | `ambiguous` array |

Code reading `response.hash` gets `undefined`; code reading `response.ids` / `response.missing` gets
`undefined` too. Both fail fast and visibly.

New UI work required: **there was previously no way to surface ambiguity**, so a client that just concatenates
`resolved` and ignores `ambiguous` will silently drop genes the user asked for. At minimum, show the
ambiguous inputs with their candidate ids and let the user choose.

### 2. Saving moved from query params to a JSON body

```diff
- POST /gene_lists?label=x&hash=123&site=sorghum&n_genes=9&isPublic=true
+ POST /gene_lists
+ Content-Type: application/json
+ { "label": "x", "site": "sorghum", "isPublic": true, "ids": [ ...resolved ids... ] }
```

### 3. ⚠️ The old save call now fails silently

`GET /gene_lists?label=x&hash=123&site=sorghum&n_genes=9&isPublic=true` still returns **`200`** — but it is
now just the *listing* endpoint ignoring unknown query params. **Nothing is saved.** Verified: the
`genelists` document count is unchanged before and after such a call.

If any client still saves via `GET`, it will appear to work and quietly lose every list. Grep the client for
`gene_lists?` with `hash=` and convert those call sites first. The `POST` form of the old query-param call is
safe by comparison — it fails loudly with `400 INVALID_CONTENT_TYPE`.

### 4. Hashes changed algorithmically (existing lists unaffected)

The hash is now `murmur3(ids.join(','))` rather than `murmur3(ids.join(''))`; the old separator-less join
collided (`["AB","C"]` and `["A","BC"]` hashed identically). Stored hashes are never recomputed, so existing
saved lists keep working. But **do not compute hashes client-side** and expect them to match — the server is
the only place that derives them.

---

## Lifecycle summary

1. **Validate:** `POST /gene_lists/validate` → resolve `ambiguous` with the user. Nothing is stored.
2. **Save:** `POST /gene_lists` with `{label, site, isPublic, ids}` → server returns the `hash`.
3. **Open:** `GET /search?q=saved_search:<hash>`; hydrate via `GET /genes?idList=…`.
4. **Delete:** `DELETE /gene_lists?listId=…` → soft (30-day trash); `&force=true` for permanent.
5. **Trash + restore:** `GET /gene_lists?…&includeDeleted=trash`, then `POST /gene_lists/restore?listId=…`.
6. **Auto-purge:** lists soft-deleted > 30 days are hard-deleted daily by
   `scripts/cleanup_expired_genelists.js` (olson cron; mongo-only).
