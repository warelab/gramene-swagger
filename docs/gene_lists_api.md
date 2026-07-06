# Gene Lists API

User "gene lists" (aka saved searches) let a signed-in user save, share, and re-open a set of genes.
This documents the current `/gene_lists*` endpoints served by this gramene-swagger instance.

- **Base path:** `https://data.sorghumbase.org/sorghum_v11` (local dev: `http://localhost:50011/sorghum_v11`)
- **Controllers:** `api/controllers/mongo.js` (list CRUD) and `api/controllers/solr.js` (`validate`)
- **Storage:** the `genelists` mongo collection lives in the shared **`userData1`** database (not the
  release-specific db), so a user's lists persist across site releases. The `site` field scopes a list to a
  site (e.g. `sorghum`). Configured in `gramene-mongodb-config/collections.js`.

## Authentication

Firebase **ID token** passed as `Authorization: Bearer <token>` (verified server-side via `firebase-admin`).

| Endpoint | Auth |
| --- | --- |
| `POST /gene_lists/validate` | none |
| `GET /gene_lists` | optional (anonymous sees only public **active** lists; `uid` defaults to 0) |
| `POST /gene_lists` (save) | required |
| `PATCH /gene_lists` (update) | required |
| `DELETE /gene_lists` (delete) | required |
| `POST /gene_lists/restore` | required |

> **All POST requests must send `Content-Type: application/json`** (even when the payload is entirely in
> query params), otherwise swagger rejects them with `400 INVALID_CONTENT_TYPE`.

## Data model — a `genelists` document

```jsonc
{
  "_id":       "<hash> <uid>",    // composite key: "<hash> <uid>"
  "hash":      3127292339,         // murmur3-x86-32 of the validated gene ids (see /validate)
  "uid":       "<firebase uid>",
  "owner":     "Andrew Olson",     // decoded token name / email / uid
  "site":      "sorghum",
  "label":     "my list",
  "n_genes":   214,
  "isPublic":  false,
  "createdAt": "2026-06-23T21:01:37.248Z",
  "deletedAt": "2026-07-06T17:00:00.000Z"  // present ONLY when soft-deleted; absent = active
}
```

The member **gene IDs are not stored** on the document. The `hash` is the link to the genes Solr core: a
list's genes are the docs matching `saved_search:<hash>` in `sorghum_genes11`. Fetch them via `/genes` or
`/search` with `q=saved_search:<hash>`.

> **Note on `saved_search` propagation:** the genes core is rebuilt per release, and cross-release/cross-site
> propagation of `saved_search` hashes is handled by a **separate (TBD) process** — the gene-lists delete /
> force / purge paths deliberately never write Solr. Only `POST /gene_lists/validate` writes `saved_search`
> (at save time).

---

## `POST /gene_lists/validate`

Validate a set of gene IDs against the core, compute the list `hash`, and tag the matching genes.

- **Auth:** none. **Body:** a JSON array of gene-id strings.

```
POST /gene_lists/validate
Content-Type: application/json

["SORBI_3001G000200", "SORBI_3001G000100"]
```

- **Response:** `{ "ids": [...validIds], "missing": [...notFoundIds], "hash": <murmur32> }`
- **Side effect:** atomically adds `hash` to `saved_search` of every valid gene (commit=true). Use the
  returned `hash` (and `ids.length` as `n_genes`) when calling `POST /gene_lists`.
- If this exact id set was already validated, it returns immediately with `missing: []`.

---

## `GET /gene_lists`

List saved lists.

| param | in | type | required | notes |
| --- | --- | --- | --- | --- |
| `site` | query | string | yes | e.g. `sorghum` |
| `isPublic` | query | boolean | yes | `true` = other users' **public** lists for the site; `false` = the **caller's own** lists (needs token) |
| `rows` | query | integer | no | default 20 |
| `includeDeleted` | query | string | no | `active` (default) or `trash` |

- `active` (default): excludes soft-deleted lists.
- `trash`: returns **the caller's own soft-deleted lists** (`deletedAt` set). Requires a valid token
  (`401` if anonymous).
- **Response:** a JSON array of `genelists` documents (soft-deleted ones include `deletedAt`).

> The handler filters by `site` + (`isPublic` or `uid`) + `deletedAt`. (The previous generic handler only
> filtered by `uid`, so results are now more correctly scoped.)

---

## `POST /gene_lists` — save a list

| param | in | type | notes |
| --- | --- | --- | --- |
| `hash` | query | integer | from `/validate` |
| `n_genes` | query | integer | list size |
| `label` | query | string | display name |
| `site` | query | string | e.g. `sorghum` |
| `isPublic` | query | boolean | shareable or private |

- **Auth:** required. Upserts `_id = "<hash> <uid>"`; sets `uid`, `owner`, and `createdAt` (on insert).
- **Response:** `{ "message": "list saved" }`.

---

## `PATCH /gene_lists` — update a list

- **Auth:** required. Query `listId` (the `_id`, required) + JSON body `{ "label"?: string, "isPublic"?: boolean }`.
- **Response:** `{ "message": "list updated", "updated": { ... } }`; `400` if no valid fields; `404` if not owned/found.

---

## `DELETE /gene_lists` — delete a list

| param | in | type | required | notes |
| --- | --- | --- | --- | --- |
| `listId` | query | string | yes | the `_id` |
| `force` | query | boolean | no | `true` = permanent now; omit/false = 30-day soft delete |

- **Auth:** required.
- **Default (soft):** sets `deletedAt`; the list is hidden from the active view but **restorable for 30 days**,
  then purged by a daily cron. → `{ "message": "list marked for deletion (restorable for 30 days)" }`
- **`force=true`:** permanent `deleteOne`. → `{ "message": "list permanently deleted" }`
- `404` if not found / not owned (soft delete also `404`s if already soft-deleted).
- **Never writes Solr `saved_search`.**

---

## `POST /gene_lists/restore` — restore a soft-deleted list

- **Auth:** required. `Content-Type: application/json`. Query `listId` (required).
- Clears `deletedAt` (only works while the doc still exists, i.e. within the 30-day window before purge).
- **Response:** `{ "message": "list restored" }`; `404` if nothing to restore (already purged / not owned).

---

## Lifecycle summary

1. **Save:** `POST /validate` → get `hash` + valid `ids` (genes get tagged `saved_search:<hash>`) → `POST /gene_lists` with `hash`, `n_genes`, `label`, `site`, `isPublic`.
2. **Open a list:** query the genes core `saved_search:<hash>` (via `/genes` or `/search`).
3. **Delete:** `DELETE /gene_lists?listId=…` → soft (30-day trash) by default; `&force=true` for permanent.
4. **Trash + restore:** `GET /gene_lists?…&includeDeleted=trash` to list deleted; `POST /gene_lists/restore?listId=…` to recover.
5. **Auto-purge:** lists soft-deleted > 30 days are hard-deleted daily by
   `scripts/cleanup_expired_genelists.js` (installed as an olson cron job; mongo-only).

## Frontend TODO (what changed)
- Delete is now **soft** (30-day trash); add a **"delete permanently"** option (`force=true`).
- Add a **trash view** (`includeDeleted=trash`) rendering `deletedAt`, plus a **restore** action
  (`POST /gene_lists/restore`).
- No change needed for `validate` / `save` / `update`; the `userData1` store move is transparent to clients.
