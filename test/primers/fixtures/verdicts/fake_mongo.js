'use strict';

// Fake gramene-mongodb-config for annotate/run tests. The collection evaluates the query subset the
// check uses: equality, $in, $lte, $gte, dotted paths through arrays; find(query, {fields}).limit(n).toArray().
// fakeMongo(docs, {fail, hang, noCollection, script}) → { mongo, collection, calls: [{query, options, limit, n}] }
// script(call) (per find, n = its 0-based index) returns 'ok', 'hang', 'fail', an Error (rejected with it) or a
// promise (returned by toArray as is); it replaces fail/hang.

function valuesAt(doc, path) {
  let cur = [doc];
  for (const part of path.split('.')) {
    const next = [];
    for (const v of cur) {
      if (v == null) continue;
      const x = v[part];
      if (Array.isArray(x)) next.push(...x);
      else if (x !== undefined) next.push(x);
      else if (Array.isArray(v)) for (const e of v) if (e && e[part] !== undefined) next.push(e[part]);
    }
    cur = next;
  }
  return cur;
}

function condMatches(value, cond) {
  if (cond !== null && typeof cond === 'object' && !Array.isArray(cond)) {
    for (const op of Object.keys(cond)) {
      const arg = cond[op];
      if (op === '$in') { if (!arg.some((x) => x === value)) return false; }
      else if (op === '$lte') { if (!(value <= arg)) return false; }
      else if (op === '$gte') { if (!(value >= arg)) return false; }
      else throw new Error('fake mongo: unsupported operator ' + op);
    }
    return true;
  }
  return value === cond;
}

function matches(doc, query) {
  return Object.keys(query).every((path) => valuesAt(doc, path).some((v) => condMatches(v, query[path])));
}

function fakeMongo(docs, opts) {
  const o = opts || {};
  const calls = [];
  const collection = {
    find(query, options) {
      const call = { query, options, limit: null, n: calls.length };
      calls.push(call);
      const cursor = {
        limit(n) {
          call.limit = n;
          return cursor;
        },
        toArray() {
          const how = o.script ? o.script(call) : o.hang ? 'hang' : o.fail ? 'fail' : 'ok';
          if (how && typeof how.then === 'function') return how;
          if (how instanceof Error) return Promise.reject(how);
          if (how === 'hang') return new Promise(() => {});
          if (how === 'fail') return Promise.reject(new Error('fake mongo failure'));
          const out = docs.filter((d) => matches(d, query));
          return Promise.resolve(call.limit ? out.slice(0, call.limit) : out);
        }
      };
      return cursor;
    }
  };
  const mongo = { genes: { mongoCollection: async () => (o.noCollection ? undefined : collection) } };
  return { mongo, collection, calls };
}

module.exports = { fakeMongo, matches, valuesAt };
