'use strict';

// Lua scripts for the Redis job store (spec §A.8.1). Registered with ioredis defineCommand.
// Scripts derive job/result keys from a prefix ARGV, which is fine on standalone Redis 5.
//
// Job docs are JSON written by JavaScript. The scripts only DECODE them (status / worker checks);
// they never re-encode, because cjson turns empty arrays into {}.

// Returns status, worker of a raw job JSON (nil, nil for missing/corrupt).
const JOB_STATUS_FN = `
local function jobStatus(raw)
  if not raw then return nil, nil end
  local ok, doc = pcall(cjson.decode, raw)
  if not ok or type(doc) ~= 'table' then return nil, nil end
  local st = doc.status
  if type(st) ~= 'string' then st = nil end
  local w = doc.worker
  if type(w) ~= 'string' then w = nil end
  return st, w
end
local function inList(csv, value)
  if type(value) ~= 'string' then return false end
  return string.find(',' .. csv .. ',', ',' .. value .. ',', 1, true) ~= nil
end
`;

// KEYS: job, targetQueue, otherQueue, finished
// ARGV: jobJson, id, ttlActiveS, maxQueued
// -> {'EXISTS', currentJson} | {'FULL', ''} | {'QUEUED', ''}
// An existing job in any status but 'error' is returned as is; an errored job is replaced and re-queued.
const SUBMIT = JOB_STATUS_FN + `
local cur = redis.call('GET', KEYS[1])
if cur then
  local st = jobStatus(cur)
  if st and st ~= 'error' then return {'EXISTS', cur} end
end
if redis.call('LLEN', KEYS[2]) + redis.call('LLEN', KEYS[3]) >= tonumber(ARGV[4]) then return {'FULL', ''} end
redis.call('SET', KEYS[1], ARGV[1], 'EX', tonumber(ARGV[3]))
redis.call('LREM', KEYS[2], 0, ARGV[2])
redis.call('RPUSH', KEYS[2], ARGV[2])
redis.call('ZREM', KEYS[4], ARGV[2])
return {'QUEUED', ''}
`;

// KEYS: queueSpec, queuePan, running, globalSlots, panSlots
// ARGV: nowMs, staleMs, globalMax, localMax, panMax, siteKey, jobKeyPrefix
// -> {id, 'spec'|'pan'} | nil
// Prunes stale host-wide slots, enforces global/local/pan caps, drops expired or non-queued ids.
const CLAIM = JOB_STATUS_FN + `
local now = tonumber(ARGV[1])
redis.call('ZREMRANGEBYSCORE', KEYS[4], '-inf', now - tonumber(ARGV[2]))
redis.call('ZREMRANGEBYSCORE', KEYS[5], '-inf', now - tonumber(ARGV[2]))
if redis.call('ZCARD', KEYS[4]) >= tonumber(ARGV[3]) or redis.call('ZCARD', KEYS[3]) >= tonumber(ARGV[4]) then return false end
local function pop(q)
  while true do
    local id = redis.call('LPOP', q)
    if not id then return nil end
    local j = redis.call('GET', ARGV[7] .. id)
    if j and jobStatus(j) == 'queued' then return id end
  end
end
local id = pop(KEYS[1])
local pan = false
if not id and redis.call('ZCARD', KEYS[5]) < tonumber(ARGV[5]) then
  id = pop(KEYS[2])
  pan = id ~= nil
end
if not id then return false end
local member = ARGV[6] .. ':' .. id
redis.call('ZADD', KEYS[3], now, id)
redis.call('ZADD', KEYS[4], now, member)
if pan then
  redis.call('ZADD', KEYS[5], now, member)
  return {id, 'pan'}
end
return {id, 'spec'}
`;

// KEYS: worker   ARGV: token, ttlMs
// -> 1 refreshed | 2 re-acquired (key had expired) | 0 held by another worker
const LOCK_REFRESH = `
local v = redis.call('GET', KEYS[1])
if v == ARGV[1] then
  redis.call('PEXPIRE', KEYS[1], tonumber(ARGV[2]))
  return 1
end
if not v then
  redis.call('SET', KEYS[1], ARGV[1], 'PX', tonumber(ARGV[2]))
  return 2
end
return 0
`;

// KEYS: worker   ARGV: token   -> 1 released | 0 not ours
const LOCK_RELEASE = `
if redis.call('GET', KEYS[1]) == ARGV[1] then
  redis.call('DEL', KEYS[1])
  return 1
end
return 0
`;

// Compare-and-set a job doc.
// KEYS: job   ARGV: jobJson, ttlS, expectStatusesCsv, expectWorker ('' = any)   -> 1 | 0
const SET_JOB = JOB_STATUS_FN + `
local st, w = jobStatus(redis.call('GET', KEYS[1]))
if not inList(ARGV[3], st) then return 0 end
if ARGV[4] ~= '' and w ~= ARGV[4] then return 0 end
redis.call('SET', KEYS[1], ARGV[1], 'EX', tonumber(ARGV[2]))
return 1
`;

// KEYS: job, partial   ARGV: gzipBuffer, ttlS, expectWorker   -> 1 | 0 (job not running for this worker)
const SET_PARTIAL = JOB_STATUS_FN + `
local st, w = jobStatus(redis.call('GET', KEYS[1]))
if st ~= 'running' then return 0 end
if ARGV[3] ~= '' and w ~= ARGV[3] then return 0 end
redis.call('SET', KEYS[2], ARGV[1], 'EX', tonumber(ARGV[2]))
return 1
`;

// KEYS: job, running, globalSlots, panSlots
// ARGV: id, member, nowMs, ttlActiveS, pan ('1'|'0'), expectWorker   -> 1 | 0
const HEARTBEAT = JOB_STATUS_FN + `
local st, w = jobStatus(redis.call('GET', KEYS[1]))
if st ~= 'running' then return 0 end
if ARGV[6] ~= '' and w ~= ARGV[6] then return 0 end
redis.call('ZADD', KEYS[2], tonumber(ARGV[3]), ARGV[1])
redis.call('ZADD', KEYS[3], tonumber(ARGV[3]), ARGV[2])
if ARGV[5] == '1' then redis.call('ZADD', KEYS[4], tonumber(ARGV[3]), ARGV[2]) end
redis.call('EXPIRE', KEYS[1], tonumber(ARGV[4]))
return 1
`;

// complete (hasResult '1') or fail (hasResult '0'), then trim `finished`.
// KEYS: job, result, partial, running, globalSlots, panSlots, finished
// ARGV: id, member, jobJson, ttlS, nowMs, resultGzip, hasResult, expectWorker, maxFinished,
//       jobKeyPrefix, resultKeyPrefix, expectStatusesCsv
// -> {1, trimmed} | {0, 0}
const FINISH = JOB_STATUS_FN + `
local st, w = jobStatus(redis.call('GET', KEYS[1]))
if not inList(ARGV[12], st) then return {0, 0} end
if ARGV[8] ~= '' and w ~= ARGV[8] then return {0, 0} end
if ARGV[7] == '1' then
  redis.call('SET', KEYS[2], ARGV[6], 'EX', tonumber(ARGV[4]))
else
  redis.call('DEL', KEYS[2])
end
redis.call('SET', KEYS[1], ARGV[3], 'EX', tonumber(ARGV[4]))
redis.call('DEL', KEYS[3])
redis.call('ZREM', KEYS[4], ARGV[1])
redis.call('ZREM', KEYS[5], ARGV[2])
redis.call('ZREM', KEYS[6], ARGV[2])
redis.call('ZADD', KEYS[7], tonumber(ARGV[5]), ARGV[1])
local trimmed = 0
local excess = redis.call('ZCARD', KEYS[7]) - tonumber(ARGV[9])
if excess > 0 then
  local old = redis.call('ZRANGE', KEYS[7], 0, excess - 1)
  for _, oid in ipairs(old) do
    local jk = ARGV[10] .. oid
    local ost = jobStatus(redis.call('GET', jk))
    if ost == nil or ost == 'done' or ost == 'error' then
      redis.call('DEL', jk, ARGV[11] .. oid)
      trimmed = trimmed + 1
    end
    redis.call('ZREM', KEYS[7], oid)
  end
end
return {1, trimmed}
`;

// Put a job back on its queue (front or back) and release its slots.
// KEYS: job, queue, running, globalSlots, panSlots, partial
// ARGV: id, member, jobJson, ttlS, expectWorker, front ('1'|'0'), expectStatusesCsv   -> 1 | 0
const REQUEUE = JOB_STATUS_FN + `
local st, w = jobStatus(redis.call('GET', KEYS[1]))
if not inList(ARGV[7], st) then return 0 end
if ARGV[5] ~= '' and w ~= ARGV[5] then return 0 end
redis.call('SET', KEYS[1], ARGV[3], 'EX', tonumber(ARGV[4]))
redis.call('LREM', KEYS[2], 0, ARGV[1])
if ARGV[6] == '1' then
  redis.call('LPUSH', KEYS[2], ARGV[1])
else
  redis.call('RPUSH', KEYS[2], ARGV[1])
end
redis.call('ZREM', KEYS[3], ARGV[1])
redis.call('ZREM', KEYS[4], ARGV[2])
redis.call('ZREM', KEYS[5], ARGV[2])
redis.call('DEL', KEYS[6])
return 1
`;

// Queued-key TTL sweep: EXPIRE every queued job to ttl_active; drop ids whose job is gone or not queued.
// KEYS: queueSpec, queuePan   ARGV: jobKeyPrefix, ttlActiveS   -> {refreshed, removed}
const SWEEP = JOB_STATUS_FN + `
local refreshed = 0
local removed = 0
for qi = 1, 2 do
  local ids = redis.call('LRANGE', KEYS[qi], 0, -1)
  for _, id in ipairs(ids) do
    local jk = ARGV[1] .. id
    if jobStatus(redis.call('GET', jk)) == 'queued' then
      redis.call('EXPIRE', jk, tonumber(ARGV[2]))
      refreshed = refreshed + 1
    else
      removed = removed + redis.call('LREM', KEYS[qi], 0, id)
    end
  end
end
return {refreshed, removed}
`;

const COMMANDS = Object.freeze({
  primersSubmit: { numberOfKeys: 4, lua: SUBMIT },
  primersClaim: { numberOfKeys: 5, lua: CLAIM },
  primersLockRefresh: { numberOfKeys: 1, lua: LOCK_REFRESH },
  primersLockRelease: { numberOfKeys: 1, lua: LOCK_RELEASE },
  primersSetJob: { numberOfKeys: 1, lua: SET_JOB },
  primersSetPartial: { numberOfKeys: 2, lua: SET_PARTIAL },
  primersHeartbeat: { numberOfKeys: 4, lua: HEARTBEAT },
  primersFinish: { numberOfKeys: 7, lua: FINISH },
  primersRequeue: { numberOfKeys: 6, lua: REQUEUE },
  primersSweep: { numberOfKeys: 2, lua: SWEEP }
});

// Register every script on an ioredis client (idempotent per client).
function define(client) {
  if (client.__primersLuaDefined) return client;
  Object.keys(COMMANDS).forEach(function (name) {
    client.defineCommand(name, { numberOfKeys: COMMANDS[name].numberOfKeys, lua: COMMANDS[name].lua });
  });
  Object.defineProperty(client, '__primersLuaDefined', { value: true });
  return client;
}

module.exports = {
  COMMANDS,
  define,
  SUBMIT,
  CLAIM,
  LOCK_REFRESH,
  LOCK_RELEASE,
  SET_JOB,
  SET_PARTIAL,
  HEARTBEAT,
  FINISH,
  REQUEUE,
  SWEEP
};
