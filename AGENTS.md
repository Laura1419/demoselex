<!-- worksgood-managed-guide:v1:start -->
# WorksGood (`wg`) project guide

This file is the **layer-2** project guide for agents working in this
WorksGood task graph. It is NOT the universal chat-agent / worker-agent
contract — that is bundled inside the `wg` binary and emitted by:

```
wg agent-guide
```

Run `wg agent-guide` at session start (or read its output from a previous
session) to get the universal role contract: chat agent vs dispatcher vs worker
distinction, `## Validation` requirement, smoke-gate, cycle handling, git
hygiene, worktree isolation, "no built-in Task tool" rules, etc.

This file only covers things specific to this project. Add project-specific
build commands, test commands, architecture notes, and service recipes here.

**At the start of each session, run `wg quickstart` in your terminal to orient yourself.**
Use `wg service start` to dispatch work — do not manually claim tasks.

This guide is written to both `CLAUDE.md` and `AGENTS.md` and kept in
lock-step. The two files exist because Claude Code and Codex CLI look for
different filenames, but they should never drift in content. Any divergence is
a bug. Update both together.
<!-- worksgood-managed-guide:v1:end -->
## Project context: trajectory-replication research

Read `docs/MISSION.md` first; it is authoritative. All paths below are relative to your
current working directory (the repository root). Use relative paths, not absolute ones.

- Living research doc sections live in `docs/` (`A-` .. `J-`). Update the section your task owns.
- Append a short entry to `docs/JOURNAL.md` after each work session.
- Record provenance for every run with `python3 scripts/provenance.py` (see `docs/PROVENANCE.md`).
- Lean code in `lean/`; Python analysis in `sim/` and `analysis/`; figures in `fig/`.
- Raw data in `data/raw/` is immutable. Never edit it. Derived data in `data/derived/`.
- Commit your work with a conventional-commit message naming the task id, e.g.
  `feat(traj-prototype): time-aware clustering prototype`.
- Do not mark a task done until its `## Validation` section actually passes.
- Honest negative results are valuable: if a hypothesis fails, say so in the doc and journal.

## Repository conventions (mandatory)

- Default branch is **`main`**. Never assume `master`.
- The completion contract enforces `git diff --check refs/heads/main..HEAD`. A
  `pre-commit` hook at `.githooks/pre-commit` strips trailing whitespace, so keep it
  enabled (`git config core.hooksPath .githooks` is already set once per repo and shared
  by worktrees).
- Do not add trailing whitespace; do not add conflict markers.
- Every commit message: `<type>(<task-id>): <summary>`.

- **Never run `git merge`, `git rebase`, `git pull`, or `git stash pop`.** Just edit,
  `git add`, `git commit` on the current branch. The coordinator handles integration.
- Keep edits scoped to the files your task owns. Never rewrite shared files owned by
  other tasks (`docs/JOURNAL.md`, `docs/PROVENANCE.md`, `docs/DECISIONS.md`): append only.
