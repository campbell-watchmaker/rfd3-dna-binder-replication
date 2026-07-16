# Contributing / workflow

`main` is protected: direct pushes and force-pushes are rejected (verified — a
direct push attempt is bounced with `GH006: Protected branch update failed ...
Changes must be made through a pull request`), and every change merges through a
pull request, with `enforce_admins` on so this isn't bypassable by accident.

There is currently one maintainer, and GitHub does not allow self-approval of a
PR — requiring an approving review would make every PR permanently unmergeable.
So `required_approving_review_count` is set to 0: **a PR is still required**,
but merging isn't gated on a review that's structurally impossible to obtain
solo. If a second maintainer joins, raise this back to 1.

1. Branch off `main` for each unit of work (`git checkout -b <topic>`).
2. Commit, push the branch, open a PR describing what changed and why.
3. CI (`.github/workflows/tests.yml`) runs the test suite in `tests/` on every PR —
   currently a geometric/chemical self-consistency check on the DNA target-prep
   math (`tests/test_target_prep.py`).
4. Read the diff and merge via GitHub once CI is green and the change looks
   right — not before, and not automatically.

Method or data-processing changes (a new relaxation strategy, a different filter
threshold, a substituted tool) get a PR description that states what changed and
why, the same way `docs/replication_log.md` documents substitutions from the
original papers.
