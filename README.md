# foam-extend — fork

Mirror of the official foam-extend repository, plus fixes for defects we hit
running `dynamicTopoFvMesh` in parallel.

Upstream lives on SourceForge, not here:
**<https://sourceforge.net/p/foam-extend/foam-extend-5.0/>**
(`https://git.code.sf.net/p/foam-extend/foam-extend-5.0`)

This fork exists to keep branches somewhere they can be linked and diffed. It is
not where the fixes get reviewed — that happens in the SourceForge ticket
tracker.

## Branches

| branch | what it is |
|---|---|
| `master`, `nextRelease`, `newTmpHandling`, `oversetUpdate` | untouched mirrors of upstream |
| `main` | this branch: the fork's own README and the sync workflow. No foam-extend source. |
| `fix/dynamicTopoFvMesh-null-xfer-master` | the null-`Xfer` fix on top of `master` |
| `fix/dynamicTopoFvMesh-null-xfer` | the same fix on top of `oversetUpdate` |

The mirror branches are kept byte-identical to upstream, which is why the fork's
own files live on a separate `main` branch. Putting them on `master` would make
it diverge and turn every sync into a merge.

## Staying in sync

`.github/workflows/sync-upstream.yml` runs daily and can be triggered by hand
from the Actions tab. It clones upstream, then pushes every branch that exists
*there* to the matching branch here, along with tags.

Two properties worth knowing:

- **It cannot delete or overwrite fork-only branches.** It pushes only refs that
  exist upstream, and `fix/*` does not.
- **It never force-pushes.** If upstream rewrites history, the push is rejected
  and the job fails instead of quietly rewriting the mirror. A failed run means
  upstream diverged and wants a look, not that the job is broken.

## The fixes

Both defects are in `dynamicTopoFvMesh` and both sit on code paths that only
parallel runs reach, which is why they went unnoticed for years — a serial run
returns before it gets there.

**1. Null `Xfer` references passed to `polyMesh::resetPrimitives`** —
SIGSEGV on the first mesh update of any parallel run, at any rank count.
`moveCoupledSubMeshes()` builds three `Xfer` arguments by dereferencing the null
pointer. `polyMesh::resetPrimitives` used to test them with `if (&fcs)`; commit
`078f96f80` (2015) correctly replaced that undefined-behaviour test with
`if (!fcs().empty())`, which dereferences — but never updated this caller.

Reported as **[ticket #97](https://sourceforge.net/p/foam-extend/tickets/97/)**.
Verified against the shipped `circCylinder3d` tutorial at 4 ranks and on a
56832-cell tetrahedral case at 2, 4, 8 and 16.

**2. Division by zero in `lengthScaleEstimator::readLengthScaleInfo`** —
`sumLength /= nTouchedNgb` with nothing guaranteeing the count is non-zero.
With `fieldRefinement`, the sweep seeds level 1 into interior cells; such a cell
on a processor patch has no neighbour at a strictly lower level and no
fixed-scale patch face, so the count is zero and the run dies with SIGFPE.

Not yet reported — the fix still has to be written and tested.

## Licence

foam-extend is GPL-3.0. Nothing here changes that.
