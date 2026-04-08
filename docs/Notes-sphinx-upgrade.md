# MOM6 Sphinx toolchain upgrade plan

## Background

The MOM6 documentation build (`docs/`) currently depends on four forked
packages, all maintained by one person (`jr3cermak`) and all last touched
in 2020:

| Package in `requirements.txt` | Upstream | Fork purpose |
|---|---|---|
| `jr3cermak/sphinx @ v3.2.1mom6.4` | `sphinx-doc/sphinx` (now 8.x) | Patch `wrap_displaymath()` to stop double-wrapping `eqnarray`/`align` blocks; LaTeX author `\and` fix; HTML equation post-processing hook |
| `jr3cermak/sphinxcontrib-autodoc_doxygen @ 0.7.13` | `rmcgibbo/sphinxcontrib-autodoc_doxygen` (dormant since 2021) | Repurpose a C++ Doxygen-to-Sphinx bridge into a Fortran one for MOM6 |
| `jr3cermak/sphinx-fortran @ 1.2.2` | `VACUMM/sphinx-fortran` (still maintained, last commit 2025-10) | Add Sphinx 3 compatibility (2020) |
| `jr3cermak/flint @ 0.0.1` | `marshallward/flint` (active through 2022) | Add `setup.py` for pip-installable packaging plus a small Project search API |

This chain forces the whole build to sit on Sphinx 3.2.1, which in turn
drags along a pile of ceiling-pinned transitive dependencies
(`jinja2<3.1`, `sphinxcontrib_applehelp<1.0.8`, `alabaster<0.7.14`,
`setuptools<82.0.0`, etc.). That arrangement is increasingly fragile, and
nothing in the chain is actively maintained except the upstream of
`sphinx-fortran`.

This document is the plan for getting MOM6 off those four forks and onto
a supportable toolchain.

---

## Goals and non-goals

**Goals**

- Build the MOM6 documentation against a current Sphinx (8.x) on a modern
  Python with a current `sphinx-rtd-theme`, without patched Sphinx.
- Eliminate all four `jr3cermak/*` dependencies.
- Leave the rendered documentation visually and structurally equivalent
  to what we publish today, including math, cross-references, tables,
  figures, `:callto:`/`:calledfrom:` links, etc.
- Make future maintenance tractable for MOM6 developers without requiring
  knowledge of a separate GitHub account's release process.

**Non-goals**

- Switching the Doxygen-to-Sphinx bridge to Breathe. Breathe has known
  gaps and conflicts around Fortran and is not an acceptable substitute.
- Replacing Doxygen as the Fortran source parser.
- Restructuring or rewriting the actual documentation content.

---

## Survey of the four forks

Before settling on a plan, we traced each fork to its upstream, read the
commit history on the `jr3cermak` branches, and diffed the code against
upstream to understand the scope of changes. Summary of findings:

### 1. `jr3cermak/sphinx` (forked Sphinx itself)

Seven commits on top of Sphinx 3.2.1 with three functional changes and
some debug leftovers:

- **`sphinx/util/math.py`** - skip Sphinx's default `\begin{split}` /
  `\begin{equation}` wrapping when the LaTeX already contains
  `\begin{equation}`, `\begin{eqnarray}`, or `\begin{align}`. Without
  this, MOM6's math-heavy sources produce broken LaTeX.
- **`sphinx/builders/latex/__init__.py`** - convert comma- and
  "and"-separated author lists into `\and` for the LaTeX builder.
- **`sphinx/cmd/build.py`** - honor an `UPDATEHTMLEQS` environment
  variable to run a local `postProcessEquations.py` after the HTML
  build.
- Commented-out `import pdb; pdb.set_trace()` lines in
  `sphinx/ext/autodoc/__init__.py` and
  `sphinx/transforms/post_transforms/__init__.py`.

Upstream Sphinx 8.x has **not** fixed the equation-wrapping issue; the
`wrap_displaymath()` function is structurally unchanged. The official
workaround is `:nowrap:` on every affected `.. math::` directive, which
would require touching many MOM6 source files.

### 2. `jr3cermak/sphinxcontrib-autodoc_doxygen`

Twenty-three commits ahead of a dormant upstream. The upstream is a
C++ documentation tool (written for OpenMM) that wires Doxygen XML into
Sphinx via the `cpp` domain. The fork effectively rewrites it into a
Fortran documentation tool for MOM6:

- `xmlutils.py` grows from 129 lines to 1059 lines. Almost every visitor
  method is new or rewritten: math blocks with label management,
  `\eqref`/`\eqref2`/`\eqref4` handling, sections (`sect1`..`sect4`),
  nested ordered and itemized lists, tables rendered as RST grid tables,
  figures with captions, footnotes, citations, emphasis, superscript,
  anchors, hyperlinks, `\htmlonly` and `\latexonly` blocks, and build-mode
  aware output (`html` vs `latex`/`latexpdf`).
- `autodoc.py` comments out `DoxygenClassDocumenter` and replaces it with
  two new documenters for Fortran concepts: `DoxygenModuleDocumenter` and
  `DoxygenTypeDocumenter`. `DoxygenMethodDocumenter` is rewritten to
  distinguish `subroutine` vs `function`, handle `result()` clauses, and
  emit `:callto:`/`:calledfrom:` field lists from `<references>` and
  `<referencedby>` XML.
- The entire extension is domain-switched from `cpp` to `f`. All role
  strings (`:cpp:any:` -> `:f:func:`/`:f:mod:`/`:f:type:`) and directive
  emissions (`.. cpp:class::` -> `.. f:module::`, etc.) are changed.
- `autosummary/` is adapted to the new documenters and gains
  auto-discovery of modules and Doxygen pages from the XML tree, with
  new `doxynamespace.rst` and `doxypage.rst` templates.
- The fork carries an optional `flint` dependency to patch up places
  where Doxygen's Fortran parsing is incomplete.

The upstream shares roughly the file layout and a handful of helper
names. Essentially none of the actual rendering logic is shared.

### 3. `jr3cermak/sphinx-fortran`

Fourteen commits ahead of upstream tag 1.1.1, six of them authored by
Cermak. The substantive change is Sphinx 3 API compatibility (logger
calls, directive registration, etc.), plus a `sphinx<4` pin.

Since 2020, the upstream `VACUMM/sphinx-fortran` has continued to
evolve - Sphinx logging API migration, parallel read support, double
precision type handling, removal of the `future` dependency, modern
`setup.cfg`/`setup.py` cleanup through October 2025 - but has not cut a
PyPI release beyond 1.1.1. The fork's Sphinx 3 fixes are subsumed by
what upstream master does today.

### 4. `jr3cermak/flint`

Two commits ahead of upstream. Adds a `setup.py`, a `find_member` /
`search_units` API on `Project`, multi-path parsing, and a verbose
flag. Upstream (`marshallward/flint`, by a MOM6 colleague) has
continued through 2022 and is 22 commits ahead of the fork point on
other axes.

---

## Plan

The plan treats each fork as a separate problem with its own minimal
resolution. The pieces are independent and can land in separate PRs.

### Piece 1 - Vendor `sphinxcontrib-autodoc_doxygen` as a local extension

**Decision: vendor, do not monkey-patch, do not keep as a dependency.**

We considered three options:

1. **Keep depending on the `jr3cermak` fork.** Rejected: unmaintained,
   pinned to Sphinx < 4, hosted on a personal account we do not control,
   and blocks upgrading the rest of the toolchain.
2. **Write a thin wrapper extension that imports the upstream package
   and monkey-patches the pieces we need.** Rejected. Monkey-patching is
   the right tool when a small, stable dependency needs a small nudge.
   Here the fork rewrote ~90% of `xmlutils.py`, replaced both
   documenters, and changed every hardcoded `cpp` domain string to `f`.
   A monkey-patch layer would end up containing all of that code anyway,
   just behind an indirection. And there is no upstream evolution to
   track: upstream has been dormant since mid-2021.
3. **Vendor the extension into the MOM6 repo as a local Sphinx
   extension.** Chosen. The code is small (~1500 lines total), we
   already effectively own it, and putting it in-tree makes it
   debuggable, `git blame`-able, and editable by any MOM6 developer
   without a separate release dance.

**Target layout**

```
docs/
  _ext/
    autodoc_doxygen/
      __init__.py
      autodoc.py
      xmlutils.py
      autosummary/
        __init__.py
        generate.py
        templates/
          doxymodule.rst   # renamed from doxynamespace.rst
          doxypage.rst
```

Note the package is renamed from `sphinxcontrib.autodoc_doxygen` to
plain `autodoc_doxygen`. The `sphinxcontrib` namespace is for
PyPI-distributed extensions; a vendored in-tree extension should not
claim that namespace.

**Changes to `conf.py`**

```python
sys.path.insert(0, os.path.abspath('_ext'))
extensions = [
    'sphinxcontrib.bibtex',
    'sphinx.ext.ifconfig',
    'autodoc_doxygen',          # was 'sphinxcontrib.autodoc_doxygen'
    'sphinxfortran.fortran_domain',
]
```

**Work items**

1. Copy `sphinxcontrib/autodoc_doxygen/` from the fork (tag `0.7.13`)
   into `docs/_ext/autodoc_doxygen/`.
2. Drop dead code that accumulated during the original porting:
   - the commented-out `DoxygenClassDocumenter` registration in
     `__init__.py`
   - `import pdb; pdb.set_trace()` lines throughout
   - the `visit_ref_angus` alternative implementation
   - the `_import_by_name_original` dead function
   - unused `six` imports and Python 2 idioms
3. Port to modern Sphinx (8.x). Known API surface to check:
   - `sphinx.ext.autodoc.Documenter` constructor signature and
     `document_members` signature
   - `sphinx.ext.autosummary.Autosummary.get_items` / `get_table`
     return types
   - `env.ref_context['cpp:parent_key']` - whether this still exists in
     modern Sphinx; the fork already has a fallback path using
     `cpp:parent_symbol`
   - `self.bridge.result` internals used by `DoxygenAutosummary`
   - Any `app.add_autodocumenter` signature changes
4. Decide what to do about the optional `flint` import. Per the source
   audit it is imported but not actually called. Either remove the
   import or make the integration real. See Piece 4.
5. Remove `sphinxcontrib-autodoc_doxygen` from `requirements.txt`.

**Risk**

The Sphinx 3 -> 8 port of the documenters is the only part with real
uncertainty. Autodoc's extension points are mostly stable, but
`Documenter` has had signature changes and the autosummary internals
have moved. Mitigation: start with a fresh venv on current Sphinx, get
the extension to import cleanly, then fix documenters one directive at
a time by running small `.rst` fixtures.

### Piece 2 - Replace `jr3cermak/sphinx-fortran` with upstream

`VACUMM/sphinx-fortran` master supports modern Sphinx and has evolved
past the point where the fork's Sphinx 3 compatibility patch is needed.
The upstream has not cut a new PyPI release, so we pin to a commit.

**Change to `requirements.txt`**

```
-git+https://github.com/jr3cermak/sphinx-fortran.git@1.2.2#egg=sphinx-fortran
+git+https://github.com/VACUMM/sphinx-fortran.git@<pinned-sha>#egg=sphinx-fortran
```

**Risk**

Low. This is a straight dependency swap, and `sphinx-fortran` is
consumed only through the `sphinxfortran.fortran_domain` extension
listed in `conf.py`. Verify during integration that the Fortran domain
roles (`f:mod`, `f:func`, `f:type`) still resolve the cross-references
emitted by the vendored `autodoc_doxygen`.

### Piece 3 - Handle the Sphinx math wrapping fix as a local patch

This is the one place where a monkey-patch is the right tool. The
change is small, the target function (`sphinx.util.math.wrap_displaymath`)
has been stable for years, and the alternative - `:nowrap:` on every
affected directive across the MOM6 source tree - is much worse.

**Approach**

Add a tiny local extension (or a few lines in `conf.py`'s `setup()`)
that replaces `sphinx.util.math.wrap_displaymath` with a version that
detects pre-existing `equation`/`eqnarray`/`align` environments and
returns the content unwrapped.

For the LaTeX author `\and` fix, the simplest thing is to continue
doing what `conf.py` already does: build `latex_authors` ourselves and
pass it to `latex_documents`. That code is already in `conf.py` and
does not require patching Sphinx.

For the `UPDATEHTMLEQS` post-processing hook, move that logic into the
`Makefile` as a post-step after `sphinx-build`, rather than patching
`sphinx.cmd.build`. The hook is only used when the environment variable
is set, so this is purely a build-system concern.

**Upstream contribution (optional)**

It is worth submitting the `wrap_displaymath` fix upstream once we have
it working locally. Sphinx issue #3785 has been open since 2017 and
covers adjacent territory. This is not a blocker for the upgrade.

### Piece 4 - Decide the fate of `flint`

`flint` is imported by the fork's `autodoc_doxygen` but per the audit
is not actually called in the rendering path. Two options:

- **Remove the import.** Simplest. Drop `flint` from `requirements.txt`
  entirely. Verify nothing in the extension or in any `.rst` source
  relies on `flint`-supplied output.
- **Actually use `flint`.** If Doxygen's Fortran parsing genuinely
  leaves gaps around `result()` clauses or argument documentation that
  MOM6 needs, wire `flint` in properly during the `autodoc_doxygen`
  port. In that case, switch from `jr3cermak/flint` to upstream
  `marshallward/flint` (Marshall is a MOM6 colleague and can accept the
  `setup.py` / search-API additions upstream).

Decision deferred until Piece 1 is underway and we can empirically see
whether any MOM6 documentation output degrades without `flint`.

---

## Target `requirements.txt` after the upgrade

Approximately:

```
sphinx
sphinx-rtd-theme
sphinxcontrib-bibtex
git+https://github.com/VACUMM/sphinx-fortran.git@<pinned-sha>#egg=sphinx-fortran
numpy
```

Gone: the four `jr3cermak/*` git dependencies, all of the
`sphinxcontrib_*` version ceilings, `jinja2<3.1`, `alabaster<0.7.14`,
`setuptools<82.0.0`, `six`, `future`. Those pins exist only because
Sphinx 3.2.1 drags them in.

---

## Ordering and validation

1. **Piece 1** first, in its own branch. Do the vendor-and-port in a
   fresh Python 3.11+ virtualenv with current Sphinx. This is the big
   one and blocks the others conceptually (modern Sphinx is only
   useful once `autodoc_doxygen` works on it).
2. **Piece 2** once Piece 1 imports cleanly - swap `sphinx-fortran` and
   verify Fortran-domain cross-references resolve.
3. **Piece 3** in parallel with Piece 1. The math monkey-patch can be
   developed and tested against a trivial standalone document.
4. **Piece 4** last, based on empirical observation during Piece 1.

**Validation at each step**

- `make html` completes without warnings we did not have before. The
  `docs/` build has historically been noisy, so the bar is
  "no *new* warnings", not "zero warnings".
- Spot-check rendered output for: math blocks (both inline and display,
  including `align`/`eqnarray`), `.. f:module::` pages, `:callto:` /
  `:calledfrom:` links, tables from Doxygen `<table>` XML, figures,
  and citations resolved via `sphinxcontrib-bibtex`.
- `make latexpdf` still produces a PDF. This exercises the math
  wrapping fix and the LaTeX author-list fix.
- RTD build: once a branch is ready, trigger a ReadTheDocs build from
  the branch and confirm the published site looks right.

---

## What we are explicitly *not* doing

- **Not switching to Breathe.** Breathe has Fortran gaps and conflicts
  that make it a worse choice than the vendored extension.
- **Not replacing Doxygen.** Doxygen is still the Fortran parser feeding
  the XML that `autodoc_doxygen` consumes.
- **Not restructuring documentation content.** This is a toolchain
  upgrade. RST sources under `docs/` should be touched only where a
  Sphinx 8 deprecation genuinely requires it.
- **Not chasing zero warnings.** The existing build has warnings; we
  will preserve the current signal-to-noise ratio, not try to fix
  longstanding doc issues in the same PR.
