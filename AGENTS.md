# AGENTS.md - OpenFOAM

This is OpenFOAM v2606 (API 2606), a large C++ CFD codebase. It uses its own build system (`wmake`), not CMake/Make. Almost nothing here is standard tooling — read below before guessing.

## Environment Setup

Every build or run command requires the OpenFOAM environment first:

**Compilation set up** (from repo root):
```bash
#source $WM_PROJECT_DIR/etc/bashrc
source $PWD/eb-r11-foss2025b.env
source $PWD/prefs.opt.sh
```

Without this, `wmake`, `wclean`, `lnInclude`, and all `FOAM_*`/`WM_*` variables are undefined.

## Build Commands

**Full build** (from repo root):
```bash
./Allwmake
```

**Build a single library** (e.g. finiteVolume):
```bash
wmake $targetType finiteVolume    # from src/
```

**Build a single application** (e.g. simpleFoam):
```bash
wmake $targetType simpleFoam      # from the app's Make/ directory
```

**Build with parallel jobs:**
```bash
wmake -j 8 $targetType <target>
```

**Recurse into subdirectories:**
```bash
wmake -all $targetType <dir>      # runs Allwmake in each subdir
```

**Clean build artifacts:**
```bash
wclean [target] [dir]
wcleanLnIncludeAll
```

Build order matters — libraries have hard dependencies. See `src/Allwmake` for the canonical order: OpenFOAM → Pstream → fileFormats → surfMesh → meshTools → finiteVolume → ...

## Key Directories

| Path | What it is |
|------|-----------|
| `src/` | Core libraries (OpenFOAM, finiteVolume, TurbulenceModels, etc.) |
| `src/OpenFOAM/` | Foundation library — everything depends on this |
| `applications/solvers/` | CFD solvers (simpleFoam, reactingFoam, etc.) |
| `applications/utilities/` | Pre/post-processing tools (e.g. foamFormatConvert) |
| `applications/tools/` | Development/admin tools |
| `wmake/` | Custom build system — do not replace with cmake/make |
| `etc/` | Shell config, runtime dictionaries, environment setup |
| `tutorials/` | Example cases with `Allrun`/`Allclean`/`Alltest` scripts |
| `modules/` | Optional submodules (adios, OpenQBMM, visualization) |
| `plugins/` | Community plugins (cfmesh, bindings, turbulence-community) |
| `platforms/$WM_OPTIONS/` | Build output: `bin/` and `lib/` |
| `test/` | User test/analysis scripts (not the official test suite) |
| `META-INFO/` | Internal api/build info — edit `api-info` carefully |

## wmake Build System

Each library/app has a `Make/` directory with:
- **`files`** — source files and output target (`LIBSO` for .so, `EXE` for binaries)
- **`options`** — include paths (`EXE_INC`/`LIB_INC`) and link libraries (`EXE_LIBS`/`LIB_LIBS`)

Dependencies are expressed via `lnInclude` symlinks — run `wmakeLnInclude -u <lib>` after changing header locations.

## Commit Message Conventions

From CONTRIBUTING.md — use these prefix tags:
- `BUG:` — bug fixes
- `ENH:` — new/enhanced functionality
- `DOC:` — documentation
- `COMP:` — compiler/build changes
- `CONFIG:` — configuration (e.g. config.sh)
- `STYLE:` — code style cleanup
- `TUT:` — tutorial changes
- `SUBMODULE:` — module/plugin changes
- `DEFEATURE:` — deprecation/removal
- `REVERT:` — revert commits

Format: `TAG: context: imperative message` (≤72 chars title, blank second line, imperative mood).

## Testing

- `tutorials/Allrun` — runs tutorial cases and reports results
- `tutorials/Alltest` — quick tutorial test with scheme/solver info
- `tutorials/Allclean` — clean tutorial run directories
- Individual tutorials can be tested by running their local `Allrun` script

## Common Gotchas

- **Never edit `etc/bashrc` directly** — use `etc/prefs.sh` or `~/.OpenFOAM/$FOAM_API/prefs.sh` for overrides.
- **Build output goes to `platforms/$WM_OPTIONS/`** (e.g. `platforms/linux64GccDPInt32Opt/`).
- **`wmake` expects `$WM_PROJECT_DIR` to be set** — it checks this at startup.
- **Library compile order is critical** — OpenFOAM core must be built before finiteVolume, which must be built before solvers. Follow `src/Allwmake` order.
- **`modules/` and `plugins/` are git submodules** — initialize with `git submodule update --init`.
- **Modules require `$FOAM_MODULE_PREFIX`** to be set (not `false`/`none`); otherwise they are skipped.
- **The `test/` directory contains user analysis scripts**, not the official test suite — the official tests are in `tutorials/`.
- **`$targetType`** is a wmake variable used in Allwmake scripts — don't hardcode it.

## GitLab Workflow

- MR target branch: `develop`
- Use `community-contributions` MR template for contributions
- Do not force-push after review; use `SQUASH` prefix commits instead
- Fork-based workflow: fork → branch → MR to upstream `develop`
