---
name: julia-best-practices
description: Best practices for writing Julia packages and fast Julia code, following the official Julia style guide and performance tips. Use this skill whenever the user writes, reviews, refactors, or debugs Julia code (.jl files, Project.toml, runtests.jl), creates or restructures a Julia package, searches for or adds a new Julia package or dependency (checking maintenance, version stability and compatibility), asks why Julia code is slow or allocates a lot, or mentions type stability, @code_warntype, BenchmarkTools, profiling (Profile, ProfileView, PProf), Revise, structs and keyword constructors (Base.@kwdef, @with_kw, Parameters.jl), startup latency, sysimages or PackageCompiler (including VS Code settings.json sysimage setup), multiple dispatch or functions with the same name and different methods, exports, or [compat] entries — even if they don't explicitly ask for "best practices".
---

# Julia best practices

Use this when writing or reviewing Julia code in this repo. The goal is code that a
seasoned Julia developer would accept in a registered package: idiomatic, generic,
type-stable, and tested. The rules come from the official Julia manual
(Style Guide and Performance Tips); the reasons are given so you can judge edge cases.

## 1. Package structure

A package should look like this:

```
MyPkg/
├── Project.toml        # name, uuid, version, [deps], [compat]
├── src/
│   ├── MyPkg.jl        # module MyPkg ... end — includes and exports only
│   └── feature.jl      # one file per topic, include()d from MyPkg.jl
├── test/
│   └── runtests.jl     # using Test, MyPkg; @testset blocks
├── docs/               # Documenter.jl (make.jl, src/index.md)
└── README.md
```

- Create new packages with `Pkg.generate("MyPkg")` or PkgTemplates.jl rather than by
  hand, so the UUID and layout are right.
- Give every dependency **and** `julia` a `[compat]` entry. Without them, the General
  registry rejects the package and users get silent breakage on major releases.
- Don't commit `Manifest.toml` for a package (do commit it for an application or a
  reproducible analysis project). A package must work across compatible versions.
- Keep the module file thin: `using`/`import`, `include`s, `export`s. Logic goes in
  topic files.
- Export only the stable user-facing API. On Julia ≥ 1.11, mark other supported names
  with `public`. Everything else is internal and can change freely.
- Test-only dependencies go in `[extras]` + `[targets]` (or `test/Project.toml`), not in
  `[deps]`.
### Development workflow with Revise.jl

Load Revise **before** the package, so edits to `src/` take effect in the running
session without restarting Julia (restarting costs recompilation time):

```julia
using Revise          # or put it in ~/.julia/config/startup.jl
using Pkg; Pkg.activate(".")   # or ]dev path/to/MyPkg from another env
using MyPkg
```

- Revise tracks files loaded via `using`/`include` inside the package. For a standalone
  script, use `includet("script.jl")` instead of `include`.
- Revise **cannot** apply redefinitions of `struct`s (until Julia 1.12, which supports
  redefining structs natively). When a struct changes, restart the session. While a
  struct's design is still changing, a `NamedTuple` or a temporarily renamed struct
  avoids frequent restarts.
- If Revise reports it couldn't apply a change, the session may hold stale code:
  restart instead of debugging ghosts.

## 2. Choosing and adding new packages

Every dependency is code the project now has to live with: an abandoned package or a
version conflict can block upgrades for years. So before suggesting or adding a
package, check it — don't pick one just because it's the first search hit or because
you remember its name. Tell the user what you checked.

### Is it maintained?

Look at the package's repository (usually GitHub, linked from JuliaHub or the General
registry) and check:

- **Recent activity:** commits and a release within roughly the last year. Old isn't
  automatically bad for a small, finished package, but a stale package with many open
  bug reports is a warning sign.
- **Recent Julia support:** its `[compat]` allows the Julia version the project uses,
  and CI runs on current Julia releases and passes.
- **Issues and PRs:** are maintainers responding? Are there unanswered reports like
  "doesn't work on Julia 1.x" or "please tag a new release"?
- **Bus factor and home:** packages under an organization (JuliaLang, JuliaData,
  SciML, JuliaStats, JuliaPlots, …) or with several maintainers are safer than a
  one-person project with no recent activity.
- **Deprecation notes:** the README may say "this package is deprecated, use X".
  Check for that before recommending it, because several older Julia packages have
  well-known successors.
- **Usage:** how many registered packages depend on it (JuliaHub shows dependents)
  and whether it has documentation and tests.
- **Weight:** how many dependencies it pulls in and how much it adds to load time.
  Prefer the standard library or a lightweight package when a heavy one is only needed
  for one small function.

If the evidence is mixed, present the options with trade-offs instead of choosing
silently. If you can't check (no network access), say so and list what the user
should check themselves.

### Is the version stable?

- Versions `0.x` can break in any minor release (`0.4` → `0.5`); `≥ 1.0` promises
  semver stability within a major version. Prefer `≥ 1.0` for core functionality, and
  pin compat tightly for `0.x` packages.
- Check the changelog or release notes for frequent breaking changes and for whether
  the latest release is recent and not a pre-release.
- Write a `[compat]` entry right away. Pkg's default `"1.2"` style means
  `[1.2.0, 2.0.0)`, which is usually what you want.

### Does it work with what's already installed?

Adding a package can force Pkg to downgrade other packages or fail to resolve.
Check before and after:

```julia
using Pkg
Pkg.activate(".")             # the project's own environment, never the global one
Pkg.status()                  # what's installed now, and which versions
Pkg.add("NewPkg"; preserve=Pkg.PRESERVE_ALL)   # try without touching existing versions
Pkg.status(; outdated=true)   # what is held back, and why (⌅ / ⌃ markers)
```

- If `PRESERVE_ALL` fails, read the resolver error before using a plain `Pkg.add`:
  it names the conflicting package and versions. Tell the user what would be
  downgraded or upgraded and let them decide.
- After adding, compare `Pkg.status()` to before. Watch for any package that was
  downgraded or is now held back (shown as `⌅`), and run the test suite.
- Use `Pkg.why("SomePkg")` (Julia ≥ 1.9) to see why an indirect dependency is present.
- Never add packages to the global `@v1.x` environment for project work; keep each
  project in its own environment (only small dev tools such as Revise, BenchmarkTools
  or PProf belong in the global one).
- For a package (not an application), check that its `[compat]` bounds don't
  conflict with the bounds of packages it's commonly used with. Too-tight bounds make
  your package impossible to install alongside others.

## 3. Style (official Julia style guide)

- **Names:** modules and types `CamelCase`; functions and variables `lowercase`
  (underscores only when it reads badly otherwise); constants may be `UPPER_CASE`.
- **Mutating functions end in `!`** and take the mutated argument first:
  `sort!(v)`, `mul!(C, A, B)`. Offer a non-mutating version when it's natural.
- **Write functions, not scripts.** Code at top level is slow and hard to test.
- **Don't over-type arguments.** Type annotations on arguments are for dispatch, not
  speed. Prefer `f(x::AbstractVector{<:Real})` or no annotation over `f(x::Vector{Float64})`,
  so ranges, views, `Float32`, and GPU arrays still work.
- **Use dispatch instead of `isa` chains** (see "Multiple dispatch" below).
- **Avoid type piracy:** don't add methods to functions you don't own for types you
  don't own (e.g. `Base.:+(::Vector, ::String)`). It changes behaviour for every
  other package loaded in the session.
- **Write docstrings** for exported names, starting with the signature in an indented
  block or code fence, then one sentence on what it does, then `# Arguments` /
  `# Examples` as needed. Examples written as `jldoctest` get tested by Documenter.
- Prefer `nothing`/`missing` over sentinel values, and `throw(ArgumentError(...))`
  (or a domain-specific exception type) over `error("...")` for bad input.
- Don't parenthesize conditions (`if x > 0`, not `if (x > 0)`), and use 4-space indent.

### Multiple dispatch: one function name, many methods

In Julia it is normal — and intended — to have many functions with the **same name**
that differ only in their argument types. A *function* is the name (`area`); each
definition for particular argument types is a *method* of it. When you call it, Julia
picks the most specific method matching the types of **all** arguments (hence
"multiple" dispatch, not just the first argument like in OOP languages).

```julia
abstract type Shape end
struct Circle <: Shape; r::Float64; end
struct Rect   <: Shape; w::Float64; h::Float64; end

area(c::Circle) = π * c.r^2
area(r::Rect)   = r.w * r.h

# dispatch on two arguments at once
describe(a::Circle, b::Circle) = "two circles"
describe(a::Circle, b::Shape)  = "a circle and another shape"
describe(a::Shape,  b::Shape)  = "two shapes"   # generic fallback for other combinations

total_area(shapes) = sum(area, shapes)   # works for every Shape, incl. future ones
```

Rules of thumb:

- **One concept, one name.** Write `area(::Circle)` and `area(::Rect)`, not
  `area_circle` / `area_rect` or `area(shape, kind::String)`. The name says *what*,
  the argument types say *how*.
- **Generic fallback + specializations.** Write a general method on an abstract type
  (or no annotation) that works everywhere, then add specialized methods only where a
  type allows something better (e.g. a faster algorithm for `Diagonal` matrices).
- **Extend existing functions for your types** instead of inventing new names:
  `Base.show`, `Base.length`, `Base.iterate`, `Base.:+`, `Base.isless`, … Your type then
  works with printing, `for` loops, `sum`, `sort` and other packages automatically.
  This is fine because you own the type (it's not piracy). Use `import Base: show`
  or write `Base.show(...)` explicitly.
- **All methods share one meaning.** Methods of the same function should do the same
  conceptual thing; if a method would do something different, give it a different
  name. Document the function once, on the generic method.
- **Avoid ambiguities:** if `f(::Int, ::Any)` and `f(::Any, ::Int)` both exist,
  `f(1, 2)` errors. Add the more specific `f(::Int, ::Int)`, or restructure.
  `Aqua.test_ambiguities` detects these.
- **Dispatch on types, not values.** Don't create types just to turn every runtime
  value into a method choice; a plain `if` is fine for value-based logic. `Val{x}` is
  for the rare case where the value is known at compile time.
- **Inspecting:** `methods(area)` lists all methods; `@which area(c)` shows which one a
  call uses; `@edit area(c)` opens its source.

## 4. Structs and keyword constructors

Structs plus multiple dispatch are the foundation of Julia design: types describe
*what the data is*, and methods dispatched on those types describe *what to do with
it*. So actively reach for a struct whenever it's reasonable, rather than defaulting
to loose variables, `Dict`s, or long argument lists.

### When to define a struct

Define one when:

- **Several values travel together** — the same 3+ arguments are passed through many
  functions (e.g. `n, L, dx, bc`), or a function returns several related values.
  Group them as `Grid`, `SolverOptions`, `Result`, …
- **Data has a meaning or invariants** — an `Interval` with `lo ≤ hi`, a `Probability`
  in `[0, 1]`, units of a physical quantity. The constructor enforces it once instead
  of every function re-checking.
- **Behaviour should differ by kind** — different algorithms, models or methods
  (`struct Euler end`, `struct RK4 end`) become types, and `solve(prob, ::Euler)` /
  `solve(prob, ::RK4)` dispatch on them instead of `if method == "euler"`. Define an
  abstract supertype (`abstract type AbstractIntegrator end`) as the extension point,
  so users can add their own.
- **Configuration with many options** — a keyword-constructed options struct (below)
  instead of a function with ten keyword arguments threaded through every call.
- **A `Dict{String,Any}` or `Dict{Symbol,Any}` holds fixed fields** — that's a struct
  in disguise, and it's slow (values are type `Any`) and typo-prone (`d["tol"]` vs
  `d["tolerance"]` fails at runtime, a misspelled field fails clearly).

Don't define one when a plain value or a standard type already says it all (a single
`Vector{Float64}` of samples), or for a one-off grouping used in a single place — a
`NamedTuple` `(; x, y)` is fine there and needs no definition.

### Designing structs

- Prefer immutable `struct` over `mutable struct`. Immutable structs can live on the
  stack or inline in arrays, and are easier to reason about. Use `mutable` only when
  fields genuinely change after construction and identity matters. To "change"
  an immutable struct, build a new one (Accessors.jl `@set` helps).
- Parametrize fields instead of using abstract types (see section 5, rule 3).
- Put invariants in an **inner constructor** so an invalid object can never exist:

  ```julia
  struct Interval{T<:Real}
      lo::T
      hi::T
      function Interval{T}(lo, hi) where {T<:Real}
          lo <= hi || throw(ArgumentError("lo must be ≤ hi, got $lo > $hi"))
          new{T}(lo, hi)
      end
  end
  Interval(lo::T, hi::T) where {T<:Real} = Interval{T}(lo, hi)
  Interval(lo::Real, hi::Real) = Interval(promote(lo, hi)...)
  ```

### Keyword constructors: `Base.@kwdef` vs `@with_kw`

Both generate a keyword-argument constructor with defaults, useful for configuration
and parameter structs with many fields.

- **`Base.@kwdef`** (built in, no dependency) — prefer it in packages:

  ```julia
  Base.@kwdef struct SolverOptions{T<:Real}
      tol::T = 1e-8
      maxiter::Int = 1_000
      verbose::Bool = false
  end
  SolverOptions(tol = 1e-6)
  ```

  It is exported as `@kwdef` from Julia 1.9 on; write `Base.@kwdef` for compatibility
  with older versions. It doesn't validate fields — add an outer check function or
  an inner constructor if you need that.

- **`@with_kw`** (Parameters.jl) adds extra features: defaults that refer to earlier
  fields, `@assert` invariants inside the struct body, a readable `show`, and
  `@unpack`/`@pack!` for pulling fields into local variables:

  ```julia
  using Parameters
  @with_kw struct Grid
      n::Int = 100
      L::Float64 = 1.0
      dx::Float64 = L / n        # depends on earlier fields
      @assert n > 0
  end
  @unpack n, dx = grid
  ```

  Use it when those features earn the extra dependency (common in scientific
  scripts and models). In new code, destructuring `(; n, dx) = grid` (Julia ≥ 1.7)
  replaces `@unpack`, so don't add Parameters.jl just for that.

- Either way, keep field types concrete or parametric — the macros don't fix
  `x::Real`-style slow fields.

## 5. Performance

Julia is fast when the compiler can infer concrete types. Most slow code violates one
of these:

1. **No non-constant globals in hot code.** Pass values as arguments, or declare them
   `const`. A non-const global has type `Any` everywhere it's used.
2. **Keep functions type-stable:** the return type should depend only on argument
   *types*, not values. `return x > 0 ? x : 0` is unstable for `Float64` input — use
   `zero(x)`. Likewise initialise accumulators with `zero(eltype(v))`, not `0`.
3. **Concrete struct fields.** `struct P; x::Real; end` is slow; use a type parameter:
   `struct P{T<:Real}; x::T; end`. Same for containers: `Vector{Float64}`, not
   `Vector{Real}` or `Vector{Any}`.
4. **Don't change a variable's type inside a function** (e.g. `x = 1` then `x = 1.5`).
5. **Function barriers:** when a type is only known at runtime (e.g. read from a file),
   pass the value into a separate inner function so that function compiles for the
   concrete type.
6. **Reduce allocations:** preallocate outputs and use `!` functions in loops; use
   `@views` or `view` for slices (slices copy by default); fuse broadcasts with
   `@.` / dots (`y .= a .* x .+ b`).
7. **Loop over arrays in column-major order** (first index innermost), or use
   `eachindex`/`CartesianIndices`.
8. **Closures:** a captured variable that's reassigned gets boxed and becomes
   type-unstable. Avoid reassigning it, or wrap with `let`.
9. `@inbounds` only after the indices are proven valid (e.g. via `eachindex`), and
   `@fastmath`/`@simd` only when reordering floating-point operations is acceptable.

### How to diagnose instead of guessing

- `@code_warntype f(args...)` — anything shown as `Any` or a red `Union` is a problem.
- JET.jl: `@report_opt f(args...)` finds runtime dispatch; `@report_call` finds likely
  errors.
- BenchmarkTools.jl: `@btime f($x)` / `@benchmark f($x)`. Interpolate globals with `$`,
  otherwise you're benchmarking global-variable access. The first call includes
  compilation, so never judge speed from one `@time` call.
- `@allocated f(x)` or the allocation count from `@btime` to check that a hot loop
  allocates nothing.
### Profiling: where the time actually goes

Run the function once first so compilation isn't profiled, then sample:

- **Quick look:** `@profview f(x)` in the VS Code Julia extension or with ProfileView.jl
  shows a flame graph. Red/yellow bars mark runtime dispatch and garbage collection.
- **PProf.jl** for deeper analysis (interactive flame graph, call graph, top-N table in
  the browser, and shareable profile files):

  ```julia
  using Profile, PProf
  f(x)                          # warm up / compile
  Profile.clear()
  @profile for _ in 1:100; f(x); end   # loop if one call is too short to sample
  pprof()                       # opens a web UI (default http://localhost:57599)
  ```

  Use `pprof(; web=false, out="profile.pb.gz")` to save without opening the browser
  (e.g. on a remote machine) and view later with `PProf.refresh(file="profile.pb.gz")`.
- **Allocation profiling** (Julia ≥ 1.8), to find which lines allocate:

  ```julia
  Profile.Allocs.clear()
  Profile.Allocs.@profile sample_rate=0.1 f(x)
  PProf.Allocs.pprof()
  ```

  Increase `sample_rate` (up to 1.0) for small workloads; it's slower but complete.

When the user asks for faster code, measure before and after and report the numbers.
Optimize the part the profiler points at, not what looks slow.

## 6. Latency and PackageCompiler.jl

**Only do this part when the user explicitly asks for it** (e.g. "build a sysimage",
"make startup faster", "compile this into an app", "add precompilation"). Don't add
PrecompileTools, PackageCompiler, sysimages or editor settings on your own initiative
while doing other work: they add build steps, large files and dependencies, and a
sysimage silently freezes package versions. If slow first-call latency is clearly
bothering the user, you may *mention* in one line that these options exist and ask
whether they want them — then wait for a yes.

Julia compiles on first call, so "time to first X" can be seconds. Once asked, choose
the lightest fix that solves the user's problem:

1. **PrecompileTools.jl** (inside the package) — cache compiled code for typical
   workloads at precompile time. This is the right fix for a *package's* latency,
   because every user benefits without doing anything:

   ```julia
   using PrecompileTools
   @setup_workload begin
       data = rand(10)
       @compile_workload begin
           process(data)
       end
   end
   ```

2. **PackageCompiler.jl** — for *users or deployments*, not for package code:
   - `create_sysimage([:Plots, :DataFrames]; sysimage_path="sys.so",
     precompile_execution_file="warmup.jl")` builds a custom system image with those
     packages baked in. Start Julia with `julia --sysimage sys.so`
     (`.dll` on Windows, `.dylib` on macOS). For VS Code, see below.
   - `create_app("MyApp", "MyAppCompiled")` builds a standalone executable; the
     package needs a `julia_main()::Cint` entry point.
   - `create_library(...)` builds a shared C library for calling Julia from other
     languages.

   Caveats to tell the user: a sysimage freezes the package versions in it, so it must
   be rebuilt after `Pkg.update` (and Revise can't change code baked into it);
   builds take minutes and the outputs are large (hundreds of MB); apps are built for
   one OS/architecture. The warm-up file should exercise real code paths, otherwise
   little gets precompiled.

### Using a sysimage in VS Code (`settings.json`)

This applies to **Visual Studio Code** with the Julia extension (the full Microsoft
Visual Studio IDE has no Julia support). Put project settings in the workspace file
`.vscode/settings.json` in the project root, so they apply only to that project;
user-wide settings go in the user `settings.json` (Command Palette → "Preferences:
Open User Settings (JSON)").

**Option A — let the extension build the sysimage** (simplest; the feature is marked
experimental):

1. Optionally add `.vscode/JuliaSysimage.toml` to control what's included:

   ```toml
   [sysimage]
   exclude = []            # packages from the environment NOT to bake in (e.g. ones you edit with Revise)
   statements_files = []   # precompile statement files, relative to the project folder
   execution_files = []    # warm-up scripts to run while building, relative to the project folder
   ```

2. Command Palette → "Tasks: Run Build Task" → "Julia: Build custom sysimage for
   current environment". It puts the sysimage in the environment folder.
3. Tell the extension to use it in `.vscode/settings.json`:

   ```json
   {
       "julia.environmentPath": "${workspaceFolder}",
       "julia.useCustomSysimage": true
   }
   ```

4. Restart the Julia REPL in VS Code ("Julia: Restart REPL").

**Option B — use a sysimage you built with PackageCompiler yourself:** pass it as a
startup argument:

```json
{
    "julia.environmentPath": "${workspaceFolder}",
    "julia.additionalArgs": ["--sysimage=${workspaceFolder}/sys.so"]
}
```

Then restart the REPL. Don't set both `useCustomSysimage` and a `--sysimage` argument.

Notes for either option:

- Rebuild the sysimage after changing or updating packages in the environment;
  otherwise the REPL keeps using the old baked-in versions.
- Don't bake in the package you're actively developing — Revise can't update code
  inside a sysimage. Exclude it (Option A) or leave it out of the package list (Option B).
- Add the sysimage file (`*.so`, `*.dll`, `*.dylib`) to `.gitignore`; it's large and
  machine-specific. `.vscode/settings.json` and `JuliaSysimage.toml` can be committed
  if the team shares the setup.
- Merge these keys into an existing `settings.json` rather than overwriting the file,
  and show the user the change before writing it.

## 7. Testing and quality checks

- Put tests in `test/runtests.jl` inside named `@testset`s; for larger packages,
  `include` one test file per source file.
- Test generic behaviour: at least one non-`Float64` element type (e.g. `Float32`,
  `Int`, a range or a view) for functions that claim to accept `AbstractVector`.
- Use `@test_throws` for error paths and `@inferred f(x)` for functions that must be
  type-stable.
- Add Aqua.jl (`Aqua.test_all(MyPkg)`) to catch missing compat entries, ambiguities,
  piracy, and stale deps.
- Run the suite with `julia --project -e 'using Pkg; Pkg.test()'` before saying work is
  done, and report the result honestly.

## 8. How to respond

- When writing code, apply these rules silently; mention the non-obvious choices
  briefly (e.g. "parametric field so it stays fast for Float32").
- When writing new code, design the types first: decide which structs (and abstract
  supertypes) the problem needs, then write methods on them (see section 4).
- When reviewing code, point out where a struct would make it clearer or faster
  (repeated argument groups, `Dict`s with fixed keys, string flags selecting an
  algorithm) and show the struct.
- When reviewing code, list concrete problems with the line, why it matters, and the
  fix. Lead with correctness, then performance, then style.
- Don't rewrite working code just for style unless asked; flag it instead.
