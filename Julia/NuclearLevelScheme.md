# NuclearLevelScheme

`NuclearLevelScheme.jl` draws experimental nuclear level schemes and compares
them with optional theoretical calculations using Plots.jl.

Each state is shown as a horizontal line. Its excitation energy is printed
above the left end, and its spin-parity assignment is printed above the right
end. When theoretical states are supplied, they are drawn in a second column.

## Requirements

- Julia
- [Plots.jl](https://docs.juliaplots.org/)

Install Plots.jl from the Julia REPL if necessary:

```julia
using Pkg
Pkg.add("Plots")
```

## Loading the module

From the repository root:

```julia
include("Julia/NuclearLevelScheme.jl")
using .NuclearLevelScheme
```

The module exports:

- `level_scheme`
- `@level_scheme`

## Level matrix format

A level scheme may be an `nstates × 3` or `nstates × 4` matrix with one state
per row. The recommended four-column format is:

```julia
levels = [
    0.0  0    +1  :firm
    1.2  2    +1  :tentative
    2.4  7/2  -1  missing
]
```

The columns are:

1. Excitation energy
2. Spin
3. Parity
4. Assignment status

Excitation energies must be finite, non-negative real numbers. Spins must be
non-negative integers, half-integers, or ranges containing those values. Parity
must be `+1` or `-1`.

The module prints parity as a superscript:

- `+1` becomes `⁺`
- `-1` becomes `⁻`

Assignment status must be one of:

- `:firm`: display the known assignment normally, such as `4⁺`.
- `:tentative`: display the assignment in parentheses, such as `(4⁺)`.
- `missing`: treat the assignment as unknown, show no spin-parity label, and
  exclude the state from theory matching.

Three-column matrices remain supported. Their assignments are treated as
`:firm` for backward compatibility.

## Spin ranges

Use a Julia range when several consecutive spins are possible:

```julia
levels = Matrix{Any}(undef, 2, 4)
levels[1, :] = Any[0.0, 0,   +1, :firm]
levels[2, :] = Any[1.2, 2:4, +1, :firm]
```

The second assignment is displayed as `(2⁺, 3⁺, 4⁺)`. A spin range is always
shown as tentative and parenthesized, even when its fourth-column status is
`:firm`.

Half-integer ranges are also supported:

```julia
levels = Matrix{Any}(undef, 1, 4)
levels[1, :] = Any[1.8, 3//2:1:7//2, -1, :tentative]
```

This displays `(3/2⁻, 5/2⁻, 7/2⁻)`. Every value must be non-negative and
integral or half-integral. During theory comparison, the experimental level is
connected to the lowest matching theoretical state for each possible spin in
its range.

The explicit `Matrix{Any}` construction is necessary because Julia expands a
range placed directly inside a matrix literal instead of storing the range as
one cell.

The energy unit is not imposed by the module. Use the same unit consistently
for every energy in both matrices.

## Experimental level scheme

```julia
experiment = [
    0.0  0    +1  :firm
    1.2  2    +1  :firm
    2.0  2    -1  :firm
    2.4  7/2  -1  :tentative
]

plot = level_scheme(experiment)
```

The default title beneath the column is `Experiment`. Supply a different title
with `experiment_title`:

```julia
plot = level_scheme(experiment; experiment_title="Measured levels")
```

## Unknown experimental assignments

Use `missing` assignment status when the complete experimental assignment is
unknown:

```julia
experiment = [
    0.0  0    +1  :firm
    0.8  0    +1  missing
    1.2  2    +1  :tentative
    2.4  7/2  -1  :firm
]
```

All these states are drawn. The `missing` row has no spin-parity label and is
not connected to theory. Its spin and parity cells may themselves also be
`missing`:

```julia
0.8  missing  missing  missing
```

Because matrices containing `missing` cannot be plain numeric matrices, Julia
will construct a matrix capable of holding both numbers and `missing`. No
manual element-type declaration is normally required when using the matrix
literal shown above.

## Experimental and theoretical schemes

Pass the theoretical matrix as the second positional argument:

```julia
experiment = [
    0.0  0    +1
    1.2  2    +1
    2.0  2    -1
    2.4  7/2  -1
]

theory = [
    0.0  0    +1
    1.3  2    +1
    1.8  2    -1
    2.6  7/2  -1
]

plot = level_scheme(
    experiment,
    theory;
    experiment_title="Experiment",
    theory_title="Shell model",
)
```

Every experimental state with a known spin and parity is connected to the
lowest-energy theoretical state having the same spin and parity. An
experimental state is not connected when:

- its spin is missing,
- its parity is missing, or
- its assignment status is missing, or
- no theoretical state has the same spin and parity.

Both `:firm` and `:tentative` assignments participate in matching.

If several experimental states share a spin-parity assignment, each is
connected to the same lowest-energy matching theoretical state.

## Titles from variable names

The macro form uses the input variable names as the titles beneath the columns:

```julia
@level_scheme experiment theory
```

This is equivalent to:

```julia
level_scheme(
    experiment,
    theory;
    experiment_title="experiment",
    theory_title="theory",
)
```

Use the function form when keyword customization is needed.

## Saving the figure

The returned object is a Plots.jl plot and can be saved in any format supported
by the active Plots backend:

```julia
scheme = level_scheme(experiment, theory; theory_title="Shell model")
NuclearLevelScheme.Plots.savefig(scheme, "level_scheme.png")
```

If `Plots` is also imported into the calling scope, the usual form works:

```julia
using Plots
savefig(scheme, "level_scheme.pdf")
```

## Appearance options

The plotting function supports these module-specific keyword arguments:

| Keyword | Default | Purpose |
|---|---:|---|
| `experiment_title` | `"Experiment"` | Title beneath the experimental column |
| `theory_title` | `"Theory"` | Title beneath the theoretical column |
| `line_length` | `0.72` | Length of every state line |
| `energy_digits` | `3` | Maximum decimal digits shown for energies |
| `level_linewidth` | `4` | Width of experimental and theoretical state lines |
| `connector_linewidth` | `1.25` | Width of matching lines |
| `connector_inset` | `0.08` | Clearance between connectors and state columns |
| `connector_color` | `:gray45` | Matching-line color |
| `level_color` | `:black` | Experimental state color |
| `theory_color` | `:royalblue3` | Theoretical state color |
| `size` | `(900, 650)` | Figure dimensions in pixels |

The built-in font sizes are 15 points for energies and 17 points for
spin-parity assignments and column titles.

Additional keyword arguments are forwarded to `Plots.plot`:

```julia
scheme = level_scheme(
    experiment,
    theory;
    level_linewidth=5,
    connector_linewidth=1,
    connector_inset=0.12,
    level_color=:black,
    theory_color=:darkorange,
    size=(1100, 750),
    background_color=:white,
)
```

A larger `connector_inset` makes connecting lines shorter. The value must be
non-negative.

## Input errors

`level_scheme` throws an `ArgumentError` for:

- an input that is not an `nstates × 3` or `nstates × 4` matrix,
- an empty matrix,
- a missing, negative, infinite, or nonnumeric excitation energy,
- a known spin that is negative, infinite, or not integer/half-integer,
- an empty spin range or a range containing an invalid spin,
- a known parity other than `+1` or `-1`, or
- a status other than `:firm`, `:tentative`, or `missing`,
- invalid appearance values such as a negative connector inset.

## Running the tests

From the repository root:

```sh
julia --project=. Julia/test/NuclearLevelSchemeTests.jl
```
