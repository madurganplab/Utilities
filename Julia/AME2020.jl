module AME2020

export AMEQuantity,
       AMENuclide,
       AMETables,
       read_ame2020,
       ame2020,
       nuclide,
       mass,
       mass_excess,
       qbeta,
       q_beta_value,
       Qβ,
       sn,
       neutron_separation_energy,
       Sₙ,
       s2n,
       two_neutron_separation_energy,
       S₂ₙ

"""
    AMEQuantity

AME value with its uncertainty and whether the AME table marked it as estimated
with `#`. Energies are in keV unless a function docstring says otherwise.
"""
struct AMEQuantity
    value::Float64
    uncertainty::Float64
    estimated::Bool
end

"""
    AMENuclide

One nuclide from the AME2020 mass table.

`mass_excess`, `binding_energy_per_A`, and `qbeta` are in keV.
`atomic_mass` is in atomic mass units, u.
"""
struct AMENuclide
    N::Int
    Z::Int
    A::Int
    element::String
    mass_excess::AMEQuantity
    binding_energy_per_A::AMEQuantity
    qbeta::Union{Nothing, AMEQuantity}
    atomic_mass::AMEQuantity
end

struct AMEReactionRow
    A::Int
    Z::Int
    element::String
    values::NTuple{6, Union{Nothing, AMEQuantity}}
end

"""
    AMETables

Parsed AME2020 mass, `rct1`, and `rct2` tables.
"""
struct AMETables
    masses::Dict{Tuple{Int, Int}, AMENuclide}
    by_element_A::Dict{Tuple{String, Int}, Tuple{Int, Int}}
    rct1::Dict{Tuple{Int, Int}, AMEReactionRow}
    rct2::Dict{Tuple{Int, Int}, AMEReactionRow}
end

const _DEFAULT_TABLES = Ref{Union{Nothing, AMETables}}(nothing)

const _DEFAULT_AME_DIR = normpath(joinpath(@__DIR__, "..", "AME2020"))

@inline function _substr(line::AbstractString, i::Int, j::Int)
    i > lastindex(line) && return ""
    return line[i:min(j, lastindex(line))]
end

_parse_int_field(line::AbstractString, i::Int, j::Int) =
    tryparse(Int, strip(_substr(line, i, j)))

function _parse_quantity(line::AbstractString,
                         value_i::Int, value_j::Int,
                         unc_i::Int, unc_j::Int)
    value_string = _substr(line, value_i, value_j)
    unc_string = _substr(line, unc_i, unc_j)

    occursin("*", value_string) && return nothing
    isempty(strip(value_string)) && return nothing

    estimated = occursin("#", value_string) || occursin("#", unc_string)
    value = parse(Float64, replace(strip(value_string), '#'=>'.'))
    uncertainty = isempty(strip(unc_string)) ? NaN :
                  parse(Float64, replace(strip(unc_string), '#'=>'.'))

    return AMEQuantity(value, uncertainty, estimated)
end

function _read_mass_table(path::AbstractString)
    masses = Dict{Tuple{Int, Int}, AMENuclide}()
    by_element_A = Dict{Tuple{String, Int}, Tuple{Int, Int}}()

    open(path, "r") do io
        for line in eachline(io)
            N = _parse_int_field(line, 5, 9)
            Z = _parse_int_field(line, 10, 14)
            A = _parse_int_field(line, 15, 19)
            (N === nothing || Z === nothing || A === nothing) && continue

            element = strip(_substr(line, 21, 23))
            isempty(element) && continue

            mass_excess = _parse_quantity(line, 29, 42, 43, 54)
            binding_energy_per_A = _parse_quantity(line, 55, 67, 69, 78)
            mass_excess === nothing && continue
            binding_energy_per_A === nothing && continue

            qbeta = _parse_quantity(line, 82, 94, 95, 105)

            integer_part = parse(Float64, strip(_substr(line, 107, 109)))
            mass_micro = _parse_quantity(line, 111, 123, 124, typemax(Int))
            mass_micro === nothing && continue
            atomic_mass = AMEQuantity(
                integer_part + mass_micro.value * 1e-6,
                mass_micro.uncertainty * 1e-6,
                mass_micro.estimated,
            )

            item = AMENuclide(N, Z, A, element, mass_excess,
                              binding_energy_per_A, qbeta, atomic_mass)
            masses[(Z, N)] = item
            by_element_A[(uppercase(element), A)] = (Z, N)
        end
    end

    return masses, by_element_A
end

function _read_reaction_table(path::AbstractString)
    rows = Dict{Tuple{Int, Int}, AMEReactionRow}()

    open(path, "r") do io
        for line in eachline(io)
            A = _parse_int_field(line, 2, 4)
            Z = _parse_int_field(line, 9, 11)
            (A === nothing || Z === nothing) && continue

            element = strip(_substr(line, 6, 8))
            isempty(element) && continue

            values = ntuple(6) do i
                value_start = 13 + (i - 1) * 22
                _parse_quantity(line, value_start, value_start + 11,
                                value_start + 12, value_start + 21)
            end

            N = A - Z
            rows[(Z, N)] = AMEReactionRow(A, Z, element, values)
        end
    end

    return rows
end

"""
    read_ame2020(dir=joinpath(@__DIR__, "..", "AME2020")) -> AMETables

Parse the AME2020 `mass_1.mas20.txt`, `rct1.mas20.txt`, and
`rct2_1.mas20.txt` files. Energy quantities are returned in keV; atomic masses
are returned in u.
"""
function read_ame2020(dir::AbstractString=_DEFAULT_AME_DIR)
    masses, by_element_A = _read_mass_table(joinpath(dir, "mass_1.mas20.txt"))
    rct1 = _read_reaction_table(joinpath(dir, "rct1.mas20.txt"))
    rct2 = _read_reaction_table(joinpath(dir, "rct2_1.mas20.txt"))
    return AMETables(masses, by_element_A, rct1, rct2)
end

"""
    ame2020() -> AMETables

Return cached parsed AME2020 tables from this repository's `AME2020` folder.
"""
function ame2020()
    tables = _DEFAULT_TABLES[]
    tables === nothing || return tables

    tables = read_ame2020()
    _DEFAULT_TABLES[] = tables
    return tables
end

_require_quantity(q::Nothing, name::AbstractString, key) =
    throw(KeyError("AME2020 has no $name value for $key"))
_require_quantity(q::AMEQuantity, name::AbstractString, key) = q

function _key(tables::AMETables, Z::Integer, N::Integer)
    key = (Int(Z), Int(N))
    haskey(tables.masses, key) || throw(KeyError("AME2020 has no nuclide with Z=$Z, N=$N"))
    return key
end

function _key(tables::AMETables, element::AbstractString, A::Integer)
    lookup = (uppercase(strip(element)), Int(A))
    haskey(tables.by_element_A, lookup) ||
        throw(KeyError("AME2020 has no nuclide $(strip(element))-$A"))
    return tables.by_element_A[lookup]
end

"""
    nuclide(tables, Z, N)
    nuclide(tables, element, A)
    nuclide(Z, N)
    nuclide(element, A)

Return the `AMENuclide` record for a nuclide.
"""
nuclide(tables::AMETables, Z::Integer, N::Integer) = tables.masses[_key(tables, Z, N)]
nuclide(tables::AMETables, element::AbstractString, A::Integer) =
    tables.masses[_key(tables, element, A)]
nuclide(Z::Integer, N::Integer) = nuclide(ame2020(), Z, N)
nuclide(element::AbstractString, A::Integer) = nuclide(ame2020(), element, A)

"""
    mass(...) -> AMEQuantity

Return atomic mass in u.
"""
mass(tables::AMETables, Z::Integer, N::Integer) = nuclide(tables, Z, N).atomic_mass
mass(tables::AMETables, element::AbstractString, A::Integer) =
    nuclide(tables, element, A).atomic_mass
mass(Z::Integer, N::Integer) = mass(ame2020(), Z, N)
mass(element::AbstractString, A::Integer) = mass(ame2020(), element, A)

"""
    mass_excess(...) -> AMEQuantity

Return mass excess in keV.
"""
mass_excess(tables::AMETables, Z::Integer, N::Integer) = nuclide(tables, Z, N).mass_excess
mass_excess(tables::AMETables, element::AbstractString, A::Integer) =
    nuclide(tables, element, A).mass_excess
mass_excess(Z::Integer, N::Integer) = mass_excess(ame2020(), Z, N)
mass_excess(element::AbstractString, A::Integer) = mass_excess(ame2020(), element, A)

"""
    qbeta(...) -> AMEQuantity

Return beta-decay Q value in keV.
"""
qbeta(tables::AMETables, Z::Integer, N::Integer) =
    _require_quantity(nuclide(tables, Z, N).qbeta, "Qβ", (Z=Z, N=N))
qbeta(tables::AMETables, element::AbstractString, A::Integer) =
    _require_quantity(nuclide(tables, element, A).qbeta, "Qβ", (element=element, A=A))
qbeta(Z::Integer, N::Integer) = qbeta(ame2020(), Z, N)
qbeta(element::AbstractString, A::Integer) = qbeta(ame2020(), element, A)

q_beta_value(args...) = qbeta(args...)
Qβ(args...) = qbeta(args...)

function _reaction_quantity(table::Dict{Tuple{Int, Int}, AMEReactionRow},
                            key::Tuple{Int, Int},
                            column::Int,
                            name::AbstractString)
    haskey(table, key) || throw(KeyError("AME2020 has no $name row for Z=$(key[1]), N=$(key[2])"))
    return _require_quantity(table[key].values[column], name, (Z=key[1], N=key[2]))
end

"""
    sn(...) -> AMEQuantity

Return one-neutron separation energy, S(n), in keV from `rct2`.
"""
sn(tables::AMETables, Z::Integer, N::Integer) =
    _reaction_quantity(tables.rct2, _key(tables, Z, N), 1, "S(n)")
sn(tables::AMETables, element::AbstractString, A::Integer) =
    _reaction_quantity(tables.rct2, _key(tables, element, A), 1, "S(n)")
sn(Z::Integer, N::Integer) = sn(ame2020(), Z, N)
sn(element::AbstractString, A::Integer) = sn(ame2020(), element, A)

neutron_separation_energy(args...) = sn(args...)
Sₙ(args...) = sn(args...)

"""
    s2n(...) -> AMEQuantity

Return two-neutron separation energy, S(2n), in keV from `rct1`.
"""
s2n(tables::AMETables, Z::Integer, N::Integer) =
    _reaction_quantity(tables.rct1, _key(tables, Z, N), 1, "S(2n)")
s2n(tables::AMETables, element::AbstractString, A::Integer) =
    _reaction_quantity(tables.rct1, _key(tables, element, A), 1, "S(2n)")
s2n(Z::Integer, N::Integer) = s2n(ame2020(), Z, N)
s2n(element::AbstractString, A::Integer) = s2n(ame2020(), element, A)

two_neutron_separation_energy(args...) = s2n(args...)
S₂ₙ(args...) = s2n(args...)

end
