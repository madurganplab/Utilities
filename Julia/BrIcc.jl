module BrIcc

export conversion_coefficients, BrIccResult

"""Result returned by `conversion_coefficients`.

`values` and `uncertainties` contain the numeric XML fields returned by BrIccS,
keyed by their XML path (for example `"Tot/ICC"`). `xml` is retained for
fields that are not represented by the lightweight parser.
"""
struct BrIccResult
    values::Dict{String,Float64}
    uncertainties::Dict{String,Float64}
    xml::String
end

const _DEFAULT_EXE = "/Applications/briccs"
const _DATA_DIR = normpath(joinpath(@__DIR__, "..", "BrIccData"))

function _number(s)
    t = strip(s)
    isempty(t) && return nothing
    t in ("*", "-", "--", "N/A") && return nothing
    try
        parse(Float64, replace(t, 'D' => 'E', 'd' => 'e'))
    catch
        nothing
    end
end

function _xml_numbers(xml::String)
    values = Dict{String,Float64}()
    uncertainties = Dict{String,Float64}()
    for m in eachmatch(r"(?s)<PureCC\s+Shell=\"([^\"]+)\"[^>]*>([^<]+)</PureCC>", xml)
        n = _number(m.captures[2])
        n === nothing || (values[m.captures[1]] = n)
    end
    # BrIccS emits simple XML. Capture leaf tags and attributes without adding
    # an XML package dependency to this small wrapper.
    for m in eachmatch(r"<([A-Za-z][\w:.-]*)(?:\s+[^>]*)?>([^<]+)</\1>", xml)
        key, val = m.captures
        key == "PureCC" && continue
        n = _number(val)
        n === nothing && continue
        occursin(r"(?i)(uncert|error|sigma|delta|d[a-z])", key) ?
            (uncertainties[key] = n) : (values[key] = n)
    end
    return values, uncertainties
end

"""Run BrIccS and return conversion coefficients.

`energy` is in keV, `Z` is the atomic number, and `multipolarity` is e.g.
`"M1"` or `"M1+E2"`. BrIccS must be installed at `/Applications/briccs`,
and its data files must be in the bundled `BrIccData` directory.
Use `all_subshells=true` to request the complete shell listing.
"""
function conversion_coefficients(Z::Integer, energy::Real;
        multipolarity::AbstractString="M1", mixing_ratio=nothing,
        energy_uncertainty=nothing, mixing_uncertainty=nothing,
        dataset::AbstractString="BrIccFO", all_subshells::Bool=false)
    5 <= Z <= 110 || throw(ArgumentError("Z must be between 5 and 110"))
    energy > 0 || throw(ArgumentError("energy must be positive (keV)"))
    executable = _DEFAULT_EXE
    isfile(executable) || throw(ArgumentError("BrIccS executable not found at $executable"))
    isdir(_DATA_DIR) || throw(ArgumentError("BrIcc data directory not found: $_DATA_DIR"))
    required = dataset == "BrIccNH" ?
        ("BrIccNHV22.idx", "BrIccNHV22.icc") :
        ("BrIccFOV22.idx", "BrIccFOV22.icc")
    missing = filter(name -> !isfile(joinpath(_DATA_DIR, name)), required)
    isempty(missing) || throw(ArgumentError(
        "Required BrIcc data file(s) missing from $_DATA_DIR: " * join(missing, ", ")))
    dir = _DATA_DIR

    args = ["-Z", string(Z), "-g", string(energy), "-L", String(multipolarity), "-w", dataset]
    mixing_ratio !== nothing && append!(args, ["-d", string(mixing_ratio)])
    energy_uncertainty !== nothing && append!(args, ["-e", string(energy_uncertainty)])
    mixing_uncertainty !== nothing && append!(args, ["-u", string(mixing_uncertainty)])
    all_subshells && push!(args, "-a")

    # The supplied macOS binary may lack its executable bit; stage it safely.
    tmp = mktempdir()
    staged = joinpath(tmp, "briccs")
    cp(executable, staged; force=true)
    chmod(staged, 0o700)
    try
        io = IOBuffer()
        run(pipeline(setenv(`$staged $args`, "BrIccHome" => dir), stdout=io, stderr=io))
        xml = String(take!(io))
        occursin(r"<F>|<ERROR>|<SEVERE>", xml) &&
            throw(ErrorException("BrIccS failed: " * replace(strip(xml), '\n' => ' ')))
        values, uncertainties = _xml_numbers(xml)
        return BrIccResult(values, uncertainties, xml)
    finally
        rm(tmp; recursive=true, force=true)
    end
end

end
