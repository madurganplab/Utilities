module NuclearLevelScheme

using Plots

export level_scheme, @level_scheme

"""
    level_scheme(experiment; experiment_title="Experiment", kwargs...)
    level_scheme(experiment, theory; experiment_title="Experiment",
                 theory_title="Theory", kwargs...)

Plot a nuclear level scheme from `nstates × 3` or `nstates × 4` matrices. The
first three columns contain `(excitation_energy, spin, parity)`, where parity is
numerically `+1` or `-1`. In a four-column matrix, assignment status is `:firm`,
`:tentative`, or `missing`. Tentative assignments are displayed in parentheses;
missing-status assignments are displayed without a spin-parity label.
Spin may be integral, half-integral, or a non-empty range of such values. A
range such as `2:4` is always displayed as tentative: `(2⁺, 3⁺, 4⁺)`. Use
`missing` for an unknown spin or parity; the level is still drawn and its known
quantum-number information is shown. Energies must be finite, non-negative real
numbers, and known spins must also be finite and non-negative.

When `theory` is supplied, it is drawn beside the experimental scheme. Every
experimental state with known spin and parity is connected to the lowest-energy
theoretical state with the same spin and parity. An experimental state with a
missing assignment, or without a matching theoretical state, is not connected.

Keyword arguments `line_length`, `energy_digits`, `level_linewidth`,
`connector_linewidth`, `connector_inset`, and other Plots.jl keyword arguments
may be used to customize the output. The returned value is a `Plots.Plot`.

# Examples

```julia
experiment = [
    0.0  0  +1
    1.2  2  +1
    2.1  3  -1
]

theory = [
    0.0  0  +1
    1.1  2  +1
    1.8  2  -1
    2.3  3  -1
]

level_scheme(experiment, theory; theory_title="Shell model")
```

Unknown experimental assignments can be entered directly:

```julia
experiment = [
    0.0  0        +1
    1.2  missing  missing
    2.1  3        missing
]
```

Use [`@level_scheme`](@ref) to derive the bottom labels from variable names.
"""
function level_scheme(
    experiment::AbstractMatrix;
    experiment_title::AbstractString="Experiment",
    kwargs...,
)
    experimental_states = _validate_levels(experiment, "experiment")
    return _draw_scheme(
        experimental_states,
        nothing;
        experiment_title,
        theory_title="",
        kwargs...,
    )
end

function level_scheme(
    experiment::AbstractMatrix,
    theory::AbstractMatrix;
    experiment_title::AbstractString="Experiment",
    theory_title::AbstractString="Theory",
    kwargs...,
)
    experimental_states = _validate_levels(experiment, "experiment")
    theoretical_states = _validate_levels(theory, "theory")
    return _draw_scheme(
        experimental_states,
        theoretical_states;
        experiment_title,
        theory_title,
        kwargs...,
    )
end

"""
    @level_scheme experiment [theory]

Plot level matrices and use their variable names as the titles beneath the
schemes. For expressions without a simple variable name, the expression text is
used as the title. Use `level_scheme` directly when keyword customization is
needed.
"""
macro level_scheme(experiment)
    experiment_title = string(experiment)
    return :(level_scheme(
        $(esc(experiment));
        experiment_title=$experiment_title,
    ))
end

macro level_scheme(experiment, theory)
    experiment_title = string(experiment)
    theory_title = string(theory)
    return :(level_scheme(
        $(esc(experiment)),
        $(esc(theory));
        experiment_title=$experiment_title,
        theory_title=$theory_title,
    ))
end

function _validate_levels(levels::AbstractMatrix, name::AbstractString)
    size(levels, 2) in (3, 4) || throw(ArgumentError(
        "$name levels must be an nstates × 3 or nstates × 4 matrix; got size $(size(levels))",
    ))
    size(levels, 1) > 0 || throw(ArgumentError("$name levels must not be empty"))

    states = NamedTuple{
        (:energy, :spin, :parity, :status),
        Tuple{Float64,Union{Missing,Float64,Vector{Float64}},Union{Missing,Int},Union{Missing,Symbol}},
    }[]
    for row in axes(levels, 1)
        energy, spin, parity = levels[row, 1], levels[row, 2], levels[row, 3]
        status = size(levels, 2) == 4 ? levels[row, 4] : :firm

        energy isa Real || throw(ArgumentError(
            "$name energy in row $row must be a real number",
        ))
        isfinite(energy) && energy >= 0 || throw(ArgumentError(
            "$name energy in row $row must be finite and non-negative",
        ))
        normalized_spin = _validate_spin(spin, name, row)
        if !ismissing(parity)
            (parity isa Real && (parity == 1 || parity == -1)) || throw(ArgumentError(
                "$name parity in row $row must be +1, -1, or missing",
            ))
        end
        (ismissing(status) || status === :firm || status === :tentative) || throw(ArgumentError(
            "$name status in row $row must be :firm, :tentative, or missing",
        ))

        push!(states, (
            energy=Float64(energy),
            spin=normalized_spin,
            parity=ismissing(parity) ? missing : Int(parity),
            status=status,
        ))
    end
    return states
end

function _yrast_indices(states)
    lowest = Dict{Float64,Int}()
    for (index, state) in pairs(states)
        (ismissing(state.spin) || ismissing(state.status)) && continue
        for spin in _spin_values(state.spin)
            if !haskey(lowest, spin) || state.energy < states[lowest[spin]].energy
                lowest[spin] = index
            end
        end
    end
    return unique(sort!(collect(values(lowest)); by=index -> states[index].energy))
end

function _experimental_match(theoretical_state, experimental_states)
    return _lowest_matching_state(theoretical_state, experimental_states)
end

function _lowest_matching_state(reference_state, candidate_states)
    matches = _matching_state_indices(reference_state, candidate_states)
    isempty(matches) && return nothing
    return matches[1]
end

function _matching_state_indices(reference_state, candidate_states)
    (ismissing(reference_state.spin) || ismissing(reference_state.parity) ||
     ismissing(reference_state.status)) && return Int[]
    matches = Int[]
    for spin in _spin_values(reference_state.spin)
        candidates = findall(state ->
            !ismissing(state.spin) &&
            !ismissing(state.parity) &&
            !ismissing(state.status) &&
            spin in _spin_values(state.spin) &&
            state.parity == reference_state.parity,
            candidate_states,
        )
        isempty(candidates) && continue
        push!(matches, candidates[argmin([candidate_states[index].energy for index in candidates])])
    end
    return unique(matches)
end

function _draw_scheme(
    experiment,
    theory;
    experiment_title::AbstractString,
    theory_title::AbstractString,
    line_length::Real=0.72,
    energy_digits::Integer=3,
    level_linewidth::Real=4,
    connector_linewidth::Real=1.25,
    connector_inset::Real=0.08,
    connector_color=:gray45,
    level_color=:black,
    theory_color=:royalblue3,
    size=(900, 650),
    kwargs...,
)
    line_length > 0 || throw(ArgumentError("line_length must be positive"))
    energy_digits >= 0 || throw(ArgumentError("energy_digits must be non-negative"))
    connector_inset >= 0 || throw(ArgumentError("connector_inset must be non-negative"))

    all_states = theory === nothing ? experiment : vcat(experiment, theory)
    maximum_energy = maximum(state.energy for state in all_states)
    energy_span = max(maximum_energy, 1.0)
    bottom_margin = 0.14 * energy_span
    top_margin = 0.12 * energy_span
    label_y_offset = 0.025 * energy_span

    experiment_center = theory === nothing ? 0.0 : -0.75
    theory_center = 0.75
    half_length = line_length / 2

    plot_object = plot(
        ;
        xlims=theory === nothing ? (-0.75, 0.75) : (-1.4, 1.4),
        ylims=(-bottom_margin, maximum_energy + top_margin),
        axis=false,
        ticks=false,
        grid=false,
        legend=false,
        framestyle=:none,
        size,
        kwargs...,
    )

    _draw_levels!(
        plot_object,
        experiment,
        experiment_center,
        half_length,
        label_y_offset,
        energy_digits,
        level_color,
        level_linewidth,
    )

    if theory !== nothing
        _draw_levels!(
            plot_object,
            theory,
            theory_center,
            half_length,
            label_y_offset,
            energy_digits,
            theory_color,
            level_linewidth,
        )

        for experiment_state in experiment
            for theory_index in _matching_state_indices(experiment_state, theory)
                plot!(
                    plot_object,
                    [
                        experiment_center + half_length + connector_inset,
                        theory_center - half_length - connector_inset,
                    ],
                    [experiment_state.energy, theory[theory_index].energy];
                    color=connector_color,
                    linewidth=connector_linewidth,
                    label=false,
                )
            end
        end
    end

    title_y = -0.72 * bottom_margin
    annotate!(plot_object, experiment_center, title_y, text(experiment_title, 17, :center))
    if theory !== nothing
        annotate!(plot_object, theory_center, title_y, text(theory_title, 17, :center))
    end

    return plot_object
end

function _draw_levels!(
    plot_object,
    states,
    center,
    half_length,
    label_y_offset,
    energy_digits,
    color,
    linewidth,
)
    for state in states
        left = center - half_length
        right = center + half_length
        plot!(
            plot_object,
            [left, right],
            [state.energy, state.energy];
            color,
            linewidth,
            label=false,
        )
        annotate!(
            plot_object,
            left,
            state.energy + label_y_offset,
            text(_format_energy(state.energy, energy_digits), 15, :left),
        )
        annotate!(
            plot_object,
            right,
            state.energy + label_y_offset,
            text(_format_spin_parity(state.spin, state.parity, state.status), 17, :right),
        )
    end
    return plot_object
end

function _format_energy(energy::Real, digits::Integer)
    rounded = round(energy; digits)
    return isinteger(rounded) ? string(Int(rounded)) : string(rounded)
end

function _format_spin_parity(spin, parity, status=:firm)
    ismissing(status) && return ""
    spin_labels = if ismissing(spin)
        String[]
    else
        [_format_spin(value) for value in _spin_values(spin)]
    end
    parity_label = ismissing(parity) ? "" : parity == 1 ? "⁺" : "⁻"
    labels = isempty(spin_labels) ? [parity_label] : [label * parity_label for label in spin_labels]
    label = join(labels, ", ")
    tentative = spin isa AbstractVector || status === :tentative
    return tentative && !isempty(label) ? "($label)" : label
end

function _format_spin(spin::Real)
    doubled_spin = round(Int, 2 * spin)
    return iseven(doubled_spin) ? string(doubled_spin ÷ 2) : "$(doubled_spin)/2"
end

_spin_values(spin::Real) = (spin,)
_spin_values(spin::AbstractVector) = spin

function _validate_spin(spin, name::AbstractString, row)
    if ismissing(spin)
        return missing
    elseif spin isa AbstractRange
        isempty(spin) && throw(ArgumentError("$name spin range in row $row must not be empty"))
        values = Float64[]
        for value in spin
            _validate_spin_value(value, name, row)
            push!(values, Float64(value))
        end
        return values
    else
        _validate_spin_value(spin, name, row)
        return Float64(spin)
    end
end

function _validate_spin_value(spin, name::AbstractString, row)
    spin isa Real || throw(ArgumentError(
        "$name spin in row $row must be a real number, range, or missing",
    ))
    isfinite(spin) && spin >= 0 || throw(ArgumentError(
        "$name spin in row $row must be finite and non-negative",
    ))
    isinteger(2 * spin) || throw(ArgumentError(
        "$name spin in row $row must be integral or half-integral",
    ))
    return nothing
end

end
