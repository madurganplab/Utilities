using Test

include(joinpath(@__DIR__, "..", "NuclearLevelScheme.jl"))
using .NuclearLevelScheme

@testset "NuclearLevelScheme" begin
    experiment = [
        0.0  0.0  +1.0
        1.2  2.0  +1.0
        2.0  2.0  -1.0
        2.4  3.5  -1.0
    ]
    theory = [
        0.0  0.0  +1.0
        1.3  2.0  +1.0
        1.8  2.0  -1.0
        2.6  3.5  -1.0
    ]

    @test level_scheme(experiment) isa NuclearLevelScheme.Plots.Plot
    @test level_scheme(experiment, theory) isa NuclearLevelScheme.Plots.Plot
    @test @level_scheme(experiment, theory) isa NuclearLevelScheme.Plots.Plot

    states = NuclearLevelScheme._validate_levels(theory, "theory")
    @test NuclearLevelScheme._yrast_indices(states) == [1, 2, 4]
    @test NuclearLevelScheme._experimental_match(states[2],
        NuclearLevelScheme._validate_levels(experiment, "experiment")) == 2
    @test NuclearLevelScheme._lowest_matching_state(
        NuclearLevelScheme._validate_levels(experiment, "experiment")[3],
        states,
    ) == 3

    unknown_assignments = [
        0.0  0        1
        1.0  missing  missing
        1.5  2        missing
        2.0  missing  -1
    ]
    unknown_states = NuclearLevelScheme._validate_levels(
        unknown_assignments,
        "experiment",
    )
    @test level_scheme(unknown_assignments) isa NuclearLevelScheme.Plots.Plot
    @test ismissing(unknown_states[2].spin)
    @test ismissing(unknown_states[2].parity)
    @test NuclearLevelScheme._format_spin_parity(missing, missing) == ""
    @test NuclearLevelScheme._format_spin_parity(2.0, missing) == "2"
    @test NuclearLevelScheme._format_spin_parity(missing, -1) == "⁻"
    @test NuclearLevelScheme._yrast_indices(unknown_states) == [1, 3]
    @test NuclearLevelScheme._experimental_match(states[2], unknown_states) === nothing

    @test_throws ArgumentError level_scheme(zeros(2, 2))
    @test_throws ArgumentError level_scheme(zeros(0, 3))
    @test_throws ArgumentError level_scheme([0.0 0.0 0.0])
    @test_throws ArgumentError level_scheme([0.0 0.25 1.0])
    @test_throws ArgumentError level_scheme([-1.0 0.0 1.0])
    @test_throws ArgumentError level_scheme(experiment, theory; connector_inset=-0.1)
end
