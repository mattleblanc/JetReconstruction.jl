include("common.jl")

@testset "Ghost Creation Tests" begin
    resolution = 1000
    ghosted_area = GhostedArea(resolution)

    # Create original event with preexisting PseudoJets
    event = [
        PseudoJet(1.1, 1.0, 1.3, 1.2),
        PseudoJet(2.4, 3.2, 0.8, 1.7),
        PseudoJet(0.9, 0.5, 3.0, 4.3),
        PseudoJet(0.9, 6.0, 5.2, 0.4)
    ]
    n_original = length(event)

    # Use the function being tested
    add_ghosts!(ghosted_area, event)

    # Ensure that the correct number of ghosts were generated
    @test length(event) == n_original + ghosted_area.n_ghosts

    # Ensure that the original jets have not changed
    @test event[1].px == 1.1
    @test event[2].py == 3.2
    @test event[3].pz == 3.0
    @test event[4].E == 0.4

    # Error bound for testing ghost pt values
    GHOST_PT2_THRESHOLD = 1.0e-89

    # Check that all added ghosts have the correct pt, use an error bound due to floating point error
    @test all(i -> abs(i._pt2 - 1.0e-90) <= GHOST_PT2_THRESHOLD, event[(n_original + 1):end])
end
