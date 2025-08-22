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

@testset "Ghosted Areas Calculation" begin
    data_file_path = "test/data/ghosted_area_example.hepmc.zst"
    
    # Parameters
    R = 0.4
    resolution = 100
    
    # Read the first event from the data file
    events = read_final_state_particles(data_file_path, maxevents = 1, skipevents = 0, T = PseudoJet)
    first_event = events[1]

    # Set up the GhostedArea structure and add ghosts
    ghosted_area = GhostedArea(resolution)
    add_ghosts!(ghosted_area, first_event)

    # Run jet reconstruction
    cluster_seq = jet_reconstruct(first_event, R = R, p = 1, algorithm = JetReconstruction.JetAlgorithm.Kt)
    clustered_jets = inclusive_jets(cluster_seq, ptmin = 5.0, T = PseudoJet)
    
    # Calculate the areas
    calculated_areas = ghosted_areas_calculation(ghosted_area, cluster_seq, clustered_jets)

    # Expected areas vector
    expected_areas = [0.58729773, 0.75188117, 0.99530458, 0.51896309, 0.50515609, 0.71486240, 0.29845130, 0.50935822, 
                      0.47263961, 0.50895802, 0.31195815, 0.29825120, 0.76198629, 0.92126704, 0.83932550]

    # Print the expected and calculated areas for comparison
    println("Expected areas:   ", expected_areas)
    println("Calculated areas: ", calculated_areas)

    # Check if the number of jets is what we expect
    @test length(calculated_areas) == length(expected_areas)

    # Check that each area is within 1% of the expected value
    for i in 1:length(expected_areas)
        @test (abs(calculated_areas[i] - expected_areas[i]) / abs(expected_areas[i])) * 100 <= 1

    end

end