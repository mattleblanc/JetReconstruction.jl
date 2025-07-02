# defines and stores relevant variables to ghost creation and area calculation
struct GhostedArea
    resolution::Int
    n_ghosts::Int
    rapidity_step::Float64
    phi_step::Float64
    ghost_pt::Float64

    # constructor
    # default value of ghost_pt set to 1.0e-45 as -45 is the smallest magnitude that doesn't result in the value defaulting to 0
    function GhostedArea(resolution::Int)
        new(resolution,
            resolution * resolution,
            Float64(10.0 / resolution),
            Float64(6.28 / resolution),
            1.0e-45)
    end
end

# generates and adds ghosts to an event (a set of PseudoJets)
function add_ghosts!(event::Vector{PseudoJet}, gen::GhostedArea)
    # Set aside memory for output
    n_original = length(event)
    resize!(event, n_original + gen.n_ghosts)

    # Generate random values
    rand_vals = rand(Float64, 2 * gen.n_ghosts) .- 0.5f0

    # loop without array bounds checking in order to maximize time efficiency
    @inbounds for k in 0:(gen.n_ghosts - 1)
        i = k ÷ gen.resolution
        # modular division means j cycles through 0-99 100 times (acts like an inner loop)
        j = k % gen.resolution
        # index to iterate through the random values
        index = 2k + 1

        rap = (i + rand_vals[index]) * gen.rapidity_step
        phi = (j + rand_vals[index + 1]) * gen.phi_step

        # k + 1 is due to indexing beginning at 1 in julia
        event[n_original + k + 1] = PseudoJet(pt = gen.ghost_pt, rap = rap, phi = phi, _pure_ghost = true)
    end

    return event
end

# calculate the number of ghosts in a jet
function ghosts_in_jet(cluster_seq::ClusterSequence, jet::PseudoJet)

    # Get the constituents of the first jet
    jet_constituents = JetReconstruction.constituents(jet, cluster_seq)

    # determine the number of ghosts in the jet by checking each constituent's _pure_ghost field
    num_ghosts = 0
    for c in jet_constituents
        if c._pure_ghost
            num_ghost += 1
        end
    end

    # return the number of ghosts
    return num_ghosts
end

# calculate the area of a jet
function GhostedAreaCalculation(cluster_seq::ClusterSequence, jet::PseudoJet, resolution::Int)
    # the total area is the rapidity range (-5,5) times the phi range (0, 2pi)
    total_area = 10.0 * 6.28

    gen = GhostedArea(resolution)

    # ghost density is defined as the total number of ghosts divided by the total area
    ghost_density = gen.n_ghosts / total_area

    # determine the number of ghosts in the jet
    captured_ghosts = ghosts_in_jet(cluster_seq, jet)

    # area is equal to the number of ghosts in the jet divided by the density of ghosts
    return captured_ghosts / ghost_density
end

# returns a vector that contains the number of ghosts in each jet
function ghosts_in_jets(cluster_seq::ClusterSequence, jets::Vector{PseudoJet})

    # initialize a vector representing the number of ghosts in each jet
    ghosts_vector = Int[]

    # loop through each jet
    for jet in jets

        # Get the constituents of the first jet
        jet_constituents = JetReconstruction.constituents(jet, cluster_seq)

        # determine the number of ghosts in the jet by checking each constituent's _pure_ghost field
        num_ghosts = 0
        for c in jet_constituents
            if c._pure_ghost
                num_ghost += 1
            end
        end

        # add number of ghosts to the vector
        push!(ghosts_vector, num_ghosts)
    end

    # return a vector representing the number of ghosts in each jet
    return ghosts_vector
end

# returns a vector that contains the calculated area of each jets
function GhostedAreasCalculation(cluster_seq::ClusterSequence, jets::Vector{PseudoJet}, resolution::Int)
    # the total area is the rapidity range (-5,5) times the phi range (0, 2pi)
    total_area = 10.0 * 6.28

    gen = GhostedArea(resolution)

    # ghost density is defined as the total number of ghosts divided by the total area
    ghost_density = gen.n_ghosts / total_area

    # determine the number of ghosts in each jet
    captured_ghosts_vector = ghosts_in_jets(cluster_seq, jets)

    # area is equal to the number of ghosts in the jet divided by the density of ghosts
    areas_vector =  captured_ghosts_vector ./ ghost_density

    # return a vector representing the area of each jet
    return areas_vector
end