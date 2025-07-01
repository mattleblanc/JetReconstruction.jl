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
        event[n_original + k + 1] = PseudoJet(pt = gen.ghost_pt, rap = rap, phi = phi)
    end

    return event
end

# calculated the number of ghosts in a jet
function ghosts_in_jet(cluster_seq::ClusterSequence, filter::Float64 = 1.0e-89)
    # Retrieve the exclusive jets, but as `PseudoJet` types
    jets = inclusive_jets(cluster_seq; ptmin = 5.0, T = PseudoJet)

    # Get the constituents of the first jet
    jet_constituents = JetReconstruction.constituents(jets[1], cluster_seq)

    # determine the number of ghosts through checking each constituent's pt2
    num_ghosts = 0
    for c in jet_constituents
        if c._pt2 < filter
            num_ghost += 1
        end
    end

    # return the number of ghosts
    return num_ghosts
end

# calculated the area of a jet
function GhostedAreaCalculation(event::Vector{PseudoJet}, resolution::int, cluster_seq::ClusterSequence)
    # the total area is the rapidity range (-5,5) times the phi range (0, 2pi)
    total_area = 10.0 * 6.28

    gen = GhostedArea(resolution)

    # ghost density is defined as the total number of ghosts divided by the total area
    ghost_density = gen.n_ghosts / total_area

    # add the ghosts to the event, and determine the number of ghosts in the jet
    add_ghosts!(event, gen)
    captured_ghosts = ghosts_in_jet(cluster_seq)

    # area is equal to the number of ghosts in the jet divided by the density of ghosts
    return captured_ghosts / ghost_density
end