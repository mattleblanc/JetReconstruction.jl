# defines and stores relevant variables to ghost creation and area calculation
struct GhostedArea
    resolution::Int
    n_ghosts::Int
    rap_max::Float64
    rap_min::Float64
    rapidity_step::Float64
    phi_step::Float64
    ghost_pt::Float64
    ghost_density::Float64

    # constructor
    # default value of ghost_pt set to 1.0e-45 as -45 is the smallest magnitude that doesn't result in the value defaulting to 0
    function GhostedArea(resolution::Int; rap_max::Float64 = 5.0, rap_min::Float64 = -5.0, ghost_pt::Float64 = 1.0e-45)
        new(resolution,
            resolution * resolution,
            rap_max,
            rap_min,
            abs(rap_max - rap_min) / resolution,
            2π / resolution,
            ghost_pt,
            (resolution * resolution) / (abs(rap_max - rap_min) * 2π))
    end
end

"""
    add_ghosts!(ghosted_area::GhostedArea, event::Vector{PseudoJet})

Adds ghost particles to a single event. Ghost particles are pseudo-particles
used for jet area calculations.

# Arguments
- `ghosted_area::GhostedArea`: The ghosted area configuration.
- `event::Vector{PseudoJet}`: The event to which ghosts will be added.

# Returns
The modified event with added ghost particles.
"""
function add_ghosts!(ghosted_area::GhostedArea, event::Vector{PseudoJet})
    # Set aside memory for output
    n_original = length(event)
    resize!(event, n_original + ghosted_area.n_ghosts)

    # Generate random values
    rand_vals = rand(Float64, 2 * ghosted_area.n_ghosts) .- 0.5

    # loop without array bounds checking in order to maximize time efficiency
    @inbounds for k in 0:(ghosted_area.n_ghosts - 1)
        i = k ÷ ghosted_area.resolution
        # modular division means j cycles through 0-99 100 times (acts like an inner loop)
        j = k % ghosted_area.resolution
        # index to iterate through the random values
        index = 2k + 1

        rap = ghosted_area.rap_min + (i + rand_vals[index]) * ghosted_area.rapidity_step
        phi = (j + rand_vals[index + 1]) * ghosted_area.phi_step

        # k + 1 is due to indexing beginning at 1 in julia
        event[n_original + k + 1] = PseudoJet(pt = ghosted_area.ghost_pt, rap = rap, phi = phi, cluster_hist_index = n_original + k + 1, pure_ghost = true)
    end

    return event
end

"""
    add_ghosts!(ghosted_area::GhostedArea, events::Vector{Vector{PseudoJet}})

Adds ghost particles to multiple events.

# Arguments
- `ghosted_area::GhostedArea`: The ghosted area configuration.
- `events::Vector{Vector{PseudoJet}}`: A collection of events to which ghosts will be added.

# Returns
A collection of modified events with added ghost particles.
"""
function add_ghosts!(ghosted_area::GhostedArea, events::Vector{Vector{PseudoJet}})
    return map(event -> add_ghosts!(ghosted_area, event), events)
end

"""
    ghosts_in_jet(cluster_seq::ClusterSequence, jet::PseudoJet)

Calculates the number of ghost particles in a given jet.

# Arguments
- `cluster_seq::ClusterSequence`: The clustering sequence of the event.
- `jet::PseudoJet`: The jet for which the ghost count is calculated.

# Returns
The number of ghost particles in the jet.
"""
function ghosts_in_jet(cluster_seq::ClusterSequence, jet::PseudoJet)
    # Get the constituents of the jet
    jet_constituents = JetReconstruction.constituents(jet, cluster_seq)

    # return the number of ghosts in the jet by checking each constituent's _pure_ghost field
    return count(is_pure_ghost, jet_constituents)
end

"""
    ghosted_area_calculation(ghosted_area::GhostedArea, cluster_seq::ClusterSequence, jet::PseudoJet)

Calculates the area of a jet based on the number of ghost particles it contains.

# Arguments
- `ghosted_area::GhostedArea`: The ghosted area configuration.
- `cluster_seq::ClusterSequence`: The clustering sequence of the event.
- `jet::PseudoJet`: The jet for which the area is calculated.

# Returns
The calculated area of the jet.
"""
function ghosted_area_calculation(ghosted_area::GhostedArea, cluster_seq::ClusterSequence, jet::PseudoJet)
    # area is equal to the number of ghosts in the jet divided by the density of ghosts
    return ghosts_in_jet(cluster_seq, jet) ./ ghosted_area.ghost_density
end

"""
    ghosts_in_jets(cluster_seq::ClusterSequence, jets::Vector{PseudoJet})

Calculates the number of ghost particles in each jet in a collection of jets.

# Arguments
- `cluster_seq::ClusterSequence`: The clustering sequence of the event.
- `jets::Vector{PseudoJet}`: A collection of jets.

# Returns
A vector containing the number of ghost particles in each jet.
"""
function ghosts_in_jets(cluster_seq::ClusterSequence, jets::Vector{PseudoJet})
    return map(jet -> ghosts_in_jet(cluster_seq, jet), jets)
end

"""
    ghosted_areas_calculation(ghosted_area::GhostedArea, cluster_seq::ClusterSequence, jets::Vector{PseudoJet})

Calculates the area of each jet in a collection of jets based on the number of ghost particles they contain.

# Arguments
- `ghosted_area::GhostedArea`: The ghosted area configuration.
- `cluster_seq::ClusterSequence`: The clustering sequence of the event.
- `jets::Vector{PseudoJet}`: A collection of jets.

# Returns
A vector containing the calculated area of each jet.
"""
function ghosted_areas_calculation(ghosted_area::GhostedArea, cluster_seq::ClusterSequence, jets::Vector{PseudoJet})
    # area is equal to the number of ghosts in the jet divided by the density of ghosts
    return ghosts_in_jets(cluster_seq, jets) ./ ghosted_area.ghost_density
end