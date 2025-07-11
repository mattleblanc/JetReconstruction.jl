using ArgParse
using Logging

using LorentzVectorHEP
using JetReconstruction

# Parsing for algorithm and strategy enums
include(joinpath(@__DIR__, "..", "parse-options.jl"))

function parse_command_line(args)
    s = ArgParseSettings(autofix_names = true)
    @add_arg_table! s begin
        "--maxevents", "-n"
        help = "Maximum number of events to read. -1 to read all events from the  file."
        arg_type = Int
        default = -1

        "--skip", "-s"
        help = "Number of events to skip at beginning of the file."
        arg_type = Int
        default = 0

        "--algorithm", "-A"
        help = """Algorithm to use for jet reconstruction: $(join(JetReconstruction.AllJetRecoAlgorithms, ", "))"""
        arg_type = JetAlgorithm.Algorithm

        "--power", "-p"
        help = """Power value for jet reconstruction"""
        arg_type = Float64

        "--strategy", "-S"
        help = """Strategy for the algorithm, valid values: $(join(JetReconstruction.AllJetRecoStrategies, ", "))"""
        arg_type = RecoStrategy.Strategy
        default = RecoStrategy.Best

        "--resolution"
        help = "Square root of number of ghosts used"
        arg_type = Int
        default = 1000

        "--data-file"
        help = "HepMC3 event file in HepMC3 to read."
        required = true
    end
    return parse_args(args, s; as_symbols = true)
end

# Cluster an event, run after adding ghosts
function cluster_event(event::Vector{PseudoJet}, args::Dict{Symbol, Any})

    distance = args[:distance]
    algorithm = args[:algorithm]
    p = args[:power]
    strategy = args[:strategy]

    new_event = PseudoJet[]
    for (i, pseudo_jet) in enumerate(event)
    # Reconstruct PseudoJet with cluster_hist_index for tracking
    new_pseudo_jet = PseudoJet(JetReconstruction.px(pseudo_jet),
                    JetReconstruction.py(pseudo_jet),
                    JetReconstruction.pz(pseudo_jet),
                    JetReconstruction.energy(pseudo_jet);
                    cluster_hist_index = i)
    push!(new_event, new_pseudo_jet)
    end

    # Create the clustering of the given vector of PseudoJets 
    cluster_seq = jet_reconstruct(new_event,
                    R = distance, p = p, algorithm = algorithm,
                    strategy = strategy)

    clustered_jets = inclusive_jets(cluster_seq, ptmin = 0.0)
    return (cluster_seq, clustered_jets)
end

function main()
    args = parse_command_line(ARGS)
    logger = ConsoleLogger(stdout, Logging.Info)
    global_logger(logger)

    # Only PseudoJet is supported for Area Calculation
    @assert JetReconstruction.is_pp(args[:algorithm]) "Area Calculation only supports pp algorithms and PseudoJet"
    jet_type = PseudoJet

    events = Vector{PseudoJet}[]

    args[:data_file] = normpath(joinpath(@__DIR__, args[:data_file]))

    # Reading event files
    events = read_final_state_particles(args[:data_file],
                                          maxevents = args[:maxevents],
                                          skipevents = args[:skip],
                                          T = jet_type)

    # Set up GhostedArea structure
    resolution = args[:resolution]
    ghosted_area = GhostedArea(resolution)

    # Ensure algorithm and power are consistent
    if isnothing(args[:algorithm]) && isnothing(args[:power])
        @warn "Neither algorithm nor power specified, defaulting to AntiKt"
        args[:algorithm] = JetAlgorithm.AntiKt
    end

    # all_jets: go through events and add jets to a single array
    all_jets = PseudoJet[]

    # Fill all_jets
    for event in events
        for jet in event
            push!(all_jets, jet)
        end
    end

    # Add ghosts to all_jets
    add_ghosts!(ghosted_area, all_jets)

    # Create clustering of PseudoJets
    clustered_event = cluster_event(all_jets, args)
    cluster_seq = cluster_event[1]
    clustered_jets = cluster_event[2]

    # Calculate areas and print them
    areas_vector = ghosted_areas_calculation(ghosted_area, cluster_seq, clustered_jets)
    println(areas_vector)
end

main()