using ArgParse
using Logging
using Plots
using CairoMakie

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

        "--distance", "-R"
        help = "Distance parameter for jet merging"
        arg_type = Float64
        default = 0.4

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

    # Create the clustering of the given vector of PseudoJets 
    cluster_seq = jet_reconstruct(event,
                    R = args[:distance], p = args[:power], algorithm = args[:algorithm],
                    strategy = args[:strategy])

    # Get the clustered jets from the cluster sequence
    clustered_jets = inclusive_jets(cluster_seq, ptmin = 10.0, T = PseudoJet)
    return (cluster_seq, clustered_jets)
end

function cluster_events(events::Vector{Vector{PseudoJet}}, args::Dict{Symbol, Any})
    # Cluster each event and return a vector of tuples (cluster_seq, clustered_jets)
    clustered_events = map(event -> cluster_event(event, args), events)
    return collect(zip(clustered_events))
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

    event_areas_vector = Vector{Float64}[]
    # Loop through each event, add ghosts, cluster, and calculate areas
    for event in events
        add_ghosts!(ghosted_area, event)
        cluster_seq, clustered_jets = cluster_event(event, args)

        # Plot the first event
        if event == events[1]
            plt = jetsplot(event, cluster_seq; Module = CairoMakie)
            save("clustering.png", plt)
        end
        
        areas_vector = ghosted_areas_calculation(ghosted_area, cluster_seq, clustered_jets)
        push!(event_areas_vector, areas_vector)
    end
    
    # Print the areas for each event
    println(event_areas_vector)
    
    # Make a histogram of all the ghosted areas in every event
    h1 = histogram(vcat(event_areas_vector...), bins = 20, title = "Ghosted Areas", xlabel = "Area", ylabel = "Count", fmt = :png)
    savefig(h1, "ghosted_areas.png")
end

main()