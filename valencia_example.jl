println("Starting Valencia test script...")

using LorentzVectorHEP
println("Loaded LorentzVectorHEP")

using JetReconstruction
println("Loaded JetReconstruction")

# Create 3 example particles in a cone
particles = [
    LorentzVector(50.0, 10.0, 0.0, sqrt(50.0^2 - 10.0^2)),
    LorentzVector(40.0, -10.0, 0.0, sqrt(40.0^2 - 10.0^2)),
    LorentzVector(20.0, 0.0, 50.0, sqrt(20.0^2 + 50.0^2))
]
println("Defined 3 particles")

# Use the Valencia algorithm with β = 1.0 and R = 1.0
jets = jet_reconstruct(particles; algorithm=JetAlgorithm.Valencia, p=1.0, R=1.0)

# Print result
println("Reconstructed jets:")
for jet in inclusive_jets(jets)
    println(jet)
end