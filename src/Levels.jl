module Levels

include("optics_utils.jl")

include("angular_momentum.jl")

include("states.jl")
include("basis.jl")

include("species.jl")
include("zeeman.jl")
include("hyperfine.jl")

include("multipole.jl")
include("rates.jl")
include("laser.jl")
include("polarisability.jl")

# Constructs the sr88/ca43 datasets at load time, so it must come last: the
# ReducedDipole/ImplicitPolarisability resolution the species constructors run
# (species.jl) calls multipole_rank (rates.jl) and the
# channel_scalar/tensor_polarisability helpers (polarisability.jl) at that
# point.
include("species_data.jl")

include("periodic_dynamics/PeriodicDynamics.jl")

end # module
