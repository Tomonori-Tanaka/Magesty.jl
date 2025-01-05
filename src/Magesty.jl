module Magesty

include("common/SortedContainer.jl")

include("types/AtomCell.jl")
include("types/AtomicIndices.jl")
include("types/UnitaryMatrixCl.jl")

include("System.jl")
include("Symmetry.jl")

using .Systems
using .Symmetries

struct SpinCluster
    system::System
    symmetry::Symmetry
end

end
