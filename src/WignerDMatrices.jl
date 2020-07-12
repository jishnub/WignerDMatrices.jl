module WignerDMatrices

using LinearAlgebra
using HalfIntegers
using HalfIntegerArrays
import HalfIntegerArrays: AbstractHalfIntegerMatrix

import Base: @propagate_inbounds

export WignerDMatrix
export WignerdMatrix

# Dictionary to cache the eigenvectors of Jy
const JyEigenDict = Dict{HalfInt,Matrix{ComplexF64}}()

include("specialpoints.jl")
include("wignerdtypes.jl")
include("wignerdeval.jl")
include("linearalgebra.jl")

end # module
