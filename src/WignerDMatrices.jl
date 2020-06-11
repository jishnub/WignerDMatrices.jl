module WignerDMatrices

using LinearAlgebra
using HalfIntegers

export WignerDMatrix
export WignerdMatrix
export WignerdMatrix!

# Dictionary to cache the eigenvectors of Jy
const JyEigenDict = Dict{HalfInt,Matrix{ComplexF64}}()

include("specialpoints.jl")
include("wignerdtypes.jl")
include("wignerdeval.jl")
include("linearalgebra.jl")

end # module
