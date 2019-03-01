# include("../src/Mozi.jl")

using Test
using Logging

using Mozi

const PATH=pwd()

include("./static_cantilever.jl")
include("./dynamic_cantilever.jl")
include("./buckling_cantilever.jl")

include("./static_quad.jl")
include("./dynamic_quad.jl")
include("./buckling_quad.jl")

include("./static_p_delta.jl")
include("./dynamic_p_delta.jl")
include("./buckling_p_delta.jl")

include("./static_large_deform.jl")
include("./dynamic_large_deform.jl")
include("./buckling_large_deform.jl")
