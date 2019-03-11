# include("../src/Mozi.jl")

using Test
using Logging

using Mozi

const PATH=pwd()

include("./static_cantilever.jl")
include("./dynamic_cantilever.jl")
include("./buckling_cantilever.jl")

# include("./static_quad.jl")
# include("./dynamic_quad.jl")
# include("./buckling_quad.jl")

# include("./static_p_delta.jl")
# include("./dynamic_p_delta.jl")
# include("./buckling_p_delta.jl")
#
# include("./static_large_deform.jl")
# include("./dynamic_large_deform.jl")
# include("./buckling_large_deform.jl")

include("./dynamic_cantilever.jl")
cdiff=result_nodal_time_history(assembly,"DIFF",2,0,3)
newmark=result_nodal_time_history(assembly,"NEWMARK",2,0,3)
wilson=result_nodal_time_history(assembly,"WILSON",2,0,3)
modaldecomp=result_nodal_time_history(assembly,"MODALDECOMP",2,0,3)
using PyPlot
cla()
plot(t,cdiff)
plot(t,newmark)
# plot(t,wilson)
plot(t,modaldecomp)
display(gcf())
