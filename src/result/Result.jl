module Result
include("../util/mmp.jl")

include("./nodal_result.jl")
include("./beam_result.jl")
include("./modal_result.jl")

using ..LoadCase

export result_nodal_reaction,result_nodal_displacement,result_nodal_time_history

export result_beam_force,result_beam_displacement

export result_modal_period,result_eigen_value,result_eigen_vector

end
