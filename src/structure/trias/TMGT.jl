function K_TMGT(elm::Tria)::Matrix{Float64}
    return K_GT9(elm)+K_TMT(elm)
end

function M_TMGT(elm::Tria)::Matrix{Float64}
    return M_GT9(elm)
end

# function C_TMGT(elm::Tria)::Matrix{Float64}
# end

function P_TMGT(elm::Tria,p::Vector{Float64})::Vector{Float64}
end

function W_TMGT(elm::Tria,u::Vector{Float64})::Vector{Float64}
end

function f_TMGT(elm::Tria,u::Vector{Float64})::Vector{Float64}
end
