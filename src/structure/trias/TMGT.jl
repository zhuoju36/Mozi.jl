function K_TMGT(elm::Tria)::Matrix{Float64}
    return K_GT9(elm)+K_TMT(elm)
end

function M_TMGT(elm::Tria)::Matrix{Float64}
    ρ=elm.material.ρ
    A=elm.A
    t=elm.t
    W=ρ*A*t
    M=zeros(18,18)
    for i in [1,2,3,7,8,9,13,14,15]
        M[i,i]=W/9
    end
    return M
end

function P_TMGT(elm::Tria,p::Vector{Float64})::Vector{Float64}
end

function W_TMGT(elm::Tria,u::Vector{Float64})::Vector{Float64}
end

function f_TMGT(elm::Tria,u::Vector{Float64})::Vector{Float64}
end
