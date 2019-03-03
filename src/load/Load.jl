module Load

using LinearAlgebra
using SparseArrays
using HCubature

export NodalForce,NodalDisp,
BeamForce,QuadForce

struct NodalForce
    id::String
    val::Vector{Float64}
    csys::String #This is a string to indicate whether the value is given in global or local csys.
    NodalForce(id,force,csys="global")=new(id,force,csys)
end

struct NodalDisp
    id::String
    val::Vector{Float64}
    csys::String
    NodalDisp(id,disp,csys="global")=new(id,disp,csys)
end

mutable struct BeamForce
    id::String
    f::Array{Float64}
    s::Array{Float64}
    σ₀::Array{Float64}
    ϵ₀::Array{Float64}
    BeamForce(id)=new(id,zeros(12,1),zeros(12,1),zeros(2,1),zeros(2,1))
end

mutable struct QuadForce
    id::String
    f::Array{Float64}
    s::Array{Float64}
    σ₀::Array{Float64}
    ϵ₀::Array{Float64}
    QuadForce(id)=new(id,zeros(24,1),zeros(24,1),zeros(4,1),zeros(4,1))
end

end
