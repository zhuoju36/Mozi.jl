export Beam

mutable struct Beam <: AbstractElement
    id::String
    hid::Int
    node1::Node
    node2::Node
    material::Material
    section::BeamSection

    release::Vector{Bool}

    elm_type::String
    mass_type::String

    center::Vector{Float64}
    l::Float64
    T::Matrix{Float64} #transform_matrix
end

function Beam(id,hid,node1,node2,material,section;elm_type="eular_shear",mass_type="concentrate")
    tol=1e-6
    o=node1.loc
    pt1=node2.loc
    pt2=node1.loc
    if abs(pt2[1]-pt1[1])<tol && abs(pt2[2]-pt1[2])<tol
        pt2=pt2.+[1,0,0]
    else
        pt2=pt2.+[0,0,1]
    end
    csys=CSys(o,pt1,pt2)
    T=zeros(12,12)
    T[1:3,1:3]=T[4:6,4:6]=T[7:9,7:9]=T[10:12,10:12]=csys.T
    l=norm(node1.loc-node2.loc)
    release=Bool.(zeros(12))
    Beam(string(id),hid,node1,node2,material,section,release,elm_type,mass_type,(pt1+pt2)/2,l,T)
end

for (root,dirs,files) in walkdir(joinpath(@__DIR__,"beams"))
    for file in files
        if file[end-2:end]==".jl"
            include(joinpath(@__DIR__,"beams",file))
        end
    end
end

function integrateK(beam::Beam)::Matrix{Float64}
    if beam.elm_type=="eular_shear"
        return K_eular_shear(beam::Beam)
    end
end

function integrateKσ(beam::Beam,d::Vector{Float64})::Matrix{Float64}
    if beam.elm_type=="eular_shear"
        return K2_eular_shear(beam::Beam,d)
    end
end

function static_condensation(K,P,rDOF::Vector{Int})
    if isempty(rDOF)
        return K,P
    end
    i=[!(x in rDOF) for x in 1:12]
    j=[(x in rDOF) for x in 1:12]
    Kᵢᵢ=K[i,i]
    Kᵢⱼ=K[i,j]
    Kⱼᵢ=K[j,i]
    Kⱼⱼ=K[j,j]
    Kⱼⱼ⁻¹=inv(Kⱼⱼ)
    Pⱼ=P[j]
    Pᵢ=P[i]
    Kᶜ=zero(K)
    Pᶜ=zero(P)
    Kᶜ[i,i]=Kᵢᵢ-Kᵢⱼ*Kⱼⱼ⁻¹*Kⱼᵢ
    Pᶜ[i]=Pᵢ-Kᵢⱼ*Kⱼⱼ⁻¹*Pⱼ
    return Kᶜ,Pᶜ
end

function integrateM(beam::Beam)::Matrix{Float64}
    E,ν=beam.material.E,beam.material.ν
    A,I₂,I₃,J,l=beam.section.A,beam.section.I₂,beam.section.I₃,beam.section.J,beam.l
    ρ=beam.material.ρ
    if beam.mass_type=="concentrate"
        return Matrix(I,12,12)*12*ρ*A*l/2
    elseif beam.mass_type=="coordinate"
        return Matrix(I,12,12)*12*ρ*A*l/2
    end
end

function integrateP(beam::Beam,beam_force)::Vector{Float64}
    if beam.elm_type=="eular_shear"
        return P_eular_shear(beam::Beam,beam_force.f,beam_force.ϵ₀)
    end
end

# end
