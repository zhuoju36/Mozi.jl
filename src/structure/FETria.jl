# module FETria
#
# using LinearAlgebra
# using SparseArrays
# using Logging
#
# using HCubature
#
# using ...CoordinateSystem
# using ..FEMaterial
# using ..FENode
#
# export Tria

mutable struct Tria <: AbstractElement
    id::String
    hid::Int

    node1::Node
    node2::Node
    node3::Node
    material::Material
    t::Float64

    membrane::Bool
    bending::Bool

    center::Array{Float64}
    A::Float64
    T::SparseMatrixCSC{Float64}

    Kmᵉ::SparseMatrixCSC{Float64}
    Kbᵉ::SparseMatrixCSC{Float64}

    Kᵉ::SparseMatrixCSC{Float64}
    Mᵉ::SparseMatrixCSC{Float64}
end

function Tria(id,hid,node1,node2,node3,material,t;membrane=true,bending=true)
    #Initialize local CSys,could be optimized by using a MSE plane
    o=(node1.loc+node2.loc+node3.loc)/3
    pt1 = (node1.loc+node2.loc)/2
    pt2 = (node2.loc+node3.loc)/2
    csys = CSys(o, pt1, pt2)

    #area is considered as the average of trangles generated by splitting the quand with diagonals
    A=det([[1;1;1] node2.loc-node1.loc node3.loc-node1.loc]')/2

    E=material.E
    ν=material.ν

    #3D to local 2D
    T=zeros(18,18)
    for i in 1:6
        T[3i-2:3i,3i-2:3i]=csys.T
    end
    T=sparse(T)

    Kbᵉ=spzeros(1,1)
    Kmᵉ=spzeros(1,1)

    Kᵉ=spzeros(1,1)
    Mᵉ=spzeros(1,1)

    Tria(id,hid,node1,node2,node3,material,t,membrane,bending,o,A,T,Kbᵉ,Kmᵉ,Kᵉ,Mᵉ)
end

function integrateKm!(elm::Tria)
    E₀,ν₀=elm.material.E,elm.material.ν
    center=elm.center
    t=elm.t
    A=elm.A
    T=elm.T[1:3,1:3]
    x₁,y₁,z₁=T*(elm.node1.loc-center)
    x₂,y₂,z₂=T*(elm.node2.loc-center)
    x₃,y₃,z₃=T*(elm.node3.loc-center)

    a₁=x₂*y₃-x₃*y₂
    a₂=x₃*y₁-x₁*y₃
    a₃=x₁*y₂-x₂*y₁

    b₁=y₂-y₃
    b₂=y₃-y₁
    b₃=y₁-y₂

    c₁=-x₂+x₃
    c₂=-x₃+x₁
    c₃=-x₁+x₂

    # N₁=(a₁+b₁*x+c₁*y)/2/A
    # N₂=(a₂+b₂*x+c₂*y)/2/A
    # N₃=(a₃+b₃*x+c₃*y)/2/A
    #
    # N=[N₁ 0 N₂ 0 N₃ 0;
    #     0 N₁ 0 N₂ 0 N₃]

    B₁=[b₁ 0;
        0 c₁;
        c₁ b₁]
    B₂=[b₂ 0;
        0 c₂;
        c₂ b₂]
    B₃=[b₃ 0;
        0 c₃;
        c₃ b₃]

    B=[B₁ B₂ B₃]/2A
    D₀=E₀/(1-ν₀^2)
    D=D₀*[1 ν₀ 0;
          ν₀ 1 0;
          0 0 (1-ν₀)/2]
    elm.Kmᵉ=B'*D*B*t*A
end

#WIP
function integrateKb!(elm::Tria)
    E₀,ν₀=elm.material.E,elm.material.ν
    center=elm.center
    t=elm.t
    T=elm.T[1:3,1:3]
    x₁,y₁,z₁=T*(elm.node1.loc-center)
    x₂,y₂,z₂=T*(elm.node2.loc-center)
    x₃,y₃,z₃=T*(elm.node3.loc-center)
    # elm.Kbᵉ=sparse(hcubature(BtDB,[-1,1],[1,1])[1])
    elm.Kbᵉ=spzeros(9,9)
end

function integrateK!(elm::Tria)
    membrane,bending=elm.membrane,elm.bending
    Kᵉ=spzeros(18,18)
    if membrane
        I=1:6
        J=[1,2,7,8,13,14]
        L=sparse(I,J,1.,6,18)
        Km=integrateKm!(elm)
        Kᵉ+=L'*Km*L
    end
    if bending
        I=1:9
        J=[3,4,5,9,10,11,15,16,17]
        L=sparse(I,J,1.,9,18)
        Kb=integrateKb!(elm)
        Kᵉ+=L'*Kb*L
    end
    elm.Kᵉ=Kᵉ
end

#WIP
function integrateKσ(elm::Tria,σ::Vector{Float64})
    E₀,ν₀=elm.material.E,elm.material.ν
    center=elm.center
    t=elm.t
    T=elm.T[1:3,1:3]
    x₁,y₁,z₁=T*(elm.node1.loc-center)
    x₂,y₂,z₂=T*(elm.node2.loc-center)
    x₃,y₃,z₃=T*(elm.node3.loc-center)
    K=Matrix{Float64}(undef,8,8)
    J=Matrix{Float64}(undef,2,2)

    o=(node1.loc+node2.loc+node3.loc+node4.loc)/4
    n1=node1.loc+[u[1];u[2];0]
    n2=node2.loc+[u[3];u[4];0]
    n3=node3.loc+[u[5];u[6];0]
    pt1 = n1+n2
    pt2 = n2+n3
    csys = CSys(o, pt1, pt2)
    V=csys.T
    T̄ᵉ=zeros(8,8)
    T̄ᵉ[1:2,1:2]=V
    T̄ᵉ[3:4,3:4]=V
    T̄ᵉ[5:6,5:6]=V
    T̄ᵉ[7:8,7:8]=V
end

#WIP
function integrateKu(elm::Tria,u)
    E₀,ν₀=elm.material.E,elm.material.ν
    center=elm.center
    t=elm.t
    T=elm.T[1:3,1:3]
    x₁,y₁,z₁=T*(elm.node1.loc-center)
    x₂,y₂,z₂=T*(elm.node2.loc-center)
    x₃,y₃,z₃=T*(elm.node3.loc-center)
    x₄,y₄,z₄=T*(elm.node4.loc-center)
    K=Matrix{Float64}(undef,8,8)
    J=Matrix{Float64}(undef,2,2)

    o=(node1.loc+node2.loc+node3.loc+node4.loc)/4
    n1=node1.loc+[u[1];u[2];0]
    n2=node2.loc+[u[3];u[4];0]
    n3=node3.loc+[u[5];u[6];0]
    n4=node4.loc+[u[7];u[8];0]
    pt1 = n1+n2
    pt2 = n2+n3
    csys = CSys(o, pt1, pt2)
    V=csys.T
    T̄ᵉ=zeros(8,8)
    T̄ᵉ[1:2,1:2]=V
    T̄ᵉ[3:4,3:4]=V
    T̄ᵉ[5:6,5:6]=V
    T̄ᵉ[7:8,7:8]=V
end

function integrateM!(elm::Tria)
    membrane,bending=elm.membrane,elm.bending
    ρ=elm.material.ρ
    center=elm.center
    t=elm.t
    A=elm.A
    T=elm.T[1:3,1:3]
    x₁,y₁,z₁=T*(elm.node1.loc-center)
    x₂,y₂,z₂=T*(elm.node2.loc-center)
    x₃,y₃,z₃=T*(elm.node3.loc-center)
    # M=Matrix{Float64}(ρ*t*A/18,18,18)
    M=spzeros(18,18)
    elm.Mᵉ=M
end

function integrateP(elm::Tria,elm_force)
    N=elm.N
    B=elm.B
    D=elm.D
    J=elm.J

    f,s,σ₀,ϵ₀=elmforce.f,elmforce.s,elmforce.σ₀,elmforce.ϵ₀
    a=[-1 -1]
    b=[1 1]

    function f1(x)
        N'*f*det(J)
    end

    function f2(x)
        N'*s*det(J)
    end

    function f3(x)
        x->B'*σ₀*det(J)
    end

    function f4(x)
        x->B'*D*ϵ₀*det(J)
    end

    Pᵉ_f=hTriarature(f1,a,b)
    Pᵉ_s=hTriarature(f2,[-1],[1])#ξ=-1
    Pᵉ_σ₀=hTriarature(f3,a,b)
    Pᵉ_ϵ₀=hTriarature(f4,a,b)

    Pᵉ_f,Pᵉ_s,Pᵉ_σ₀,Pᵉ_ϵ₀=calc_force(elm,elm_force) #Pᵉ=Pᵉf+Pᵉs+Pᵉσ₀+Pᵉϵ₀
    Pᵉ=Pᵉ_f+Pᵉ_s+Pᵉ_σ₀+Pᵉ_ϵ₀
end

# end
