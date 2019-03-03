module FEBeam

using LinearAlgebra
using SparseArrays
using HCubature

using ...CoordinateSystem
import ..FENode

export Beam

mutable struct Beam
    id::String
    hid::Int
    node1
    node2
    material
    section

    release::Vector{Bool}

    l::Float64
    T::Matrix{Float64} #transform_matrix
    K̄ᵉ::SparseMatrixCSC{Float64}
    M̄ᵉ::SparseMatrixCSC{Float64}
    P̄ᵉ::SparseMatrixCSC{Float64}
    Kᵉ::SparseMatrixCSC{Float64}
    Mᵉ::SparseMatrixCSC{Float64}
    Pᵉ::SparseMatrixCSC{Float64}
end

function Beam(id,hid,node1,node2,material,section)
    tol=1e-6
    o=node1.loc
    pt1=node2.loc
    pt2=node1.loc
    if abs(pt2[1]-pt1[1])<tol && abs(pt2[2]-pt1[2])<tol
        pt2=pt2+[1,0,0]
    else
        pt2=pt2+[0,0,1]
    end
    csys=CSys(o,pt1,pt2)
    T=zeros(12,12)
    T[1:3,1:3]=T[4:6,4:6]=T[7:9,7:9]=T[10:12,10:12]=csys.T
    l=norm(node1.loc-node2.loc)
    release=Bool.(zeros(12))
    K̄=spzeros(12,12)
    M̄=spzeros(12,12)
    P̄=spzeros(12,1)
    K=spzeros(12,12)
    M=spzeros(12,12)
    P=spzeros(12,1)
    Beam(string(id),hid,node1,node2,material,section,release,l,T,K̄,M̄,P̄,K,M,P)
end

function integrateK!(beam::Beam)::SparseMatrixCSC
    E,ν=beam.material.E,beam.material.ν
    A,I₂,I₃,J,l=beam.section.A,beam.section.I₂,beam.section.I₃,beam.section.J,beam.l
    As₂,As₃=beam.section.As₂,beam.section.As₃
    G=E/2/(1+ν)
    ϕ₂,ϕ₃=12E*I₃/(G*As₂*l^2),12E*I₂/(G*As₃*l^2)
    K=zeros(12,12)
    K[1,1]=E*A/l
    K[2,2]=12E*I₃/l^3/(1+ϕ₂)
    K[3,3]=12E*I₂/l^3/(1+ϕ₃)
    K[4,4]=G*J/l
    K[5,5]=(4+ϕ₃)*E*I₂/l/(1+ϕ₃)
    K[6,6]=(4+ϕ₂)*E*I₃/l/(1+ϕ₂)
    K[7,7]=E*A/l
    K[8,8]=12E*I₃/l^3/(1+ϕ₂)
    K[9,9]=12E*I₂/l^3/(1+ϕ₃)
    K[10,10]=G*J/l
    K[11,11]=(4+ϕ₃)*E*I₂/l/(1+ϕ₃)
    K[12,12]=(4+ϕ₂)*E*I₃/l/(1+ϕ₂)

    K[3,5]=K[5,3]=-6E*I₂/l^2/(1+ϕ₃)
    K[6,8]=K[8,6]=-6E*I₃/l^2/(1+ϕ₂)
    K[9,11]=K[11,9]=6E*I₂/l^2/(1+ϕ₃)

    K[2,6]=K[6,2]=6E*I₃/l^2/(1+ϕ₂)
    K[5,9]=K[9,5]=6E*I₂/l^2/(1+ϕ₃)
    K[8,12]=K[12,8]=-6E*I₃/l^2/(1+ϕ₂)

    K[7,1]=K[1,7]=-E*A/l
    K[8,2]=K[2,8]=-12E*I₃/l^3/(1+ϕ₂)
    K[9,3]=K[3,9]=-12E*I₂/l^3/(1+ϕ₃)
    K[10,4]=K[4,10]=-G*J/l
    K[11,5]=K[5,11]=(2-ϕ₃)*E*I₂/l/(1+ϕ₃)
    K[12,6]=K[6,12]=(2-ϕ₂)*E*I₃/l/(1+ϕ₂)

    K[3,11]=K[11,3]=-6E*I₂/l^2/(1+ϕ₃)

    K[2,12]=K[12,2]=6E*I₃/l^2/(1+ϕ₂)

    beam.K̄ᵉ=K
    rDOF=findall(x->x==true,beam.release)
    if length(rDOF)!=0
        beam.Kᵉ=sparse(static_condensation(K,zeros(12),rDOF)[1])
    else
        beam.Kᵉ=K
    end
    return sparse(beam.Kᵉ)
end

function integrateKσ(beam::Beam,σ)::SparseMatrixCSC
    E,ν=beam.material.E,beam.material.ν
    A,I₂,I₃,J,l=beam.section.A,beam.section.I₂,beam.section.I₃,beam.section.J,beam.l
    As₂,As₃=beam.section.As₂,beam.section.As₃
    G=E/2/(1+ν)
    T=σ*A

    ϕ₂,ϕ₃=12E*I₃/(G*As₂*l^2),12E*I₂/(G*As₃*l^2)
    K=zeros(12,12)

    K[2,2]=(6/5+2ϕ₂+ϕ₂^2)/(1+ϕ₂)^2
    K[3,3]=(6/5+2ϕ₃+ϕ₃^2)/(1+ϕ₃)^2
    K[4,4]=J/A
    K[5,5]=(2*l^2/15+l^2*ϕ₃/6+l^2*ϕ₃^2/12)/(1+ϕ₃)^2
    K[6,6]=(2*l^2/15+l^2*ϕ₂/6+l^2*ϕ₂^2/12)/(1+ϕ₂)^2
    K[8,8]=(6/5+2ϕ₂+ϕ₂^2)/(1+ϕ₂)^2
    K[9,9]=(6/5+2ϕ₃+ϕ₃^2)/(1+ϕ₃)^2
    K[10,10]=J/A
    K[11,11]=(2*l^2/15+l^2*ϕ₃/6+l^2*ϕ₃^2/12)/(1+ϕ₃)^2
    K[12,12]=(2*l^2/15+l^2*ϕ₂/6+l^2*ϕ₂^2/12)/(1+ϕ₂)^2

    K[3,5]=K[5,3]=-(l/10)/(1+ϕ₃)^2
    K[6,8]=K[8,6]=-(l/10)/(1+ϕ₂)^2
    K[9,11]=K[11,9]=(l/10)/(1+ϕ₃)^2

    K[2,6]=K[6,2]=(l/10)/(1+ϕ₂)^2
    K[5,9]=K[9,5]=(l/10)/(1+ϕ₃)^2

    K[2,8]=K[8,2]=-(6/5+2ϕ₂+ϕ₂^2)/(1+ϕ₂)^2
    K[3,9]=K[9,3]=-(6/5+2ϕ₃+ϕ₃^2)/(1+ϕ₃)^2
    K[4,10]=K[10,4]=-J/A
    K[5,11]=K[11,5]=-(l^2/30+l^2*ϕ₃/6+l^2*ϕ₃^2/12)/(1+ϕ₃)^2
    K[6,12]=K[12,6]=-(l^2/30+l^2*ϕ₂/6+l^2*ϕ₂^2/12)/(1+ϕ₂)^2

    K[3,11]=K[11,3]=-(l/10)/(1+ϕ₃)^2

    K[2,12]=K[12,2]=(l/10)/(1+ϕ₂)^2

    rDOF=findall(x->x==true,beam.release)
    if length(rDOF)!=0
        Kᵉ=sparse(static_condensation(T/l*K,zeros(12),rDOF)[1])
    else
        Kᵉ=T/l*K
    end
    return sparse(Kᵉ)
end

# function integrateK!(beam,p_delta=false,σ₀=0)::SparseMatrixCSC
#     Kᵉ=integrateKi(beam)
#     if p_delta
#         Kᵉ+=integrateKσ(beam,σ₀)
#     end
#     beam.K̄ᵉ=Kᵉ
#     rDOF=findall(x->x==true,beam.release)
#     if length(rDOF)!=0
#         beam.Kᵉ=sparse(static_condensation(Array(beam.K̄ᵉ),zeros(12),rDOF)[1])
#     else
#         beam.Kᵉ=beam.K̄ᵉ
#     end
#     return Kᵉ
# end

function static_condensation(K,P,rDOF)
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
    # Kᶜ[i,j].=0
    # Kᶜ[i,j].=0
    # Kᶜ[j,j].=0

    Pᶜ[i]=Pᵢ-Kᵢⱼ*Kⱼⱼ⁻¹*Pⱼ
    # Pᶜ[j].=0
    return Kᶜ,Pᶜ
end


function integrateM!(beam::Beam)
    E,ν=beam.material.E,beam.material.ν
    A,I₂,I₃,J,l=beam.section.A,beam.section.I₂,beam.section.I₃,beam.section.J,beam.l
    ρ=beam.material.ρ
    beam.M̄ᵉ=sparse(Matrix(I,12,12)*12*ρ*A*l/2)
    rDOF=findall(x->x==true,beam.release)
    if length(rDOF)!=0
        beam.Mᵉ=sparse(Matrix(I,12,12)*12*ρ*A*l/2)
    else
        beam.Mᵉ=beam.M̄ᵉ
    end
    return beam.Mᵉ
end

function calc_force(beam,beamforce)
    f₁,f₂=beamforce.f[1:6],beamforce.f[7:12]
    l=beam.l
    fi,v2i,v3i,ti,m2i,m3i,fj,v2j,v3j,tj,m2j,m3j=beamforce.f
    Nᵀf(x)=[
          0;
    (1 - 3*x^2 + 2*x^3)*(v2i + x*(-v2i + v2j));
    (1 - 3*x^2 + 2*x^3)*(v3i + x*(-v3i + v3j));
          0;
    l*(x - 2*x^2 + x^3)*(v3i + x*(-v3i + v3j));
    l*(x - 2*x^2 + x^3)*(v2i + x*(-v2i + v2j));
          0;
        (3*x^2 - 2*x^3)*(v2i + x*(-v2i + v2j));
        (3*x^2 - 2*x^3)*(v3i + x*(-v3i + v3j));
         0;
         l*(-x^2 + x^3)*(v3i + x*(-v3i + v3j));
         l*(-x^2 + x^3)*(v2i + x*(-v2i + v2j));
   ]*l

   Nᵀf2(x)=[
         -0.5*(-1 + x)*(fi + (1+x)/2*(-fi + fj));0;0;
         -0.5*(-1 + x)*(ti + (1+x)/2*(-ti + tj));0;0;
           0.5*(1 + x)*(fi + (1+x)/2*(-fi + fj));0;0;
           0.5*(1 + x)*(ti + (1+x)/2*(-ti + tj));0;0
  ]*l/2

    Pᵉ_f=hquadrature(Nᵀf,0,1)[1]+hquadrature(Nᵀf2,-1,1)[1]#Pᵉf=∫NᵀfdV #体积力

    Pᵉ_s=beamforce.s# Pᵉs=∫NᵀTdS #边界力

    f₁,f₂=beamforce.σ₀[1],beamforce.σ₀[2]
    Bᵀσ₀(x)=[-0.5*(f₁+0.5*(x+1)*(-f₁+f₂));
              0.5*(f₁+0.5*(x+1)*(-f₁+f₂))]*beam.l/2
    a,b=hquadrature(Bᵀσ₀,-1,1)[1]# Pᵉσ₀=-∫Bᵀσ₀dV #初应力
    Pᵉ_σ₀=[a,0,0,0,0,0,b,0,0,0,0,0]

    f₁,f₂=beamforce.ϵ₀[1],beamforce.ϵ₀[2]
    BᵀDϵ₀(x)=[-0.5*(f₁+0.5*(x+1)*(-f₁+f₂));
               0.5*(f₁+0.5*(x+1)*(-f₁+f₂))]*beam.l/2
    a,b=hquadrature(BᵀDϵ₀,-1,1)[1]*beam.material.E*beam.section.A# Pᵉϵ₀=∫BᵀDϵ₀dV #初应变
    Pᵉ_ϵ₀=[a,0,0,0,0,0,b,0,0,0,0,0]
    return Pᵉ_f,Pᵉ_s,Pᵉ_σ₀,Pᵉ_ϵ₀
end

function integrateP!(beam::Beam,beam_force)
    Pᵉ_f,Pᵉ_s,Pᵉ_σ₀,Pᵉ_ϵ₀=calc_force(beam,beam_force) #Pᵉ=Pᵉf+Pᵉs+Pᵉσ₀+Pᵉϵ₀
    beam.P̄ᵉ=Pᵉ_f+Pᵉ_s+Pᵉ_σ₀+Pᵉ_ϵ₀
    rDOF=findall(x->x==true,beam.release)
    if length(rDOF)!=0
        beam.Pᵉ=sparse(static_condensation(Array(beam.K̄ᵉ),Array(beam.P̄ᵉ),rDOF)[2])
    else
        beam.Pᵉ=beam.P̄ᵉ
    end
    return beam.Pᵉ
end

end
