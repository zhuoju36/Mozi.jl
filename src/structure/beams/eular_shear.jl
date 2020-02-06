function K_eular_shear(beam::Beam)::Matrix{Float64}
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

    return K
end

function K2_eular_shear(beam::Beam,u::Vector{Float64})::Matrix{Float64}
    E,ν=beam.material.E,beam.material.ν
    A,I₂,I₃,J,l=beam.section.A,beam.section.I₂,beam.section.I₃,beam.section.J,beam.l
    As₂,As₃=beam.section.As₂,beam.section.As₃
    G=E/2/(1+ν)
    T=E*A*(u[7]-u[1])/l
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

    K*=T/l
    return K
end

function P_eular_shear(beam::Beam,f::Vector{Float64},ϵ₀::Vector{Float64},)::Vector{Float64}
    l=beam.l
    fi,v2i,v3i,ti,m2i,m3i,fj,v2j,v3j,tj,m2j,m3j=f
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

    Pᵉᶠ=hquadrature(Nᵀf,0,1)[1]+hquadrature(Nᵀf2,-1,1)[1]#Pᵉf=∫NᵀfdV #体积力

    f₁,f₂=ϵ₀[1],ϵ₀[2]
    BᵀDϵ₀(x)=[-0.5*(f₁+0.5*(x+1)*(-f₁+f₂));
               0.5*(f₁+0.5*(x+1)*(-f₁+f₂))]*beam.l/2
    a,b=hquadrature(BᵀDϵ₀,-1,1)[1]*beam.material.E*beam.section.A# Pᵉϵ₀=∫BᵀDϵ₀dV #初应变
    Pᵉᵋ=[a,0,0,0,0,0,b,0,0,0,0,0]
    return reshape(Pᵉᶠ+Pᵉᵋ,12)
end

function f_eular_shear(elm::Beam,uᵉ::Vector{Float64})::Vector{Float64}
    rDOF=findall(x->x==true,elm.release)
    Kᵉ=K_eular_shear(elm)
    K̃ᵉ,P̃ᵉ=static_condensation(K̃ᵉ,zeros(12),rDOF)
    return K̃ᵉ*uᵉ
end

function W_eular_shear(elm::Beam,uᵉ::Vector{Float64})::Float64
    rDOF=findall(x->x==true,elm.release)
    Kᵉ=K_eular_shear(elm)
    K̃ᵉ,P̃ᵉ=static_condensation(K̃ᵉ,zeros(12),rDOF)
    return uᵉ'*(K̃ᵉ*uᵉ)
end
