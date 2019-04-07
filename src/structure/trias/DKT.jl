#Reference:
#Batoz J.L., A Study of Three-Node Triangular Plate Bending Elements
function K_DKT(elm::Tria)::Matrix{Float64}
    E₀,ν₀=elm.material.E,elm.material.ν
    center=elm.center
    t=elm.t
    T=elm.T[1:3,1:3]
    x₁,y₁,z₁=T*(elm.node1.loc.-center)
    x₂,y₂,z₂=T*(elm.node2.loc.-center)
    x₃,y₃,z₃=T*(elm.node3.loc.-center)
    K=Matrix{Float64}(undef,12,12)
    J=Matrix{Float64}(undef,2,2)
    j=Matrix{Float64}(undef,2,2)
    function BtDB(ξ,η)
        x=[x₁,x₂,x₃,x₄]
        y=[y₁,y₂,y₃,y₄]

        xᵢⱼ=[x[i]-x[j] for (i,j) in zip([1,2,3,4],[2,3,4,1])]
        yᵢⱼ=[y[i]-y[j] for (i,j) in zip([1,2,3,4],[2,3,4,1])]
        lᵢⱼ²=xᵢⱼ.^2 .+ yᵢⱼ.^2

        a₅,a₆,a₇,a₈=-xᵢⱼ ./ lᵢⱼ²
        b₅,b₆,b₇,b₈=0.75*xᵢⱼ.*yᵢⱼ ./ lᵢⱼ²
        c₅,c₆,c₇,c₈=(0.25*xᵢⱼ.^2 .- 0.5*yᵢⱼ.^2) ./ lᵢⱼ²
        d₅,d₆,d₇,d₈=-yᵢⱼ ./ lᵢⱼ²
        e₅,e₆,e₇,e₈=(-0.25*xᵢⱼ.^2 .+ 0.5*yᵢⱼ.^2) ./ lᵢⱼ²

        dHˣdξ=[1.5*(4*η*a₅ - 4*ξ*a₆ + 4*(1 - η - ξ)*a₆), -4*η*b₅, 4*η*c₅ + 4*ξ*c₆ - 4*(1 - η - ξ)*c₆ - 2*(1 - η - ξ) - 2*(0.5 - η - ξ), 1.5*(4*η*a₄ + 4*ξ*a₆ - 4*(1 - η - ξ)*a₆), -4*ξ*b₆ + 4*(1 - η - ξ)*b₆, 2*ξ - 4*η*c₄ + 4*ξ*c₆ - 4*(1 - η - ξ)*c₆ - 1 + 2*ξ, 1.5*(-4*η*a₄ - 4*η*a₅), 4*η*b₄, -4*η*c₄ + 4*η*c₅ - 2*(1 - η - ξ) - 2*(0.5 - η - ξ)]
        dHˣdη=[1.5*(4*η*a₅ - 4*ξ*a₆ - 4*(1 - η - ξ)*a₅), -4*η*b₅ + 4*(1 - η - ξ)*b₅, 4*η*c₅ + 4*ξ*c₆ - 4*(1 - η - ξ)*c₅ - 2*(1 - η - ξ) - 2*(0.5 - η - ξ), 1.5*(4*ξ*a₄ + 4*ξ*a₆), -4*ξ*b₆, -4*ξ*c₄ + 4*ξ*c₆, 1.5*(-4*η*a₅ - 4*ξ*a₄ + 4*(1 - η - ξ)*a₅), 4*ξ*b₄, 4*η*c₅ - 4*ξ*c₄ - 4*(1 - η - ξ)*c₅ - 2*(1 - η - ξ) - 2*(0.5 - η - ξ)]
        dHʸdξ=[1.5*(4*η*d₅ - 4*ξ*d₆ + 4*(1 - η - ξ)*d₆), -4*η*e₅ - 4*ξ*e₆ + 4*(1 - η - ξ)*e₆ + 2*(1 - η - ξ) + 2*(0.5 - η - ξ), 4*η*b₅, 1.5*(4*η*d₄ + 4*ξ*d₆ - 4*(1 - η - ξ)*d₆), -2*ξ + 4*η*e₄ - 4*ξ*e₆ + 4*(1 - η - ξ)*e₆ - (-1 + 2*ξ), 4*η*b₅, 1.5*(-4*η*d₄ - 4*η*d₅), 4*η*e₄ - 4*η*e₅, 4*η*b₅]
        dHʸdη=[1.5*(4*η*d₅ - 4*ξ*d₆ - 4*(1 - η - ξ)*d₅), -4*η*e₅ - 4*ξ*e₆ + 4*(1 - η - ξ)*e₅ + 2*(1 - η - ξ) + 2*(0.5 - η - ξ), -(-4*η*b₅ + 4*(1 - η - ξ)*b₅), 1.5*(4*ξ*d₄ + 4*ξ*d₆), 4*ξ*e₄ - 4*ξ*e₆, -(-4*η*b₅ + 4*(1 - η - ξ)*b₅), 1.5*(-4*η*d₅ - 4*ξ*d₄ + 4*(1 - η - ξ)*d₅), -2*η - 4*η*e₅ + 4*ξ*e₄ + 4*(1 - η - ξ)*e₅ - (-1 + 2*η), -(-4*η*b₅ + 4*(1 - η - ξ)*b₅)]

        x₁₂,x₂₃,x₃₁=xᵢⱼ
        y₁₂,y₂₃,y₃₁=xᵢⱼ
        A=(x₃₁*y₁₂-x₁₂*y₃₁)
        B=1/A*[y₃₁*dHˣdξ+y₁₂*dHˣdη;
                    -x₃₁*dHʸdξ-x₁₂*dHʸdη;
                    -x₃₁*dHˣdξ-x₁₂*dHˣdη+y₃₁*+y₁₂]
        D=D₀*[1 ν₀ 0;
              ν₀ 1 0;
              0  0 (1-ν₀)/2]
        return transpose(B)*D*B
    end
    K=quad_hammer(BtDB,[0,1,0],[0,0,1],3)
    #9x9 to 18x18
    I=1:9
    J=[3,4,5,9,10,11,15,16,17]
    L=sparse(I,J,1.,9,18)
    return L'*K*L
end

function M_DKT(elm::Tria)::Matrix{Float64}
    ρ=elm.material.ρ
    A=elm.A
    t=elm.t
    W=ρ*A*t
    M=zeros(18,18)
    for i in [3,9,15]
        M[i,i]=W/3
    end
    return M
end

function P_DKT(elm::Tria,p::Float64)::Vector{Float64}
    A=elm.A
    P=zeros(18)
    P[3]=P[9]=P[15]=p*A/3
end

function W_TMGT(elm::Tria,u::Vector{Float64})::Vector{Float64}

end

function f_TMGT(elm::Tria,u::Vector{Float64})::Vector{Float64}
end
