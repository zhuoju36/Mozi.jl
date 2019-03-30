using SymEngine
using HCubature
using LinearAlgebra

function L₂(X;x=symbols("x"),y=symbols("y"))
    res=Matrix{SymEngine.Basic}(undef,3,size(X,2))
    for i in 1:size(X,2)
        res[1,i]=diff(X[1,i],x)
        res[2,i]=diff(X[2,i],y)
        res[3,i]=diff(X[1,i],y)+diff(X[2,i],x)
    end
    return res
end

x=symbols("x₁ x₂ x₃ x₄")
y=symbols("y₁ y₂ y₃ y₄")
xᵢⱼ=[x[i]-x[j] for (i,j) in zip([1,2,3,4],[2,3,4,1])]
yᵢⱼ=[y[i]-y[j] for (i,j) in zip([1,2,3,4],[2,3,4,1])]
lᵢⱼ²=[xᵢⱼ.^2 .+ yᵢⱼ.^2]

Hˣ=Array{Basic}(undef,12,1)
Hʸ=Array{Basic}(undef,12,1)

a=symbols("a₁ a₂ a₃ a₄ a₅ a₆ a₇ a₈")
b=symbols("b₁ b₂ b₃ b₄ b₅ b₆ b₇ b₈")
c=symbols("c₁ c₂ c₃ c₄ c₅ c₆ c₇ c₈")
d=symbols("d₁ d₂ d₃ d₄ d₅ d₆ d₇ d₈")
e=symbols("e₁ e₂ e₃ e₄ e₅ e₆ e₇ e₈")

ξᵢ=[1,-1,-1,1]
ηᵢ=[1,1,-1,-1]
ξ,η=symbols("ξ η")

Nᵢ=Array{Basic}(undef,8,1)
Nᵢ[1:4]=1/4*(1 .+ ξᵢ*ξ).*(1 .+ ηᵢ*η)
Nᵢ[5]=0.5*(1-ξ^2)*(1-η)
Nᵢ[6]=0.5*(1-η^2)*(1+ξ)
Nᵢ[7]=0.5*(1-ξ^2)*(1+η)
Nᵢ[8]=0.5*(1-η^2)*(1-ξ)

Hˣ[1]=1.5*(a[5]*Nᵢ[5]-a[8]*Nᵢ[8])
Hˣ[2]=b[5]*Nᵢ[5]+b[8]*b[8]
Hˣ[3]=Nᵢ[1]-c[5]*Nᵢ[5]-c[8]*Nᵢ[8]
Hʸ[1]=1.5*(d[5]*Nᵢ[5]-d[8]*Nᵢ[8])
Hʸ[2]=-Nᵢ[1]+e[5]*Nᵢ[5]+e[8]*Nᵢ[8]
Hʸ[3]=b[5]*Nᵢ[5]-b[8]*Nᵢ[8]

Hˣ[4]=1.5*(a[6]*Nᵢ[6]-a[5]*Nᵢ[5])
Hˣ[5]=b[6]*Nᵢ[6]+b[5]*b[5]
Hˣ[6]=Nᵢ[2]-c[6]*Nᵢ[6]-c[5]*Nᵢ[5]
Hʸ[4]=1.5*(d[6]*Nᵢ[6]-d[5]*Nᵢ[5])
Hʸ[5]=-Nᵢ[2]+e[6]*Nᵢ[6]+e[5]*Nᵢ[5]
Hʸ[6]=b[6]*Nᵢ[6]-b[5]*Nᵢ[5]

Hˣ[7]=1.5*(a[7]*Nᵢ[7]-a[6]*Nᵢ[6])
Hˣ[8]=b[7]*Nᵢ[7]+b[6]*b[6]
Hˣ[9]=Nᵢ[1]-c[7]*Nᵢ[7]-c[6]*Nᵢ[6]
Hʸ[7]=1.5*(d[7]*Nᵢ[7]-d[6]*Nᵢ[6])
Hʸ[8]=-Nᵢ[1]+e[7]*Nᵢ[7]+e[6]*Nᵢ[6]
Hʸ[9]=b[7]*Nᵢ[7]-b[6]*Nᵢ[6]

Hˣ[10]=1.5*(a[8]*Nᵢ[8]-a[7]*Nᵢ[7])
Hˣ[11]=b[8]*Nᵢ[8]+b[7]*b[7]
Hˣ[12]=Nᵢ[1]-c[8]*Nᵢ[8]-c[7]*Nᵢ[7]
Hʸ[10]=1.5*(d[8]*Nᵢ[8]-d[7]*Nᵢ[7])
Hʸ[11]=-Nᵢ[1]+e[8]*Nᵢ[8]+e[7]*Nᵢ[7]
Hʸ[12]=b[8]*Nᵢ[8]-b[7]*Nᵢ[7]

dHˣdξ=diff.(Hˣ,ξ)
dHˣdη=diff.(Hˣ,η)
dHʸdξ=diff.(Hʸ,ξ)
dHʸdη=diff.(Hʸ,η)

x₁,x₂,x₃,y₁,y₂,y₃,x₄,y₄=symbols("x₁ x₂ x₃ y₁ y₂ y₃ x₄ y₄")
ξ,η=symbols("ξ η")
E₀,ν₀=symbols("E₀ ν₀")

ξ₁,η₁=-1,-1
ξ₂,η₂=1,-1
ξ₃,η₃=1,1
ξ₄,η₄=-1,1
xᵢ=[x₁,x₂,x₃,x₄]
yᵢ=[y₁,y₂,y₃,y₄]
ξᵢ=[ξ₁,ξ₂,ξ₃,ξ₄]
ηᵢ=[η₁,η₂,η₃,η₄]
Nᵢ=1/4*(1 .+ ξᵢ*ξ).*(1 .+ ηᵢ*η)

x₁,x₂,x₃,y₁,y₂,y₃,x₄,y₄=symbols("x₁ x₂ x₃ y₁ y₂ y₃ x₄ y₄")
X=[x₁ y₁;x₂ y₂;x₃ y₃;x₄ y₄]
J=Matrix{SymEngine.Basic}(undef,2,2)
J[1,1]=sum(diff.(Nᵢ,ξ).*X[:,1])
J[1,2]=sum(diff.(Nᵢ,ξ).*X[:,2])
J[2,1]=sum(diff.(Nᵢ,η).*X[:,1])
J[2,2]=sum(diff.(Nᵢ,η).*X[:,2])
open("./jacobi.jl","w+") do f
    for i in 1:size(J,1)
        for j in 1:size(J,2)
            # eval(Meta.parse("K["*string(i)*","*string(j)*"]="*string(K[i,j])))
            write(f,"J["*string(i)*","*string(j)*"]="*string(J[i,j])*"\n")
        end
    end
end
