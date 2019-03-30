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

x=symbols("x₁ x₂ x₃")
y=symbols("y₁ y₂ y₃")
xᵢⱼ=[x[i]-x[j] for (i,j) in zip([1,2,3],[2,3,1])]
yᵢⱼ=[y[i]-y[j] for (i,j) in zip([1,2,3],[2,3,1])]
lᵢⱼ²=[xᵢⱼ.^2 .+ yᵢⱼ.^2]

A,t=symbols("A,t")
a=symbols("a₁ a₂ a₃ a₄ a₅ a₆ a₇ a₈")
b=symbols("b₁ b₂ b₃ b₄ b₅ b₆ b₇ b₈")
c=symbols("c₁ c₂ c₃ c₄ c₅ c₆ c₇ c₈")
d=symbols("d₁ d₂ d₃ d₄ d₅ d₆ d₇ d₈")
e=symbols("e₁ e₂ e₃ e₄ e₅ e₆ e₇ e₈")
E₀,ν₀=symbols("E₀ ν₀")


D₀=E₀*t^3/12(1-ν₀^2)
D=D₀*[1 ν₀ 0;
      ν₀ 1 0;
      0  0 (1-ν₀)/2]

Hˣ=Array{Basic}(undef,9,1)
Hʸ=Array{Basic}(undef,9,1)
ξᵢ=[1,-1,-1,1]
ηᵢ=[1,1,-1,-1]
ξ,η=symbols("ξ η")
Nᵢ[1]=2*(1-ξ-η)*(0.5-ξ-η)
Nᵢ[2]=ξ*(2ξ-1)
Nᵢ[3]=η*(2η-1)
Nᵢ[4]=4*ξ*η
Nᵢ[5]=4*η*(1-ξ-η)
Nᵢ[6]=4*ξ*(1-ξ-η)

Hˣ[1]=1.5*(a[6]*Nᵢ[6]-a[5]*Nᵢ[5])
Hˣ[2]=b[5]*Nᵢ[5]+b[6]*b[6]
Hˣ[3]=Nᵢ[1]-c[5]*Nᵢ[5]-c[6]*Nᵢ[6]
Hʸ[1]=1.5*(d[6]*Nᵢ[6]-d[5]*Nᵢ[5])
Hʸ[2]=-Nᵢ[1]+e[5]*Nᵢ[5]+e[6]*Nᵢ[6]
Hʸ[3]=-Hˣ[2]

Hˣ[4]=1.5*(a[4]*Nᵢ[4]-a[6]*Nᵢ[6])
Hˣ[5]=b[6]*Nᵢ[6]+b[4]*b[4]
Hˣ[6]=Nᵢ[2]-c[6]*Nᵢ[6]-c[4]*Nᵢ[4]
Hʸ[4]=1.5*(d[4]*Nᵢ[4]-d[6]*Nᵢ[6])
Hʸ[5]=-Nᵢ[2]+e[6]*Nᵢ[6]+e[4]*Nᵢ[4]
Hʸ[6]=-Hˣ[2]

Hˣ[7]=1.5*(a[5]*Nᵢ[5]-a[4]*Nᵢ[4])
Hˣ[8]=b[4]*Nᵢ[4]+b[5]*b[5]
Hˣ[9]=Nᵢ[1]-c[4]*Nᵢ[4]-c[5]*Nᵢ[5]
Hʸ[7]=1.5*(d[5]*Nᵢ[5]-d[4]*Nᵢ[4])
Hʸ[8]=-Nᵢ[3]+e[4]*Nᵢ[4]+e[5]*Nᵢ[5]
Hʸ[9]=-Hˣ[2]

dHˣdξ=diff.(Hˣ,ξ)
dHˣdη=diff.(Hˣ,η)
dHʸdξ=diff.(Hʸ,ξ)
dHʸdη=diff.(Hʸ,η)

@show dHˣdξ
@show dHˣdη
@show dHʸdξ
@show dHʸdη

x₁₂,x₃₁=xᵢⱼ[1],xᵢⱼ[3]
y₁₂,y₃₁=yᵢⱼ[1],yᵢⱼ[3]

y₃₁.*dHˣdξ
B=[y₃₁*transpose(dHˣdξ)+y₁₂*transpose(dHˣdη);
  -x₃₁*transpose(dHʸdξ)-y₁₂*transpose(dHʸdη);
  -x₃₁*transpose(dHˣdξ)-x₁₂*transpose(dHˣdη)+y₃₁*transpose(dHʸdξ)+y₁₂*transpose(dHʸdη)]/2/A

K=transpose(B)*D*B

open("./k_DKT.jl","w+") do f
  for i in 1:size(K,1)
      for j in i:size(K,2)
          if i==j
              write(f,"K["*string(i)*","*string(j)*"]="*string(K[i,j])*"\n")
          elseif string(K[i,j])!="0"
              write(f,"K["*string(i)*","*string(j)*"]="*"K["*string(j)*","*string(i)*"]="*string(K[i,j])*"\n")
          end
      end
  end
end
