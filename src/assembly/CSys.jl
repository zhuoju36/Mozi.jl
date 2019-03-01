module CSysModule

using LinearAlgebra

export CSys,get_transform_matrix

struct CSys
    O::Array{Float64}
    x::Array{Float64}
    y::Array{Float64}
    z::Array{Float64}
    T::Matrix{Float64} #transorm_matrix
end

function CSys(o::Array,p₁::Array,p₂::Array)
    v₁=p₁-o
    v₂=p₂-o
    if abs(v₁⋅v₂/BLAS.nrm2(v₁)/BLAS.nrm2(v₂))==1
        error("Two vectors should not be parallel!")
    end
    x=v₁/BLAS.nrm2(v₁)
    z=v₁×v₂
    z=z/BLAS.nrm2(z)
    y=z×x
    T=[reshape(x,1,3);
        reshape(y,1,3);
        reshape(z,1,3)]
    CSys(o,x,y,z,T)
end

end
