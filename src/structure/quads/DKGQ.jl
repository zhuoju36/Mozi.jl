function K_DKGQ(elm::Quad)::SparseMatrixCSC{Float64}
    Kᵉ=spzeros(24,24)
    I=1:12
    J=[1,2,6,7,8,12,13,14,18,19,20,24]
    L=sparse(I,J,1.,12,24)
    Km=K_GQ12(elm)
    Kᵉ+=L'*Km*L

    I=1:12
    J=[3,4,5,9,10,11,15,16,17,21,22,23]
    L=sparse(I,J,1.,12,24)
    Kb=K_DKQ(elm)
    Kᵉ+=L'*Kb*L
    return sparse(Kᵉ)
end
