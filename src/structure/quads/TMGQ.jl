function K_TMGQ(elm::Quad)::SparseMatrixCSC{Float64}
    return K_GQ12(elm)+K_TMQ(elm)
end
