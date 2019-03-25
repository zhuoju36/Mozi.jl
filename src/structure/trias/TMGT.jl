function K_TMGT(elm::Tria)::SparseMatrixCSC{Float64}
    return K_GT9(elm)+K_TMT(elm)
end
