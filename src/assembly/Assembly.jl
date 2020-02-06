
module FEAssembly

using SparseArrays
using LinearAlgebra
import Base.Filesystem

using ..FEStructure
using ..FESparse
using ..LoadCase

export assemble!,assemble_cable_link!,
clear_result!

mutable struct LCNode
    case
    children::Vector{LCNode}
end

mutable struct Assembly
    structure::Structure
    lcset::LoadCaseSet

    node_count::Int64
    beam_count::Int64
    quad_count::Int64
    tria_count::Int64

    nDOF::Int64
    nfreeDOF::Int64
    restrainedDOFs::Vector{Int64}

    lc_tree::LCNode
    working_path::String
end

function DFSfind(lcnode::LCNode,lc::String)
    if lcnode.case !="root" && lcnode.case.id==lc
        return lcnode
    elseif isempty(lcnode.children)
        return nothing
    else
        for child in lcnode.children
            res=DFSfind(child,lc)
            if res isa LCNode
                return res
            end
        end
    end
end

function idxmap(node::Node)
    i=node.hid
    collect(6*i-5:6*i)
end

function idxmap(elm::Beam)
    i = elm.node1.hid
    j = elm.node2.hid
    [6i-5:6i;6j-5:6j]
end

function idxmap(elm::Quad)
    i = elm.node1.hid
    j = elm.node2.hid
    k = elm.node3.hid
    l = elm.node4.hid
    [6i-5:6i;6j-5:6j;6k-5:6k;6l-5:6l]
end

function idxmap(elm::Tria)
    i = elm.node1.hid
    j = elm.node2.hid
    k = elm.node3.hid
    [6i-5:6i;6j-5:6j;6k-5:6k]
end


function disperse(A::Matrix,i::Vector{Int},N::Int)
    m,n=size(A)
    I = repeat(i,outer=n)
    J = repeat(i,inner=m)
    V = reshape(A,m*n)
    SparseMatrixCOO(N,N,I,J,V)
end

function disperse(A::Vector,i::Vector{Int},N::Int)
    # m=length(A)
    # I = i
    # J = [1 for i in 1:m]
    # V = A
    dropzeros!(SparseVector(N,i,A))
end

function appendcoo!(a::SparseMatrixCOO,b::SparseMatrixCOO)
    append!(a.rowptr,b.rowptr)
    append!(a.colptr,b.colptr)
    append!(a.nzval,b.nzval)
    return a
end

function appendcoo!(a::SparseVector,b::SparseVector)
    append!(a.nzind,b.nzind)
    append!(a.nzval,b.nzval)
    return a
end

assembleK(elm,N)=disperse(elm.T'*integrateK(elm)*elm.T,idxmap(elm),N)
assembleM(elm,N)=disperse(elm.T'*integrateM(elm)*elm.T,idxmap(elm),N)
assembleP(elm,N,f)=disperse(elm.T'*integrateP(elm,f),idxmap(elm),N)

function getvalues(dict::Dict,keys::Vector)
    v=[]
    for k in keys
        push!(v,dict[k])
    end
    return v
end

function calc_Pg(structure::Structure,gaccel::Float64,nDOF::Int)
    #calculate gravity
    _gset=LoadCaseSet()
    add_static_case!(_gset,"__Gravity__",1.0)
    Pg=spzeros(nDOF)
    for elm in values(structure.beams)
        #load on repeat id will be added together
        f=[0,0,-gaccel*elm.material.ρ*elm.section.A,0,0,0,
           0,0,-gaccel*elm.material.ρ*elm.section.A,0,0,0]
        T=sparse(elm.T)
        fᵉ=T*f
        fi1,fi2,fi3,mi1,mi2,mi3,fi1,fi2,fi3,mi1,mi2,mi3=fᵉ
        add_beam_distributed!(_gset,"__Gravity__",elm.id,fi1,fi2,fi3,mi1,mi2,mi3,fi1,fi2,fi3,mi1,mi2,mi3)
    end
    #add quad / tria gravity
    gforce=collect(values(_gset.statics["__Gravity__"].beam_forces))
    if !isempty(gforce)
        Pg=reduce(appendcoo!,assembleP.(getvalues(structure.beams,getproperty.(gforce,:id)),nDOF,gforce),init=Pg)
    end
    gforce=collect(values(_gset.statics["__Gravity__"].quad_forces))
    if !isempty(gforce)
        Pg=reduce(appendcoo!,assembleP.(getvalues(structure.quads,getproperty.(gforce,:id)),nDOF,gforce),init=Pg)
    end
    gforce=collect(values(_gset.statics["__Gravity__"].tria_forces))
    if !isempty(gforce)
        Pg=reduce(appendcoo!,assembleP.(getvalues(structure.trias,getproperty.(gforce,:id)),nDOF,gforce),init=Pg)
    end
    return to_array(Pg)
end

"""
    assemble!(structure,lcset;path=pwd())
集成结构与工况
# 参数
- `structure::Structure`: Structure类型实例
- `lcset::LoadCaseSet`: LoadCaseSet类型实例
# 关键字参数
- `mass_source::String`:`weight` or `loadcases`
- `mass_cases::Array{String}`: 质量荷载工况
- `mass_case_factors::Array{Number}`:工况叠加系数
- `path`: 工作路径
# 返回
- `assembly`: Assembly对象
"""
function assemble!(structure,lcset;mass_source="weight",mass_cases=[],mass_cases_factors=[],path=pwd())
    nDOF=length(structure.nodes)*6
    restrainedDOFs=[]
    node_count=length(structure.nodes)
    beam_count=length(structure.beams)
    quad_count=length(structure.quads)
    tria_count=length(structure.trias)

    K=spzeros_coo(nDOF,nDOF)
    M=spzeros_coo(nDOF,nDOF)
    C=sparse(0.02*I,nDOF,nDOF)

    lc_tree=LCNode("root",[])
    lcs=merge(lcset.statics,lcset.modals,lcset.bucklings,lcset.time_histories,lcset.response_spectrums)

    while !isempty(lcs)
        for k in keys(lcs)
            if lcs[k].plc==""
                lcnode=LCNode(pop!(lcs,k),[])
                push!(lc_tree.children,lcnode)#root
            else
                lcnode=DFSfind(lc_tree,lcs[k].plc)
                if lcnode isa LCNode
                    child=LCNode(pop!(lcs,k),[])
                    push!(lcnode.children,child)
                else
                    continue
                end
            end
        end
    end

    for node in values(structure.nodes)
        idx=node.hid*6-6
        for r in node.restraints
            idx+=1
            if r
                push!(restrainedDOFs,idx)
            end
        end
    end

    sort!(restrainedDOFs)
    mask=[(i in restrainedDOFs) ? false : true for i in 1:nDOF]

@time begin
    if node_count>0
        K=reduce(appendcoo!,assembleK.(values(structure.nodes),nDOF))
        M=reduce(appendcoo!,assembleM.(values(structure.nodes),nDOF))
    end
end

@time begin
    if beam_count>0
        K=reduce(appendcoo!,assembleK.(values(structure.beams),nDOF))
        M=reduce(appendcoo!,assembleM.(values(structure.beams),nDOF))
    end
end

@time begin
    if quad_count>0
        K=reduce(appendcoo!,assembleK.(values(structure.quads),nDOF))
        M=reduce(appendcoo!,assembleM.(values(structure.quads),nDOF))
    end
end

@time begin
    if tria_count>0
        K=reduce(appendcoo!,assembleK.(values(structure.trias),nDOF))
        M=reduce(appendcoo!,assembleM.(values(structure.trias),nDOF))
    end
end

    if structure.damp=="constant"
        C=sparse(structure.ζ₁*I,nDOF,nDOF)
    elseif structure.damp=="rayleigh"
        C=structure.ζ₁*K+structure.ζ₂*M
    else
        throw("Structural damping error!")
    end

    structure.K=to_csc(K)
    structure.C=C

    Pg=calc_Pg(structure,9.81,nDOF)
    # for load in values(_gset.statics["__Gravity__"].beam_forces)
    #     elm=structure.beams[load.id]
    #     Kᵉ=integrateK(elm)
    #     Pᵉ=integrateP(elm,load)
    #     rDOF=findall(x->x==true,elm.release)
    #     K̃ᵉ,P̃ᵉ=FEStructure.static_condensation(Kᵉ,Pᵉ,rDOF)
    #     idx=idxmap(elm)
    #     Pg+=disperse(reshape(elm.T'*P̃ᵉ,12),idx,nDOF)
    # end

    cases=merge(lcset.statics,lcset.bucklings,lcset.time_histories)
    for l in keys(cases)
        loadcase=cases[l]
        P=zeros(nDOF)
        if haskey(lcset.statics,l)
            P.+=Pg*loadcase.gfactor
        end

        for load in values(loadcase.beam_forces)
            elm=structure.beams[load.id]
            Kᵉ=integrateK(elm)
            Pᵉ=integrateP(elm,load)
            rDOF=findall(x->x==true,elm.release)
            K̃ᵉ,P̃ᵉ=FEStructure.static_condensation(Kᵉ,Pᵉ,rDOF)
            idx=idxmap(elm)
            P+=Array(disperse(reshape(elm.T'*P̃ᵉ,12),idx,nDOF))
            # println(P)
            # println("///")
        end
        for load in values(loadcase.nodal_forces)
            node=structure.nodes[load.id]
            Pⁿ=reshape(node.T'*load.val,6)
            idx=idxmap(node)
            P+=Array(disperse(Pⁿ,idx,nDOF))
        end
        cases[l].P=P
    end

    if mass_source=="weight"
        structure.M=to_csc(M)
    elseif mass_source=="loadcases"
        M=spzeros_coo(nDOF,nDOF)
        for (lc,fac) in (mass_cases, mass_case_factors)
            M+=Diagonal(lcset.statics[lc].P)*fac
        end
        structure.M=sparse(M)
    else
        throw("mass_source should be only weight or loadcases")
    end
    return Assembly(structure,lcset,node_count,beam_count,quad_count,tria_count,
                    nDOF,nDOF-length(restrainedDOFs),restrainedDOFs,lc_tree,path)
end

function clear_result!(assembly)
    nDOF=assembly.nDOF
    assembly.structure.K=spzeros(nDOF,nDOF)
    assembly.structure.K̄=spzeros(nDOF,nDOF)
    assembly.structure.M=spzeros(nDOF,nDOF)
    assembly.structure.M̄=spzeros(nDOF,nDOF)
    assembly.structure.C=spzeros(nDOF,nDOF)
    assembly.structure.C̄=spzeros(nDOF,nDOF)
    for i in keys(assembly.lcset.statics)
        assembly.lcset.statics[i].P=zeros(nDOF)
        assembly.lcset.statics[i].P̄=zeros(nDOF)
    end
    assembly.restrainedDOFs=[]
    result_dir=joinpath(assembly.working_path,".analysis")
    rm(result_dir,recursive=true)
    return
end

end
