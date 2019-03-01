
module AssemblyModule

using SparseArrays
using LinearAlgebra
import Base.Filesystem

using ..FEStructure
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

    node_count::Int
    beam_count::Int
    quad_count::Int

    nDOF::Int
    nfreeDOF::Int
    restrainedDOFs::Vector

    lc_tree::LCNode
    working_path::String
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

    K=spzeros(nDOF,nDOF)
    K̄=spzeros(nDOF,nDOF)
    M=spzeros(nDOF,nDOF)
    M̄=spzeros(nDOF,nDOF)
    C=sparse(0.02*I,nDOF,nDOF)
    C̄=sparse(0.02*I,nDOF,nDOF)

    lc_tree=LCNode("root",[])
    lcs=merge(lcset.statics,lcset.modals,lcset.bucklings,lcset.time_histories,lcset.response_spectrums)
    function DFSfind(lcnode,lc::String)
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
# @show 2
# @time begin
    for elm in values(structure.beams)
        i = elm.node1.hid
        j = elm.node2.hid
        T=sparse(elm.T)

        I=collect(1:12)
        J=[6i-5:6i;6j-5:6j]
        G=sparse(I,J,1.0,12,nDOF)

        Kᵉ=integrateK!(elm)
        Mᵉ=integrateM!(elm)

        A=T*G
        K+=A'*Kᵉ*A
        M+=A'*Mᵉ*A
    end

    for elm in values(structure.quads)
        i = elm.node1.hid
        j = elm.node2.hid
        k = elm.node3.hid
        l = elm.node4.hid

        T=sparse(elm.T)

        I=collect(1:24)
        J=[6i-5:6i;6j-5:6j;6k-5:6k;6l-5:6l]
        G=sparse(I,J,1.0,24,nDOF)

        Kᵉ=integrateK!(elm)
        Mᵉ=integrateM!(elm)

        A=T*G
        K+=A'*Kᵉ*A
        M+=A'*Mᵉ*A
    end
# end
# @show 3
    for node in values(structure.nodes)
        T=node.T
        spring=T'*node.spring
        mass=T'*node.mass
        idx=node.hid*6-6
        Kⁿ=sparse(Diagonal(spring))
        Mⁿ=sparse(Diagonal(mass))

        i=node.hid
        I=collect(1:6)
        J=collect(6i-5:6i)
        G=sparse(I,J,1.0,6,nDOF)

        A=T*G
        K+=A'*Kⁿ*A
        M+=A'*Mⁿ*A
    end

    if structure.damp=="constant"
        C=sparse(structure.ζ₁*I,nDOF,nDOF)
    elseif structure.damp=="rayleigh"
        C=structure.ζ₁*K+structure.ζ₂*M
    else
        throw("Structural damping error!")
    end
# @show 4
    structure.K=K
    structure.C=C

    #calculate gravity
    _gset=LoadCaseSet()
    add_static_case!(_gset,"__Gravity__",1.0)
    Pg=spzeros(nDOF,1)
    for elm in values(structure.beams)
        #load on repeat id will be added together
        f=[0,0,-9.81*elm.material.ρ*elm.section.A,0,0,0,
           0,0,-9.81*elm.material.ρ*elm.section.A,0,0,0]
        T=sparse(elm.T)
        fᵉ=T*f
        fi1,fi2,fi3,mi1,mi2,mi3,fi1,fi2,fi3,mi1,mi2,mi3=fᵉ
        add_beam_distributed!(_gset,"__Gravity__",elm.id,fi1,fi2,fi3,mi1,mi2,mi3,fi1,fi2,fi3,mi1,mi2,mi3)
    end
    for load in values(_gset.statics["__Gravity__"].beam_forces)
        elm=structure.beams[load.id]
        i = elm.node1.hid
        j = elm.node2.hid
        T=sparse(elm.T)
        I=collect(1:12)
        J=[6i-5:6i;6j-5:6j]
        G=sparse(I,J,1.0,12,nDOF)
        Pᵉ=integrateP!(elm,load)
        Pg+=G'*T'*Pᵉ
    end
# @show 5
    #calculate loadcase loads
    cases=merge(lcset.statics,lcset.bucklings,lcset.time_histories)
    for l in keys(cases)
        loadcase=cases[l]
        P=spzeros(nDOF,1)
        if l in keys(lcset.statics)
            P+=Pg*loadcase.gfactor
        end
        for load in values(loadcase.beam_forces)
            elm=structure.beams[load.id]
            i = elm.node1.hid
            j = elm.node2.hid
            T=sparse(elm.T)
            I=collect(1:12)
            J=[6i-5:6i;6j-5:6j]
            G=sparse(I,J,1.0,12,nDOF)
            Pᵉ=integrateP!(elm,load)
            P+=G'*T'*Pᵉ
        end
# @show 6
        for load in values(loadcase.nodal_forces)
            node=structure.nodes[load.id]
            i = node.hid
            T=sparse(node.T)
            I=collect(1:6)
            J=collect(6i-5:6i)
            G=sparse(I,J,1.0,6,nDOF)
            Pⁿ=reshape(load.val,6,1)
            P+=G'*T'*Pⁿ
        end
        cases[l].P=Array(P)[:,1]
    end
# @show 7
    if mass_source=="weight"
        structure.M=M
    elseif mass_source=="loadcases"
        M=zero(M)
        for (lc,fac) in (mass_cases, mass_case_factors)
            M+=Diagonal(lcset.statics[lc].P)*fac
        end
        structure.M=M
    else
        throw("mass_source should be only weight or loadcases")
    end
    return Assembly(structure,lcset,node_count,beam_count,quad_count,
    nDOF,nDOF-length(restrainedDOFs),restrainedDOFs,lc_tree,path)
end


function assemble_cable_link!(assembly)
    structure=assembly.structure
    loadcase=assembly.loadcase
    nDOF=assembly.nDOF
    K=assembly.K
    M=assembly.M
    C=assembly.C
    P=spzeros(nDOF,1)

    if !assembly.is_assembled #only assemble K & M when the structure it is not yet assembled
        for elm in structure.beams
            i = elm.node1.hid
            j = elm.node2.hid
            T=sparse(elm.T)
            G=zeros(12,nDOF)
            G[1:6,6i-5:6i]=Matrix(1.0I,6,6)
            G[7:12,6j-5:6j]=Matrix(1.0I,6,6)
            G=sparse(G)
            Kᵉ=integrateK!(elm)
            Mᵉ=integrateM!(elm)

            K+=G'*T'*Kᵉ*T*G
            M+=G'*T'*Mᵉ*T*G
        end

        for node in structure.nodes
            T=node.T
            spring=T*node.spring
            mass=T*node.mass
            res=node.restraints
            idx=node.hid*6-6
            for i in 1:6
                if spring[i]!=0
                    K[idx+i,idx+i]+=spring[i]
                end
                if mass[i]!=0
                    M[idx+i,idx+i]+=mass[i]
                end
            end
        end

        assembly.K=K
        assembly.M=M
        assembly.is_assembled=true
    end

    if loadcase.gfactor!=0.
        for elm in structure.cables
            #load on repeat id will be added together
            f=[0,0,-loadcase.gfactor*9.81*elm.ρ*elm.A,
               0,0,-loadcase.gfactor*9.81*elm.ρ*elm.A]
            T=elm.T
            fᵉ=T'*f
            fi1,fi2,fi3,fi1,fi2,fi3=fᵉ
            add_cable_distributed!(loadcase,elm.id,fi1,fi2,fi3,fi1,fi2,fi3)
        end
        for elm in structure.links
            #load on repeat id will be added together
            f=[0,0,-loadcase.gfactor*9.81*elm.ρ*elm.A,
               0,0,-loadcase.gfactor*9.81*elm.ρ*elm.A]
            T=elm.T
            fᵉ=T'*f
            fi1,fi2,fi3,fi1,fi2,fi3=fᵉ
            add_link_distributed!(loadcase,elm.id,fi1,fi2,fi3,fi1,fi2,fi3)
        end
        for elm in structure.beams
            #load on repeat id will be added together
            f=[0,0,-loadcase.gfactor*9.81*elm.ρ*elm.A,0,0,0,
               0,0,-loadcase.gfactor*9.81*elm.ρ*elm.A,0,0,0]
            T=elm.T
            fᵉ=T'*f
            fi1,fi2,fi3,mi1,mi2,mi3,fi1,fi2,fi3,mi1,mi2,mi3=fᵉ
            add_beam_distributed!(loadcase,elm.id,fi1,fi2,fi3,mi1,mi2,mi3,fi1,fi2,fi3,mi1,mi2,mi3)
        end
    end

    for load in loadcase.cable_forces
        idx=findfirst(x->x.id==load.id,structure.cables)
        if idx isa Nothing
            throw("load on an invalid cable with id "*string(load.id)*"!")
        end
        elm=structure.cables[idx]
        i = elm.node1.hid
        j = elm.node2.hid
        L=spzeros(6,12)
        L[1:3,1:3]=Matrix(1.0I,3,3)
        L[4:6,7:9]=Matrix(1.0I,3,3)
        T=sparse(L'*elm.T*L)
        G=spzeros(12,nDOF)
        G[1:6,6i-5:6i]=Matrix(1.0I,6,6)
        G[7:12,6j-5:6j]=Matrix(1.0I,6,6)
        G=sparse(G)
        Pᵉ=integrateP!(elm,load)
        P+=G'*T'*Pᵉ
    end

    for load in loadcase.link_forces
        idx=findfirst(x->x.id==load.id,structure.links)
        if idx isa Nothing
            throw("load on an invalid link with id "*string(load.id)*"!")
        end
        elm=structure.links[idx]
        i = elm.node1.hid
        j = elm.node2.hid
        L=spzeros(6,12)
        L[1:3,1:3]=Matrix(1.0I,3,3)
        L[4:6,7:9]=Matrix(1.0I,3,3)
        T=sparse(L'*elm.T*L)
        G=spzeros(12,nDOF)
        G[1:6,6i-5:6i]=Matrix(1.0I,6,6)
        G[7:12,6j-5:6j]=Matrix(1.0I,6,6)
        G=sparse(G)
        Pᵉ=integrateP!(elm,load)
        P+=G'*T'*Pᵉ
    end

    for load in loadcase.beam_forces
        idx=findfirst(x->x.id==load.id,structure.beams)
        if idx isa Nothing
            throw("load on an invalid beam with id "*string(load.id)*"!")
        end
        elm=structure.beams[idx]
        i = elm.node1.hid
        j = elm.node2.hid
        T=sparse(elm.T)
        G=zeros(12,nDOF)
        G[1:6,6i-5:6i]=Matrix(1.0I,6,6)
        G[7:12,6j-5:6j]=Matrix(1.0I,6,6)
        G=sparse(G)
        Pᵉ=integrateP!(elm,load)
        P+=G'*T'*Pᵉ
    end

    for load in loadcase.nodal_forces
        idx=findfirst(x->x.id==load.id,structure.nodes)
        if idx isa Nothing
            throw("load on an invalid node with id "*string(load.id)*"!")
        end
        node=structure.nodes[idx]
        i = node.hid
        T=sparse(node.T)
        G=zeros(6,nDOF)
        G[1:6,6i-5:6i]=Matrix(1.0I,6,6)
        G=sparse(G)
        Pⁿ=reshape(load.val,6,1)
        P+=G'*T'*Pⁿ
    end

    assembly.P=P

    #find free nodes
    all_nodes=Set()
    beam_nodes=Set()
    for node in structure.nodes
        push!(all_nodes,node.hid)
    end
    for elm in structure.beams
        i=elm.node1.hid
        j=elm.node2.hid
        push!(beam_nodes,i,j)
    end
    free_nodes=collect(setdiff(all_nodes,beam_nodes))

    free_DOFs=[]
    for hid in free_nodes
        for i in 1:6
            push!(free_DOFs,hid*6-6+i)
        end
    end
    for node in structure.nodes
        idx=node.hid*6-6
        for r in node.restraints
            idx+=1
            if r
                push!(assembly.restrainedDOFs,idx)
            end
        end
    end

    restrainedDOFs=collect(Set([free_DOFs;assembly.restrainedDOFs]))
    assembly.restrainedDOFs=[]
    assembly.K̄=assembly.K
    assembly.P̄=assembly.P
    for i in restrainedDOFs
        assembly.K̄[i,i]+=1e20
        assembly.P̄[i]+=1e20*assembly.d[i]
    end
    return
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
    # result_dir=chmod(result_dir, 777, recursive=true)
    rm(result_dir,recursive=true)
    # for f in Filesystem.readdir(result_dir)
    #     file=Filesystem.joinpath(result_dir,f)
    #     Filesystem.rm(file,force=true)
    # end
    return
end

end
