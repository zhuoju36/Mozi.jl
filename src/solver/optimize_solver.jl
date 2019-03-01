using LinearAlgebra
import JuMP
import Clp
import Convex, SCS

export optimize_cable_force, optimize_cable_force_flexible

function Laplcian()
    """Provides a scripting component.
    Inputs:
        x: The x script variable
        y: The y script variable
    Output:
        a: The a output variable"""

__author__ = "Dell"
__version__ = "2019.02.19"

import rhinoscriptsyntax as rs

tol=1e-3

def dist(x,y):
    s=0
    s+=(x.X-y.X)**2
    s+=(x.Y-y.Y)**2
    s+=(x.Z-y.Z)**2
    return s**0.5

def find(b,a):
    for i in a:
        if dist(i,b)<tol:
            return True
    else:
        return False

def get(b,a):
    for i in range(len(a)):
        if dist(a[i],b)<tol:
            return i

nodes=[]
links=[]
neibours={}

for i in x+y:
    if not find(i.PointAtStart,nodes):
        neibours[len(nodes)]=[]
        nodes.append(i.PointAtStart)

    if not find(i.PointAtEnd,nodes):
        neibours[len(nodes)]=[]
        nodes.append(i.PointAtEnd)


for i in x:
    n1=get(i.PointAtStart,nodes)
    n2=get(i.PointAtEnd,nodes)
    for w in range(50):
        neibours[n1].append(nodes[n2])
        neibours[n2].append(nodes[n1])

for i in y:
    n1=get(i.PointAtStart,nodes)
    n2=get(i.PointAtEnd,nodes)
    neibours[n1].append(nodes[n2])
    neibours[n2].append(nodes[n1])


def update_loc(j):
    m=len(neibours[j])
    x=sum([n.X for n in neibours[j]])
    y=sum([n.Y for n in neibours[j]])
    z=sum([n.Z for n in neibours[j]])
    nodes[j].X=x/m
    nodes[j].Y=y/m
    nodes[j].Z=z/m

def is_boundary(node,z):
    for i in z:
        if dist(node,i)<tol:
            return True
    else:
        return False
maxiter=100

for i in range(maxiter):
    for j in range(len(nodes)):
        if not is_boundary(nodes[j],z):
            update_loc(j)

a=nodes

end

function find_boundary_nodes(structure)
    boundary_nodes=Set()
    restrainted_nodes=Set()
    cablelink_nodes=Set()
    beam_nodes=Set()
    for node in structure.nodes
        if true in node.restraints
            push!(restrainted_nodes,node.hid)
        end
    end
    for elm in structure.links
        i=elm.node1.hid
        j=elm.node2.hid
        push!(cablelink_nodes,i,j)
    end
    for elm in structure.cables
        i=elm.node1.hid
        j=elm.node2.hid
        push!(cablelink_nodes,i,j)
    end
    for elm in structure.beams
        i=elm.node1.hid
        j=elm.node2.hid
        push!(beam_nodes,i,j)
    end
    boundary_nodes=union(intersect(cablelink_nodes,beam_nodes), intersect(cablelink_nodes,restrainted_nodes))
    inner_nodes=setdiff(cablelink_nodes,boundary_nodes)
    return collect(boundary_nodes),collect(inner_nodes)
end

# function optimize_cable_force_flexible(assembly)
#     structure=assembly.structure
#     loadcase=assembly.loadcase
#     if !assembly.is_assembled
#         throw("The structure should be assembled with a loadcase first!")
#     end
#     boundary_nodes,inner_nodes=SolverModule.find_boundary_nodes(structure)
#     cablelink_nodes=[inner_nodes;boundary_nodes]
#     cables=structure.cables
#     links=structure.links
#     m₁=length(cables)
#     m₂=length(links)
#     m=m₁+m₂
#     n=length(cablelink_nodes)
#     n₂=length(boundary_nodes)
#     n₁=n-n₂
#     fᵉ=zeros(m)
#
#     Ts=[]
#     ns=[]
#     es=[]
#
#     for cable in structure.cables
#         push!(Ts,cable.T)
#         push!(es,(cable.node1.hid,cable.node2.hid))
#     end
#
#     for link in structure.links
#         push!(Ts,link.T)
#         push!(es,(link.node1.hid,link.node2.hid))
#     end
#
#     fₑ₁=zeros(3n₁)
#     fₑ₂=zeros(3n₂)
#     i=1
#     for hid in inner_nodes
#         fₑ₁[3i-2:3i]=assembly.P[6*hid-5:6*hid-3]
#         i+=1
#     end
#     i=1
#     for hid in boundary_nodes
#         fₑ₂[3i-2:3i]=assembly.P[6*hid-5:6*hid-3]
#         i+=1
#     end
#
#     S=zeros(3n,m)
#     for i=1:m
#         T=Ts[i]
#         L=[1;0;0;-1;0;0]
#         j=findfirst(x->x==es[i][1],cablelink_nodes)
#         k=findfirst(x->x==es[i][2],cablelink_nodes)
#         G=zeros(6,3n)
#         G[1:3,3j-2:3j]=Matrix(1.0I,3,3)
#         G[4:6 3k-2:3k]=Matrix(1.0I,3,3)
#         S[:,i]=G'*T'*L
#     end
#     S₁=S[1:3n₁,:]
#     S₂=S[3n₁+1:end,:]
#
#     w=zeros(3n₂)
#     k=0
#     for i in boundary_nodes
#         for j in 1:3
#             k+=1
#             w[k]=assembly.K[6*i-6+j,6*i-6+j]
#         end
#     end
#
#     for i in 1:length(boundary_nodes)
#         idx=boundary_nodes[i]
#         vec=zeros(3)
#         mask1=[c.node1.hid==idx for c in assembly.structure.cables]
#         cb1=assembly.structure.cables[mask1]
#         x1=[x.node2.loc-x.node1.loc for x in cb1]
#         if length(x1)!=0
#             vec+=sum(x1)
#         end
#         mask2=[c.node2.hid==idx for c in assembly.structure.cables]
#         cb2=assembly.structure.cables[mask2]
#         x2=[x.node1.loc-x.node2.loc for x in cb2]
#         if length(x2)!=0
#             vec+=sum(x2)
#         end
#         # mask1=[c.node1.hid==idx for c in assembly.structure.links]
#         # cb1=assembly.structure.links[mask1]
#         # x1=[x.node1.loc-x.node2.loc for x in cb1]
#         # if length(x1)!=0
#         #     vec+=sum(x1)
#         # end
#         # mask2=[c.node2.hid==idx for c in assembly.structure.links]
#         # cb2=assembly.structure.links[mask2]
#         # x2=[x.node2.loc-x.node1.loc for x in cb2]
#         # if length(x2)!=0
#         #     vec+=sum(x2)
#         # end
#
#         val=norm(vec)
#         if val!=0
#             vec./=norm(vec)
#         end
#         w[i*3-2]*=vec[1]
#         w[i*3-1]*=vec[2]
#         w[i*3]*=vec[3]
#     end
#
#     for i in length(w)
#         if w[i]!=0
#             w[i]=1/w[i]
#         end
#     end
#     w./=norm(w)
#
#     cᵀ=w'*[S₂ Matrix(1.0I,3n₂,3n₂)]
#     m₂=m-m₁
#     Z=zeros(m+3n₂,m+3n₂)
#     Z[1:m₁,1:m₁]=Matrix(1.0I,m₁,m₁)
#
#
#
#     println([x.id for x in assembly.structure.cables])
#     println([force/elm.A/elm.E for (force,elm) in zip(JuMP.getvalue(x)[1:m₁],assembly.structure.cables)])
#
#     println([x.id for x in assembly.structure.links])
#     println([force/elm.A/elm.E for (force,elm) in zip(JuMP.getvalue(x)[m₁+1:m],assembly.structure.links)])
# end

function optimize_cable_force(assembly)
    structure=assembly.structure
    loadcase=assembly.loadcase
    if !assembly.is_assembled
        assemble!(assembly)
    end
    boundary_nodes,inner_nodes=find_boundary_nodes(structure)
    cablelink_nodes=[inner_nodes;boundary_nodes]
    cables=structure.cables
    links=structure.links
    m₁=length(cables)
    m₂=length(links)
    m=m₁+m₂
    n=length(cablelink_nodes)
    n₂=length(boundary_nodes)
    n₁=n-n₂
    fᵉ=zeros(m)

    Ts=[] #transform matrix of link
    ns=[]
    es=[]

    for cable in structure.cables
        push!(Ts,cable.T)
        push!(es,(cable.node1.hid,cable.node2.hid))
    end

    for link in structure.links
        push!(Ts,link.T)
        push!(es,(link.node1.hid,link.node2.hid))
    end

    fₑ₁=zeros(3n₁)
    fₑ₂=zeros(3n₂)
    i=1
    for hid in inner_nodes
        idx=findfirst(x->x.id==hid,loadcase.nodal_forces)
        if typeof(idx)==Nothing continue end
        fₑ₁[3i-2:3i]=loadcase.nodal_forces[idx].val[1:3]
        i+=1
    end
    i=1
    for hid in boundary_nodes
        idx=findfirst(x->x.id==hid,loadcase.nodal_forces)
        if typeof(idx)==Nothing continue end
        fₑ₂[3i-2:3i]=loadcase.nodal_forces[idx].val[1:3]
        i+=1
    end

    S=zeros(3n,m)
    for i=1:m
        T=Ts[i]
        L=[1;0;0;-1;0;0]
        j=findfirst(x->x==es[i][1],cablelink_nodes)
        k=findfirst(x->x==es[i][2],cablelink_nodes)
        # T=zeros(6,6)
        # T[1:3,1:3]=T[4:6,4:6]=V
        G=zeros(6,3n)
        G[1:3,3j-2:3j]=Matrix(1.0I,3,3)
        G[4:6,3k-2:3k]=Matrix(1.0I,3,3)
        S[:,i]=G'*T'*L
        # println(structure.nodes[cablelink_nodes[j]].id,structure.nodes[cablelink_nodes[k]].id)
        # println(G*T*L*1)
    end
    S₁=S[1:3n₁,:]
    S₂=S[3n₁+1:end,:]
    w=ones(3n₂)
    cᵀ=w'*[S₂ Matrix(1.0I,3n₂,3n₂)]
    # S₁⋅fᵉ+f₁=zeros(n₁,1)
    m₂=m-m₁
    Z=zeros(m+3n₂,m+3n₂)
    Z[1:m₁,1:m₁]=Matrix(1.0I,m₁,m₁)
    # Z⋅x>=0
    # x>=x_max
    # x<=x_max
    opti_model = JuMP.Model(solver=Clp.ClpSolver())
    # 定义变量，注意这里使用了宏（macro），宏的调用也是Julia&JuMP高效编译/元编程(metaprogramming)的重要技巧
    x_min=-10000*ones(m+3n₂)
    x_max=10000*ones(m+3n₂)
    JuMP.@variable(opti_model, -10000 <= x[1:m+3n₂] <=10000)
    # 定义不等式约束
    JuMP.@constraints(opti_model, begin
        S₁*x[1:m]+fₑ₁.==0
        x[1:m₁] .>= 0
        x[1:m] .>= -1000
        x[1:m] .<= 1000
        x[m+1:end] .== fₑ₂
    end)
    # JuMP.@constraint(opti_structure, x[1:m₁].>=0)
    # JuMP.@constraint(opti_structure, x.<=10000)
    # JuMP.@constraint(opti_structure, x.>=-10000)

    # 定义目标函数
    JuMP.@objective(opti_model, Min, dot(cᵀ,x))
    # 求解
    JuMP.solve(opti_model)
    # 返回最优目标函数值，最优解（原问题），最优解（对偶问题）
    b=JuMP.getvalue(x)
    # println([x.id for x in assembly.structure.cables])
    # println([force/elm.A/elm.E for (force,elm) in zip(JuMP.getvalue(x)[1:m₁],assembly.structure.cables)])
    #
    # println([x.id for x in assembly.structure.links])
    # println([force/elm.A/elm.E for (force,elm) in zip(JuMP.getvalue(x)[m₁+1:m],assembly.structure.links)])
    # return JuMP.getobjectivevalue(opti_model), JuMP.getvalue(x)[1:m]#, JuMP.getdual(constr)

end
