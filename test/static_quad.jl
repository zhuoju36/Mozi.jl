
include("../src/Mozi.jl")
using Test
using Logging

using .Mozi

const PATH=pwd()

# @showbanner "Basic quad membrane test"
st=Structure()
lcset=LoadCaseSet()

add_uniaxial_metal!(st,"steel",2e11,0.3,7849.0474)

# add_node!(st,1,0,0,0)
# add_node!(st,2,6,0,0)
# add_node!(st,3,12,0,0)
# add_node!(st,4,18,0,0)
# add_node!(st,5,0,6,0)
# add_node!(st,6,6,6,0)
# add_node!(st,7,12,6,0)
# add_node!(st,8,18,6,0)
#
# add_quad!(st,1,1,2,6,5,"steel",1e-3)
# add_quad!(st,2,2,3,7,6,"steel",1e-3)
# add_quad!(st,3,3,4,8,7,"steel",1e-3)
# add_static_case!(lcset,"DL",0)
# add_nodal_force!(lcset,"DL",8,0,-1e5,0,0,0,0)
#
# set_nodal_restraint!(st,1,true,true,true,true,true,true)
# set_nodal_restraint!(st,5,true,true,true,true,true,true)
#
# for i in [2,3,4,6,7,8]
#     set_nodal_restraint!(st,i,false,false,true,true,true,false)
# end
#
# assembly=assemble!(st,lcset,path=PATH)
#
# solve(assembly)
#
# r=result_nodal_displacement(assembly,"DL",8)
# @show r

add_node!(st,1,1,1,0)
add_node!(st,2,-1,1,0)
add_node!(st,3,-1,-1,0)
add_node!(st,4,1,-1,0)

add_quad!(st,1,1,2,3,4,"steel",1e-3)

add_static_case!(lcset,"DL",0)
add_nodal_force!(lcset,"DL",4,0,-1e5,0,0,0,0)

set_nodal_restraint!(st,1,true,true,true,true,true,true)
set_nodal_restraint!(st,2,true,true,true,true,true,true)

for i in [3,4]
    set_nodal_restraint!(st,i,false,false,true,true,true,false)
end

assembly=assemble!(st,lcset,path=PATH)
Kᵉ=assembly.structure.quads["1"].Kᵉ
T=assembly.structure.quads["1"].T


# solve(assembly)
#
# r=result_nodal_displacement(assembly,"DL",4)
#
# @show r

# @showbanner Basic quad bending test
st=Structure()
lcset=LoadCaseSet()
#
add_uniaxial_metal!(st,"steel",2e11,0.3,7849.0474)

add_node!(st,1,1,1,0)
add_node!(st,2,-1,1,0)
add_node!(st,3,-1,-1,0)
add_node!(st,4,1,-1,0)

add_quad!(st,1,1,2,3,4,"steel",1e-3)

add_static_case!(lcset,"DL",0)
add_nodal_force!(lcset,"DL",4,0,0,-10,0,0,0)

set_nodal_restraint!(st,1,true,true,true,true,true,true)
set_nodal_restraint!(st,2,true,true,true,true,true,true)


for i in [3,4]
    set_nodal_restraint!(st,i,true,true,false,false,false,true)
end

assembly=assemble!(st,lcset,path=PATH)

solve(assembly)

r=result_nodal_displacement(assembly,"DL",4)
@show assembly.structure.quads["1"].Kᵉ
@show r
@test r≈[0,0,-1.0643,0.79525,0.22807,0] atol=1e-3

# st=Structure()
# lcset=LoadCaseSet()
# #
# add_uniaxial_metal!(st,"steel",2e11,0.3,7849.0474)
# add_node!(st,1,0,0,0)
# add_node!(st,2,6,0,0)
# add_node!(st,3,12,0,0)
# add_node!(st,4,18,0,0)
# add_node!(st,5,0,6,0)
# add_node!(st,6,6,6,0)
# add_node!(st,7,12,6,0)
# add_node!(st,8,18,6,0)
#
# add_quad!(st,1,1,2,6,5,"steel",1e-3;membrane=true,plate=false)
# add_quad!(st,2,2,3,7,6,"steel",1e-3;membrane=true,plate=false)
# add_quad!(st,3,3,4,8,7,"steel",1e-3;membrane=true,plate=false)
# add_static_case!(lcset,"DL",0)
# add_nodal_force!(lcset,"DL",8,0,-1e5,0,0,0,0)
#
# set_nodal_restraint!(st,1,true,true,true,true,true,true)
# set_nodal_restraint!(st,5,true,true,true,true,true,true)
