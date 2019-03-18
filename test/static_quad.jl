
include("../src/Mozi.jl")
using Test
using Logging

using .Mozi

const PATH=pwd()
macro showbanner(word,total=99)
    n=length(word)
    m=(total-n)÷2
    for i in 1:m
        print("-")
    end
    print(word)
    for i in 1:total-m-n
        print("-")
    end
    println()
end

@showbanner "Basic quad membrane test"
st=Structure()
lcset=LoadCaseSet()

add_uniaxial_metal!(st,"steel",2e11,0.3,7849.0474)

add_node!(st,1,1,1,0)
add_node!(st,2,-1,1,0)
add_node!(st,3,-1,-1,0)
add_node!(st,4,1,-1,0)

add_quad!(st,1,1,2,3,4,"steel",1e-3)

add_static_case!(lcset,"DL",0)
add_nodal_force!(lcset,"DL",4,0,-1e5,0,0,0,0)

set_nodal_restraint!(st,1,true,true,true,true,true,true)
set_nodal_restraint!(st,2,true,true,true,true,true,true)

assembly=assemble!(st,lcset,path=PATH)
solve(assembly)

r=result_nodal_displacement(assembly,"DL",4)
@show r

@showbanner "Basic quad bending test"
lcset=LoadCaseSet()

add_static_case!(lcset,"DL",0)
add_nodal_force!(lcset,"DL",4,0,0,-10,0,0,0)

assembly=assemble!(st,lcset,path=PATH)
solve(assembly)

r=result_nodal_displacement(assembly,"DL",4)
@test r≈[0.0, 0.0, -1.02971, 0.744574, 0.218848, 0.0] atol=1e-3

@showbanner "Isoparam quad membrane test"
st=Structure()
lcset=LoadCaseSet()

add_uniaxial_metal!(st,"steel",2e11,0.3,7849.0474)

add_node!(st,1,2,1,0)
add_node!(st,2,-1,1,0)
add_node!(st,3,-1,-0.5,0)
add_node!(st,4,2,-0.5,0)
add_quad!(st,1,4,3,2,1,"steel",1e-3)

set_nodal_restraint!(st,1,true,true,true,true,true,true)
set_nodal_restraint!(st,2,true,true,true,true,true,true)

add_static_case!(lcset,"DL",0)
add_nodal_force!(lcset,"DL",4,0,-1e5,0,0,0,0)
assembly=assemble!(st,lcset,path=PATH)
solve(assembly)

r=result_nodal_displacement(assembly,"DL",4)
@show r

@showbanner "Isoparam quad bending test"
lcset=LoadCaseSet()

add_static_case!(lcset,"DL",0)
add_nodal_force!(lcset,"DL",4,0,0,-10,0,0,0)

assembly=assemble!(st,lcset,path=PATH)
solve(assembly)

r=result_nodal_displacement(assembly,"DL",4)
@show r

@showbanner "Quad cantilever test"
st=Structure()
lcset=LoadCaseSet()
#
add_uniaxial_metal!(st,"steel",2e11,0.3,7849.0474)
add_node!(st,1,0,0,0)
add_node!(st,2,6,0,0)
add_node!(st,3,12,0,0)
add_node!(st,4,18,0,0)
add_node!(st,5,0,6,0)
add_node!(st,6,6,6,0)
add_node!(st,7,12,6,0)
add_node!(st,8,18,6,0)

add_quad!(st,1,1,2,6,5,"steel",1e-3)
add_quad!(st,2,2,3,7,6,"steel",1e-3)
add_quad!(st,3,3,4,8,7,"steel",1e-3)
add_static_case!(lcset,"DL",0)
add_nodal_force!(lcset,"DL",8,0,-1e5,0,0,0,0)

set_nodal_restraint!(st,1,true,true,true,true,true,true)
set_nodal_restraint!(st,5,true,true,true,true,true,true)

assembly=assemble!(st,lcset,path=PATH)
solve(assembly)

r=result_nodal_displacement(assembly,"DL",8)
@show r
