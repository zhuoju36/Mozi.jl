include("../src/Mozi.jl")

using Test
using Logging

using .Mozi

const PATH=pwd()

@info "---------- P-Δ test ----------"
st=Structure()
lcset=LoadCaseSet()

add_uniaxial_metal!(st,"steel",2e11,0.2,7849.0474)
add_uniaxial_metal!(st,"cable",1.9e11,0.3,7849.0474)

add_beam_section!(st,"frame",3,1500,3000,20000,30000)
add_beam_section!(st,"string",5,60)

add_node!(st,1,0,0,0)
add_node!(st,2,10,0,0)
add_node!(st,3,10,10,0)
add_node!(st,4,20,10,0)
add_node!(st,5,10,0,-5)
add_node!(st,6,10,10,-5)

add_beam!(st,"b1",1,2,"steel","frame")
add_beam!(st,"b2",2,3,"steel","frame")
add_beam!(st,"b3",3,4,"steel","frame")
add_beam!(st,"c1",1,5,"cable","string")
add_beam!(st,"c2",5,6,"cable","string")
add_beam!(st,"c3",6,4,"cable","string")
set_beam_release!(st,"c1",false,false,false,false,true,true,false,false,false,false,true,true)

set_nodal_restraint!(st,1,true,true,true,true,true,true)
set_nodal_restraint!(st,6,true,true,true,true,true,true)

add_static_case!(lcset,"D",1)
add_static_case!(lcset,"D>P",0,nl_type="2nd",plc="D")
add_static_case!(lcset,"D>P>L",0,nl_type="2nd",plc="D>P")


add_beam_distributed!(lcset,"D>P>L","b1",0,0,3000*5,0,0,0,0,0,3000*5,0,0,0)
add_beam_distributed!(lcset,"D>P>L","b2",0,0,3000*5,0,0,0,0,0,3000*5,0,0,0)
add_beam_distributed!(lcset,"D>P>L","b3",0,0,3000*5,0,0,0,0,0,3000*5,0,0,0)

assembly=assemble!(st,lcset,path=PATH)
solve(assembly)

r=result_nodal_displacement(assembly,"D>P>L",2)
@show r
