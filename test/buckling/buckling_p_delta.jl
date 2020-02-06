include("../../src/Mozi.jl")

using .Mozi
using Test

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

@info "---------- P-Δ test ----------"
st=Structure()
lcset=LoadCaseSet()

add_uniaxial_metal!(st,"steel",2e11,0.2,7849.0474)
add_uniaxial_metal!(st,"cable",1.9e11,0.3,7849.0474)

add_beam_section!(st,"frame",3,1500,3000,20,30)
add_beam_section!(st,"string",5,60)

add_node!(st,1,0,0,0)
add_node!(st,2,0,0,1)
add_node!(st,3,1,0,0)
add_node!(st,4,0,1,0)
add_node!(st,5,-1,0,0)
add_node!(st,6,0,-1,0)

add_beam!(st,"b1",1,2,"steel","frame")
add_beam!(st,"c1",3,2,"cable","string")
add_beam!(st,"c2",4,2,"cable","string")
add_beam!(st,"c3",5,2,"cable","string")
add_beam!(st,"c4",6,2,"cable","string")

# set_beam_release!(st,"c1",false,false,false,false,true,true,false,false,false,false,true,true)
# set_beam_release!(st,"c2",false,false,false,false,true,true,false,false,false,false,true,true)
# set_beam_release!(st,"c3",false,false,false,false,true,true,false,false,false,false,true,true)
# set_beam_release!(st,"c4",false,false,false,false,true,true,false,false,false,false,true,true)

set_nodal_restraint!(st,1,true,true,true,true,true,true)
set_nodal_restraint!(st,3,true,true,true,true,true,true)
set_nodal_restraint!(st,4,true,true,true,true,true,true)
set_nodal_restraint!(st,5,true,true,true,true,true,true)
set_nodal_restraint!(st,6,true,true,true,true,true,true)


add_static_case!(lcset,"D",0)
add_static_case!(lcset,"D>P",0,nl_type="2nd",plc="D")
add_static_case!(lcset,"D>P>L",0,nl_type="2nd",plc="D>P")

add_beam_strain!(lcset,"D","c1",-0.03)
add_beam_strain!(lcset,"D","c2",-0.03)
add_beam_strain!(lcset,"D","c3",-0.03)
add_beam_strain!(lcset,"D","c4",-0.03)

# add_nodal_force!(lcset,"D",2,0,0,-1e6,0,0,0)

assembly=assemble!(st,lcset,path=PATH)
solve(assembly)

r=result_nodal_displacement(assembly,"D",2)
# r=result_nodal_displacement(assembly,"D>P",2)
# r=result_nodal_displacement(assembly,"D>P>L",2)
println(r)
