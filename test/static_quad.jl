@showbanner "Basic quad membrane test"
st=Structure()
lcset=LoadCaseSet()

add_uniaxial_metal!(st,"steel",2e11,0.3,7849.0474)

add_node!(st,1,0,0,0)
add_node!(st,2,6,0,0)
add_node!(st,3,12,0,0)
add_node!(st,4,18,0,0)
add_node!(st,5,0,6,0)
add_node!(st,6,6,6,0)
add_node!(st,7,12,6,0)
add_node!(st,8,18,6,0)

add_quad!(st,1,1,2,6,5,"steel",1e-3;membrane=true,plate=false)
add_quad!(st,2,2,3,7,6,"steel",1e-3;membrane=true,plate=false)
add_quad!(st,3,3,4,8,7,"steel",1e-3;membrane=true,plate=false)

add_static_case!(lcset,"DL",0)
add_nodal_force!(lcset,"DL",8,0,-1e5,0,0,0,0)

set_nodal_restraint!(st,1,true,true,true,true,true,true)
set_nodal_restraint!(st,5,true,true,true,true,true,true)
for i in [2,3,4,6,7,8]
    set_nodal_restraint!(st,i,false,false,true,true,true,true)
end

assembly=assemble!(st,lcset,path=PATH)

solve(assembly)

r=result_nodal_displacement(assembly,"DL",8)
@show r
# passed=@test r≈[-2.83076e-10,-9.04954e-10,8926.15,53556.9,-80335.3,-3.29473e-8] atol=1e-1
# @show passed
