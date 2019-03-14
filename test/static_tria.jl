@showbanner "Basic tria membrane test"
st=Structure()
lcset=LoadCaseSet()

add_uniaxial_metal!(st,"steel",2e5,0.2,7849.0474)
add_static_case!(lcset,"DL",0)

add_node!(st,1,0,0,0)
add_node!(st,2,6,0,0)
add_node!(st,3,12,0,0)
add_node!(st,4,18,0,0)
add_node!(st,5,0,6,0)
add_node!(st,6,6,6,0)
add_node!(st,7,12,6,0)
add_node!(st,8,18,6,0)

add_tria!(st,1,1,6,5,"steel",1e-3)
add_tria!(st,2,1,2,6,"steel",1e-3)
add_tria!(st,3,2,7,6,"steel",1e-3)
add_tria!(st,4,2,3,7,"steel",1e-3)
add_tria!(st,5,3,8,7,"steel",1e-3)
add_tria!(st,6,3,4,8,"steel",1e-3)

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
