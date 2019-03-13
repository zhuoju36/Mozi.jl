st=Structure()
lcset=LoadCaseSet()

add_uniaxial_metal!(st,"steel",2e5,0.2,7849.0474)

# add_node!(st,1,0,0,0)
# add_node!(st,2,6,0,0)
# add_node!(st,3,12,0,0)
# add_node!(st,4,18,0,0)
# add_node!(st,5,0,6,0)
# add_node!(st,6,6,6,0)
# add_node!(st,7,12,6,0)
# add_node!(st,8,18,6,0)
add_node!(st,1,0,0,0)
add_node!(st,2,3,0,0)
add_node!(st,3,0,4,0)

add_tria!(st,1,1,2,3,"steel",1)

assembly=assemble!(st,lcset,path=pwd())

# solve(assembly)
