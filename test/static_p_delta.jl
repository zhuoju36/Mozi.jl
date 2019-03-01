@info "---------- P-Δ test ----------"
st=Structure()
lcset=LoadCaseSet()

add_uniaxial_metal!(st,"steel",2e11,0.2,7849.0474)
add_general_section!(st,"frame",4.26e-3,3.301e-6,6.572e-5,9.651e-8,1e-3,1e-3,0,0)

add_node!(st,1,0,0,0)
add_node!(st,2,18,12,12)
add_beam!(st,1,1,2,"steel","frame")

set_nodal_restraint!(st,1,true,true,true,true,true,true)

add_static_case!(lcset,"DL",0.5)
add_static_case!(lcset,"DL+SD",0.3,nl_type="2nd",plc="DL")
add_static_case!(lcset,"DL+SD+LL",0.3,nl_type="2nd",plc="DL+SD")

assembly=assemble!(st,lcset,path=PATH)
solve(assembly)

r=result_nodal_displacement(assembly,"DL+SD+LL",2)
passed=@test r≈[0.454025, 0.302683, -0.98385, -0.0336002, 0.0504003, -5.94366e-13] atol=1e-3
@show passed

r=result_nodal_reaction(assembly,"DL+SD+LL",1)
passed=@test r≈[-2.83076e-10,-9.04954e-10,8926.15,53556.9,-80335.3,-3.29473e-8] atol=10
@show passed
