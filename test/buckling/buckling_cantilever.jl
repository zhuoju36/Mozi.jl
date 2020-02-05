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

st=Structure()
add_uniaxial_metal!(st,"steel",2e11,0.2,7849.0474)
add_general_section!(st,"frame",4.26e-3,3.301e-6,6.572e-5,9.651e-8,1e-3,1e-3,0,0)

N=8
for i in 0:N
    add_node!(st,i,0,0,10i/N)
end

for i in 0:N-1
    add_beam!(st,i,i,i+1,"steel","frame")
end

set_nodal_restraint!(st,0,true,true,true,true,true,true)


lcset=LoadCaseSet()

add_static_case!(lcset,"DL",1.)

@time begin
assembly=assemble!(st,lcset,path=PATH)
end

@time begin
solve(assembly)
end

r=result_nodal_reaction(assembly,"DL",0)
# @test r≈[-2.83076e-10,-9.04954e-10,8926.15,53556.9,-80335.3,-3.29473e-8] rtol=1e-2
r=result_nodal_displacement(assembly,"DL",N)
# @test r≈[0.454025, 0.302683, -0.98385, -0.0336002, 0.0504003, -5.77725e-13] rtol=1e-2
println(r)
