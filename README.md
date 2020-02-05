# Mozi.jl

[![](https://travis-ci.org/zhuoju36/Mozi.jl.svg?branch=master)](https://travis-ci.org/zhuoju36/Mozi.jl)
[![codecov](https://codecov.io/gh/zhuoju36/Mozi.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/zhuoju36/Mozi.jl)
[![MIT Licence](https://badges.frapsoft.com/os/mit/mit.svg?v=103)](https://opensource.org/licenses/mit-license.php)

An FEM solver core for structural engineering
一个开源的结构工程有限元分析框架

# 快速入门

## 建立结构
st=Structure()
add_uniaxial_metal!(st,"steel",2e11,0.2,7849.0474)
add_general_section!(st,"frame",4.26e-3,3.301e-6,6.572e-5,9.651e-8,1e-3,1e-3,0,0)

N=6000
for i in 0:N
    add_node!(st,i,3i/1000,2i/1000,2i/1000)
end

for i in 0:N-1
    add_beam!(st,i,i,i+1,"steel","frame")
end

set_nodal_restraint!(st,0,true,true,true,true,true,true)

## 设置边界荷载
lcset=LoadCaseSet()
add_static_case!(lcset,"DL",1.1)


## 整装
assembly=assemble!(st,lcset,path=PATH)

## 分析求解
solve(assembly)

## 后处理
r=result_nodal_reaction(assembly,"DL",0)
r=result_nodal_displacement(assembly,"DL",N)
