module LoadCase

include("./Load.jl")

using SparseArrays

using ..Enums
using .Load

export LoadCaseSet,AbstractLoadCase

# export StaticCase,ModalCase,BucklingCase,TimeHistoryCase,ResponseSpectrumCase

export add_static_case!, add_modal_case!, add_buckling_case!, add_time_history_case!

export add_nodal_disp!,add_nodal_force!,
add_cable_distributed!,add_cable_strain!,add_cable_stress!,
add_link_distributed!,add_link_strain!,add_link_stress!,
add_beam_distributed!,add_beam_strain!
export set_modal_params!,set_buckling_params!,set_time_history_params!
export is_static, is_modal, is_buckling, is_time_history, is_response_spectrum, get_all_id

abstract type AbstractLoadCase end

mutable struct StaticCase <: AbstractLoadCase
    id::String
    hid::Int

    nl_type::String
    plc::String

    gfactor::Float64
    nodal_forces::Dict{String,NodalForce}
    nodal_disps::Dict{String,NodalDisp}
    beam_forces::Dict{String,BeamForce}
    quad_forces::Dict{String,QuadForce}
    tria_forces::Dict{String,TriaForce}

    P::Vector{Float64}
    P̄::Vector{Float64}

    StaticCase(id,hid,nl_type="1st",plc="",gfactor=0.)=new(
        string(id),
        hid,
        nl_type,
        plc,
        gfactor,
        Dict{String,NodalForce}(),
        Dict{String,NodalDisp}(),
        Dict{String,BeamForce}(),
        Dict{String,QuadForce}(),
        Dict{String,TriaForce}(),
        zeros(1),
        zeros(1))
end

mutable struct ModalCase <: AbstractLoadCase
    id::String
    hid::Int

    modal_type::String
    plc::String

    nev::Int
    tolev::Number
    maxiterev::Int

    P::Vector
    P̄::Vector

    ModalCase(id,hid,modal_type,plc="")=new(
        string(id),
        hid,
        modal_type,
        plc,
        6,
        1e-9,
        100,
        zeros(1),
        zeros(1))
end

mutable struct BucklingCase <: AbstractLoadCase
    id::String
    hid::Int

    plc::String

    gfactor::Float64
    nodal_forces::Dict{String,NodalForce}
    nodal_disps::Dict{String,NodalDisp}
    beam_forces::Dict{String,BeamForce}
    quad_forces::Dict{String,QuadForce}

    shift::Float64
    positive_only::Bool

    nev::Int
    tolev::Number
    maxiterev::Int

    P::Vector
    P̄::Vector

    BucklingCase(id,hid,nl_type="1st",plc="",gfactor=0.)=new(
        string(id),
        hid,
        plc,
        gfactor,
        Dict{String,NodalForce}(),
        Dict{String,NodalDisp}(),
        Dict{String,BeamForce}(),
        Dict{String,QuadForce}(),
        0.,
        false,
        6,
        1e-12,
        100,
        zeros(1),
        zeros(1))
end

mutable struct TimeHistoryCase <: AbstractLoadCase
    id::String
    hid::Int

    t::Array{Float64}
    f::Array{Float64}

    plc::String

    nodal_forces::Dict{String,NodalForce}
    nodal_disps::Dict{String,NodalDisp}
    beam_forces::Dict{String,BeamForce}
    quad_forces::Dict{String,QuadForce}

    algorithm::String
    α::Float64
    β::Float64
    γ::Float64
    θ::Float64
    modal_case::String

    P::Vector
    P̄::Vector
    TimeHistoryCase(id,hid,t,f,plc="")=new(
        string(id),
        hid,
        t,f,
        plc,
        Dict{String,NodalForce}(),
        Dict{String,NodalDisp}(),
        Dict{String,BeamForce}(),
        Dict{String,QuadForce}(),
        "newmark",
        0.85,
        0.25,
        0.5,
        1.4,
        "",
        zeros(1),
        zeros(1))
end

mutable struct ResponseSpectrumCase <: AbstractLoadCase
    id::String
    hid::Int

    t::Array{Float64}
    α::Array{Float64}

    modal_case::String

    plc::String

    dirfactors::Array{Float64}
end

mutable struct LoadCaseSet
    statics::Dict{String,StaticCase}
    modals::Dict{String,ModalCase}
    bucklings::Dict{String,BucklingCase}
    time_histories::Dict{String,TimeHistoryCase}
    response_spectrums::Dict{String,ResponseSpectrumCase}

    lc_tree::Array
    LoadCaseSet()=new(Dict{String,StaticCase}(),
    Dict{String,ModalCase}(),
    Dict{String,BucklingCase}(),
    Dict{String,TimeHistoryCase}(),
    Dict{String,ResponseSpectrumCase}(),
    [])
end

get_all_id(lcset::LoadCaseSet)=union(keys(lcset.statics),keys(lcset.modals),keys(lcset.bucklings),keys(lcset.time_histories),keys(lcset.response_spectrums))
is_static(loadcase)=loadcase isa StaticCase
is_modal(loadcase)=loadcase isa ModalCase
is_buckling(loadcase)=loadcase isa BucklingCase
is_time_history(loadcase)=loadcase isa TimeHistoryCase
is_response_spectrum(loadcase)=loadcase isa ResponseSpectrumCase

include("./interface.jl")

end
