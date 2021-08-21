using ModelingToolkit
using DifferentialEquations
using LinearAlgebra

@parameters t

const D = Differential(t)

struct Mass
    initial_pos::Vector{Float64}
    system::ODESystem
    forces::Vector{Num}
end

function Mass(; name, m = 1.0, ipos = [0., 0.], ivel = [0., 0.])
    ps = @parameters m=m
    @variables pos[1:2](t)=ipos vel[1:2](t)=ivel
    eqs = collect(D.(pos) .~ vel)
    sys = ODESystem(eqs, t, [pos..., vel...], ps; name)
    return Mass(
        ipos,
        sys,
        Num[0, 0],
    )
end

function Base.getproperty(mass::Mass, prop::Symbol)
    if prop == :initial_pos
        return getfield(mass, prop)
    elseif prop == :system
        return getfield(mass, prop)
    elseif prop == :forces
        return getfield(mass, prop)
    else
        return getproperty(mass.system, prop)
    end
end

struct Spring
    system::ODESystem
    equations::Vector{Equation}
    force::Vector{Num}
end

function Spring(; name, k = 1e4, l = 1.0)
    ps = @parameters k=k l=l
    @variables x(t), dir[1:2](t)
    sys = ODESystem(Equation[], t, [x, dir...], ps; name)
    return Spring(
        sys,
        Equation[],
        -sys.k .* collect(sys.dir) .* (sys.x - sys.l)  ./ sys.x,
    )
end

function Base.getproperty(spr::Spring, prop::Symbol)
    if prop == :system
        return getfield(spr, prop)
    elseif prop == :equations
        return getfield(spr, prop)
    elseif prop == :force
        return getfield(spr, prop)
    else
        return getproperty(spr.system, prop)
    end
end

function (spr::Spring)(a_pos, b_pos)
    push!(spr.equations, [
        spr.system.x ~ norm(collect(a_pos .- b_pos))
        collect(spr.system.dir .~ collect(a_pos .- b_pos))
    ]...)
end

function (spr::Spring)(a::Mass, b::Mass)
    spr(a.system.pos, b.system.pos)
    a.forces .+= spr.force
    b.forces .-= spr.force
end

function (spr::Spring)(a::Mass, b::Vector{Float64})
    spr(a.system.pos, b)
    a.forces .+= spr.force
end

function create_system(masses::Vector{Mass}, springs::Vector{Spring}; name, g = [0., -9.81])
    eqs = Equation[]
    for spr in springs
        push!(eqs, spr.equations...)
    end

    for m in masses
        push!(eqs, collect(D.(m.system.vel) .~ m.forces ./ m.system.m .+ g)...)
    end
    
    model = compose(ODESystem(eqs, t; name = Symbol(:_, name)), (m.system for m in masses)..., (s.system for s in springs)...; name)
    structural_simplify(model)
end
