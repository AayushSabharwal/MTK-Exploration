using ModelingToolkit
using DifferentialEquations
using LinearAlgebra

@parameters t

D = Differential(t)

const g = [0., -9.81] # gravity

function Mass(; name, m = 1.0, xy = [0., 0.], uv = [0., 0.])
    ps = @parameters m=m    # every Mass has a mass...
    sts = @variables pos[1:2](t)=xy v[1:2](t)=uv # position and velocity
    eqs = collect(D.(pos) .~ v) # velocity is rate of change of position
    ODESystem(eqs, t, [pos..., v...], ps; name)
end

function Spring(; name, k = 1e4, l = 1.)
    ps = @parameters k=k l=l    # spring constant, and spring length
    @variables x(t), dir[1:2](t)    # spring extension, and direction
    ODESystem(Equation[], t, [x, dir...], ps; name)
end

# connect a spring between two positions
function connect_spring(spring, a, b)
    [
        spring.x ~ norm(collect(a .- b))
        collect(spring.dir .~ collect(a .- b))
    ]
end

# The (vector) force a spring exerts
spring_force(spring) = -spring.k .* collect(spring.dir) .* (spring.x - spring.l)  ./ spring.x

# masses and springs are named tuples
# connections are tuples (spring_index, mass_index, mass_index/point)
function get_system(masses, springs, connections)
    nmass = length(masses)
    mass = map(i -> Mass(; masses[i]..., name = Symbol(:mass_, i)), 1:nmass)
    nspr = length(springs)
    spr = map(i -> Spring(; springs[i]..., name = Symbol(:spr_, i)), 1:nspr)
    eqs = Equation[]
    forces = fill(zeros(Num, 2), nmass)
    for conn in connections
        forces[conn[2]] .+= spring_force(spr[conn[1]])
        if conn[3] isa Int
            forces[conn[3]] .-= spring_force(spr[conn[1]])
        end
    end
    for i in 1:nmass
        forces[i] = forces[i] / mass[i].m + g
        push!(eqs, collect(D.(mass[i].v) .~ forces[i])...)
    end
    for conn in connections
        if conn[3] isa Int
            push!(eqs, connect_spring(spr[conn[1]], mass[conn[2]].pos, mass[conn[3]].pos)...)
        else
            push!(eqs, connect_spring(spr[conn[1]], mass[conn[2]].pos, conn[3])...)
        end
    end
    return mass, spr, eqs
end

# use the equations and subsystems to create and return a solution and the overall system
function get_sol(eqs, subs...; tspan = (0., 15.))
    @named model = compose(ODESystem(eqs, t; name = :inner), subs...)
    sys = structural_simplify(model)
    prob = ODEProblem(sys, [], tspan)
    solve(prob, Rosenbrock23()), sys
end
