using ModelingToolkit
using DifferentialEquations
using LinearAlgebra

@parameters t

D = Differential(t)

const g = [0., -9.81] # gravity


@connector function Mass(; name, m = 1.0, xy = [0., 0.], uv = [0., 0.])
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

# connect a spring between a mass and a point
function connect_spring(spring, mass, point::Vector{Float64})
    [
        spring.x ~ norm(collect(mass.pos .- point)) # extension is the distance
        collect(spring.dir .~ collect(mass.pos .- point))   # direction is the vector
    ]
end

# connect a spring between two masses
function connect_spring(spring, a, b)
    [
        spring.x ~ norm(collect(a.pos .- b.pos))
        collect(spring.dir .~ collect(a.pos .- b.pos))
    ]
end

# The (vector) force a spring exerts
spring_force(spring) = -spring.k .* collect(spring.dir) .* (spring.x - spring.l)  ./ spring.x

# use the equations and subsystems to create and return a solution and the overall system
function get_sol(eqs, subs...; tspan = (0., 15.))
    @named model = compose(ODESystem(eqs, t; name = :inner), subs...)
    sys = structural_simplify(model)
    prob = ODEProblem(sys, [], tspan)
    solve(prob, Rosenbrock23()), sys
end
