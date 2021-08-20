using ModelingToolkit, DifferentialEquations

@parameters t

# a point in the circuit with a given voltage and current
@connector function Pin(; name)
    sts = @variables v(t)=1.0 i(t)=1.0
    ODESystem(Equation[], t, sts, []; name)
end

# ground pin, where voltage is 0
function Ground(; name)
    @named g = Pin()
    eqs = [g.v ~ 0.]
    compose(ODESystem(eqs, t, [], []; name), g)
end

# abstraction for all 2-pin components
function OnePort(; name)
    # two pins
    @named p = Pin()
    @named n = Pin()
    # component has a voltage across it, and current through
    sts = @variables v(t)=1.0 i(t)=1.0
    eqs = [
        v ~ p.v - n.v    # KVL
        0. ~ p.i + n.i   # KCL
        i ~ p.i          # Current through component is current through +ve pin
    ]
    compose(ODESystem(eqs, t, sts, []; name), p, n)
end

function Resistor(; name, R = 1.0)
    # 2 pin component
    @named oneport = OnePort()
    @unpack v, i = oneport
    ps = @parameters R=R
    eqs = [
        v ~ i * R
    ]
    extend(ODESystem(eqs, t, [], ps; name), oneport)
end

function Capacitor(; name, C = 1.0)
    @named oneport = OnePort()
    @unpack v, i = oneport
    ps = @parameters C=C
    D = Differential(t)
    eqs = [
        D(v) ~ i / C
    ]
    extend(ODESystem(eqs, t, [], ps; name), oneport)
end

function Inductor(; name, L = 1.0)
    @named oneport = OnePort()
    @unpack v, i = oneport
    ps = @parameters L=L
    D = Differential(t)
    eqs = [
        D(i) ~ v / L
    ]
    extend(ODESystem(eqs, t, [], ps; name), oneport)
end

function DCVoltageSource(; name, V = 1.0)
    @named oneport = OnePort()
    @unpack v = oneport
    ps = @parameters V=V
    eqs = [
        V ~ v
    ]
    extend(ODESystem(eqs, t, [], ps; name), oneport)
end

# V = Vm sin(ωt+ϕ)
function ACVoltageSource(; name, Vm = 1.0, ω = 2π, ϕ = 0)
    @named oneport = OnePort()
    @unpack v, i = oneport
    ps = @parameters Vm=Vm ω=ω ϕ=ϕ
    eqs = [
        v ~ Vm * sin(ω * t + ϕ)
    ]
    extend(ODESystem(eqs, t, [], ps; name), oneport)
end

star_connect(comps...) = connect([comp.n for comp in comps]...)

function ModelingToolkit.connect(::Type{Pin}, ps...)
    eqs = [
        0. ~ sum(p->p.i, ps)    # KCL
    ]
    # KVL
    for i in 1:length(ps)-1
        push!(eqs, ps[i].v ~ ps[i+1].v)
    end
    return eqs
end

function series_connection(comps...)
    eqs = Equation[]
    for i in 1:length(comps)-1
        eqs = vcat(eqs, connect(comps[i].p, comps[i+1].n))
    end
    return eqs
end

function parallel_connection(comps...)
    return [
        connect((comp.p for comp in comps)...)
        connect((comp.n for comp in comps)...)
    ]
end
