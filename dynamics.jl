# Now we can discus dynamics.
#
# We will use an example of a driven cavity in which we might have some decay.
#
# This time we need the entire package, including solvers for equations of dynamics.

using QuantumOptics

# The system will be a single-mode cavity, to good approximation a harmonic oscillator

cutoff = 64
B = FockBasis(cutoff)

ω = 1.0
n = number(B)
H₀ = ω * n
α = 4.0*exp(im*π/3)
ψ₀ = coherentstate(B, α)
a = destroy(B)
p = a+a' # TODO position and momentum operator
q = im*(a'-a)

# # Unitary dynamics

ts = 0:0.01:3*2π
_, ψs = timeevolution.schroedinger(ts,ψ₀,H₀);

#

using CairoMakie

lines(ts, real.(expect.((n,), ψs)))
current_figure()

#

lines(ts, norm.(ψs))
current_figure()

#

f = Figure()
Axis(f[1,1], aspect=DataAspect())
lines!(real.(expect.((p,), ψs)), real.(expect.((q,), ψs)), color=ts)
current_figure()

# ## Time-dependent drives

function Hdynamic(t, psi)
    H₀ + 4*sin(40ω*t)*p
end

_, ψs = timeevolution.schroedinger_dynamic(ts,ψ₀,Hdynamic;)

f = Figure()
Axis(f[1,1], aspect=DataAspect())
lines!(real.(expect.((p,), ψs)), real.(expect.((q,), ψs)), color=ts)
current_figure()

# This is pretty darn slow though, as it recreates temporary objects for each step of the solver.

using BenchmarkTools
@benchmark _, ψs = timeevolution.schroedinger_dynamic(ts,ψ₀,Hdynamic;)

# ## "Lazy" ops

Hlazy = H₀ + TimeDependentSum(t->4*sin(40ω*t), p)
@benchmark _, ψs = timeevolution.schroedinger_dynamic(ts,ψ₀,Hlazy;)

# TODO more examples of lazy operators that are not just about time-dependence
# TODO example of multiple dispatch by creating your own Lazy operator

# # Accuracy and Precision

# Let's make a very "oscillatory" drive, to make the solver sweat a bit.

Hlazy = H₀ + TimeDependentSum(t->40*sin(4000ω*t), p)
ts = 0:1.0:10
_, ψs = timeevolution.schroedinger_dynamic(ts,ψ₀,Hlazy;)

# the result should be normalized

ψs[end] |> norm

# we can try a few different solvers

using OrdinaryDiffEqLowOrderRK: DP5
using OrdinaryDiffEqTsit5: Tsit5
using OrdinaryDiffEqVerner: Vern8

#

@time _, ψs = timeevolution.schroedinger_dynamic(ts,ψ₀,Hlazy;);
norm(ψs[end])

#

@time _, ψs = timeevolution.schroedinger_dynamic(ts,ψ₀,Hlazy; alg=Tsit5());
norm(ψs[end])

#

@time _, ψs = timeevolution.schroedinger_dynamic(ts,ψ₀,Hlazy; alg=Vern8());
norm(ψs[end])

#

@time _, ψs = timeevolution.schroedinger_dynamic(ts,ψ₀,Hlazy; alg=Vern8(), reltol=0, abstol=1e-12);
norm(ψs[end])

# Helpful discussion on choosing a good solver, with plently of links to more resources can be found on [the julia discourse forum](https://discourse.julialang.org/t/setting-abstol-and-reltol-when-solving-the-schrodinger-equation-with-ordinarydiffeq/125534).

# # Non-unitary dynamics with the Lindblad master equation
#

ts = 0:0.01:3*2π
_, ρs = timeevolution.master(ts,ψ₀,H₀,[√0.1 * a]);

#

f = Figure()
Axis(f[1,1], aspect=DataAspect())
lines!(real.(expect.((p,), ρs)), real.(expect.((q,), ρs)), color=ts)
current_figure()

#

# TODO unhelpful error message if J is not is not a vector (but rather just an operator) in _, ρs = timeevolution.master(ts,ψ₀,H₀,a)
# TODO the docs mark these as internal
# TODO Lazy time dependent ops with master_dynamic

# # Non-unitary dynamics with Monte Carlo Trajectories

ts = 0:0.01:3*2π
_, ψs = timeevolution.mcwf(ts,ψ₀,H₀,[√0.1 * a]);

#

f = Figure()
Axis(f[1,1], aspect=DataAspect())
lines!(real.(expect.((p,), ψs)), real.(expect.((q,), ψs)), color=ts)
current_figure()

# To really see something interesting with the MC trajectories we want non-Gaussian dynamics

ts = 0:0.1:3*2π
ψ₀ = fockstate(B,20)
_, ρs = timeevolution.master(ts,ψ₀,H₀,[√0.1 * a]);

lines(ts, real.(expect.((n,), ρs)))

#

sols = []
for i in 1:40
    _, ψs = timeevolution.mcwf(ts,ψ₀,H₀,[√0.1 * a]);
    push!(sols, ψs)
    lines!(ts, real.(expect.((n,), ψs)), color=(:gray,0.1))
end
current_figure()

#

sols_avg = sum([dm.(psis) for psis in sols])/length(sols)
lines!(ts, real.(expect.((n,), sols_avg)))
current_figure()

# ## Extracting more from a single trajectory
#
# There is an even more efficient way to think of this, as a way to have rigorous low bounds on performance.
#
# We will use a more realistic example, a bosonic code in which we perform some logical unitary.

α = 4.0
l0 = (coherentstate(B, α) + coherentstate(B, -α))/√2
l1 = (coherentstate(B, im*α) + coherentstate(B, -im*α))/√2;

#

quads = -10:0.1:10

f = Figure()
Axis(f[1,1], aspect=DataAspect())
w = wigner(l0, quads, quads)
m = maximum(abs.(w))
heatmap!(quads, quads, w, colormap=:redsblues, colorrange=(-m,m))
Axis(f[1,2], aspect=DataAspect())
w = wigner(l1, quads, quads)
m = maximum(abs.(w))
heatmap!(quads, quads, w, colormap=:redsblues, colorrange=(-m,m))
current_figure()

#

H = n
ts = 0:0.01:pi/2
_,ψs = timeevolution.schroedinger(ts,l0,H)

pl0 = projector(l0)
pl1 = projector(l1)

lines(ts, real.(expect.((pl0,),ψs)))
lines!(ts, real.(expect.((pl1,),ψs)))
current_figure()

# Now with loss

@time _,ρs = timeevolution.master(ts,l0,H,[√0.01 * a]);
lines!(ts, real.(expect.((pl0,),ρs)), linestyle=:dash)
lines!(ts, real.(expect.((pl1,),ρs)), linestyle=:dash)
current_figure()

# And finally, with a non-hermitian Hamiltonian

@time _,ψnhs = timeevolution.schroedinger(ts,l0,H - 0.01im * a * a');
lines!(ts, real.(expect.((pl0,),ψnhs)), linestyle=:dot)
lines!(ts, real.(expect.((pl1,),ψnhs)), linestyle=:dot)
current_figure()

# TODO Nakajima Zwanzig solver
