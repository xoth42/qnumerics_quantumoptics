using QuantumOptics

cutoff = 64
B = FockBasis(cutoff)
ω = 1.0
# number([T=ComplexF64,] b::FockBasis)

#   Number operator for the
#   specified Fock space
#   with optional data type
#   T.
n = number(B) 
# hubert space
H = ω * n

α = 4 # amplitude: exitation of average number 2
ψ=coherentstate(B,α)

# check if correct
destroy(B)*ψ ≈ α *ψ

# take a look at this...
using CairoMakie
# make a wigner diagram

#quadratre
quad = -10:0.1:10
w = wigner(ψ,quad,quad)
fig = Figure()
ax = Axis(fig[1,1], aspect=DataAspect())
heatmap!(quad,quad,w)
fig

# make it show the moving values
ts = 0:0.1:3*2π
_, ψs = timeevolution.schroedinger(ts,ψ,H);

fig = Figure()
ax = Axis(fig[1,1], aspect=DataAspect())
heatmap!(quad,quad,w)
w = wigner(ψs[end-2],quad,quad)

ax = Axis(fig[1,2], aspect=DataAspect())
heatmap!(quad,quad,w)
fig

##

a = destroy(B)
# position
x = a+ a'
# momentum
p = im*(a'-a)


fig = Figure()
ax = Axis(fig[1,1], aspect=DataAspect())
    # alternatively, x_expectation_vals = [expect(x,state) for state in ψs]
lines!(
    real.(expect.((x,),ψs)), # broadcast to get expected value for the psi's
    real.(expect.((p,),ψs)),
    color=ts
)
fig


##

lines(norm.(ψs)) # shows some error, not perfectly 1

##

# Time dependent-Hamiltonian
# (silly way)
function  Hdynamic(t,ψ)
    return H+4*sin(20*ω*t)*p
end


_, ψs = timeevolution.schroedinger_dynamic(ts,ψ,Hdynamic)


fig = Figure(size=(300,300))
ax = Axis(fig[1,1], aspect=DataAspect())
    # alternatively, x_expectation_vals = [expect(x,state) for state in ψs]
lines!(
    real.(expect.((x,),ψs)), # broadcast to get expected value for the psi's
    real.(expect.((p,),ψs)),
    color=ts
)
fig

##
# the gc usage is kind of bad for the Hdynamic
using BenchmarkTools
@benchmark _, ψs = timeevolution.schroedinger_dynamic(ts,ψ,Hdynamic)
# ^ this is really bad (really bad) 
# mem 146 Mib, 100ms
# we can use some julia tools to speed this up, just like what we had done with 'mul!'

# sum vectors, not matricies, 'lazy operators' -> to avoid temporary large matricies
Hlazy = H + TimeDependentSum(t -> 4*sin(20*ω*t),p) # this uses mul! under the hood
# like => function mul!(buffer, ::TimeDependentSum, input)
#   mul!(buffer,p,input)
#   buffer .= 4*sin(20*ω*t) .* buffer
# end
# also there is SciMLOperators

@benchmark _, ψs = timeevolution.schroedinger_dynamic(ts,ψ,Hlazy)
# 3.7 Mib, 36ms

# Some other solvers we can use
using OrdinaryDiffEqLowOrderRK: DP5
using OrdinaryDiffEqTsit5: Tsit5
using OrdinaryDiffEqVerner: Vern8

# 'make it annoying to the solver'
Hlazy = H + TimeDependentSum(t -> 80*sin(800*ω*t),p) 
ts = 0:0.5:3*3*2π

# norm(ψs[end]) should be one, only not one if there is error and it is drifting
@time _, ψs = timeevolution.schroedinger_dynamic(ts,ψ,Hlazy);
@time _, ψs = timeevolution.schroedinger_dynamic(ts,ψ,Hlazy);
norm(ψs[end])


@time _, ψs = timeevolution.schroedinger_dynamic(ts,ψ,Hlazy;alg=Tsit5());
@time _, ψs = timeevolution.schroedinger_dynamic(ts,ψ,Hlazy;alg=Tsit5());
norm(ψs[end])


@time _, ψs = timeevolution.schroedinger_dynamic(ts,ψ,Hlazy;alg=Vern8(),reltol=0,abstol=1e-9);
@time _, ψs = timeevolution.schroedinger_dynamic(ts,ψ,Hlazy;alg=Vern8(),reltol=0,abstol=1e-9);
norm(ψs[end])


# Lindbald?? master eqn? -> open system simulation, only need H of system and colapse operator (difficult to derive)
H 
ts = 0:0.01:3*2π
# colapse ops, some sort of decay?
# decay rate 
γ = 0.2
c = √γ *  a
# ρs because answers will be desntiy matrices now. Colapse operator gets squared in the linbald op. Can have multiple colapse operators
_,ρs = timeevolution.master(ts,ψ,H,[c])


fig = Figure(size=(400,400))
ax = Axis(fig[1,1], aspect=DataAspect())
    # alternatively, x_expectation_vals = [expect(x,state) for state in ψs]
lines!(
    real.(expect.((x,),ρs)), # broadcast to get expected value for the psi's
    real.(expect.((p,),ρs)),
    color=ts
)
fig
# the most computationally expensive way to make a spiral
# decay is expnential, rotation is constant

w = wigner(ρs[end-40],quad,quad)

ax = Axis(fig[1,2], aspect=DataAspect())
heatmap!(quad,quad,w)
fig


## 

# No master eqn, instead we use a trick to sample the kets from a trajectory
# The colapse operator goes to -> a probablistic symbol with a possible descrete error/effect 
