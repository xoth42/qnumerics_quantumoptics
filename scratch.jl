using QuantumOpticsBase

# Define a basis for the harmoic oscillator
# max dimention we will deal with
cutoff = 25
F = FockBasis(cutoff)

# Basis for the spin
S = SpinBasis(1//2)

# Number of photons for the cavity ?
# introduce some states
ψ = fockstate(F,5) # 5 excitations, 5 photons, less than the cutoff of the fockbasis 

s = spindown(S)

# more generally:
G = GenericBasis(100)
basisstate(G,3)

# to figure out how to make a state in some basis
tensor(F,S) |> typeof

mystate = ψ ⊗ s # a new bigger state 

# how big?
dm(mystate)

# it's very sparse, so let's use a sparce Array
sparse(dm(mystate)) 
sparse(dm(mystate))  |> typeof


# other operators
projector(ψ) ⊗ projector(s)
projector(ψ ⊗ s)

projector(ψ)⊗ identityoperator(S)

create(F) == destroy(F)' # ' -> adjoint()

sigmap(S) * spindown(S) == spinup(S)

# why sparse vs. dense ops?
using BenchmarkTools
B = GenericBasis(1000)
ψ2 = basisstate(B,101)
P = projector(ψ2)  
typeof(P)

@benchmark P*ψ2 # 250 μs

# now try sparse
Ps = sparse(P)
@benchmark Ps*ψ2 # 5 μs


# but let's try not to allocate memory on the fly, let's pre-allocate
using LinearAlgebra: mul!
buffer = copy(ψ2)
@benchmark mul!(buffer, Ps, ψ2) # 2 μs, no memory!!


ψ = fockstate(F,4) ⊗ spinup(S)

# what we were doing before to sim statevectors
myop = create(F) ⊗ identityoperator(S) 

# instead, use embed:embed(basis1[, basis2], indices::Vector, operators::Vector)

#   Tensor product of operators
#   where missing indices are filled
#   up with identity operators.

# create(F) lives in the hubert space of F, F⊗S is the larger hilbert space
embed(F⊗S,1,create(F))  == create(F) ⊗ identityoperator(S)


# bell pair
l = spindown(S)
h = spinup(S)
bell = (l⊗l + h⊗h) / √2

# density matrix of bell pair
dm(bell)

# do a partial trace, and measure the first subsystem
ptrace(dm(bell),1)