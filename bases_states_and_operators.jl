]#=

Instead simply using Array types or sparse matrices,
we can use one of the many available state vector packages,
which wraps these types in extra metadata, simplifying the
tracking of subsystems of varying dimensions in the larger system.

We will use QuantumOptics because of its simplicity,
and reasonable completeness. It happens to be relatively fast too,
and makes it easy to create weird custom operators.

Other equally good options are:

- qutip in python (everyone knows it, fueled the whole field)
- QuantumToolbox in julia (close clone of qutip in Julia)
- dinamiqs in jax (best in autodiff)

=#

# We will start with a smaller library
# that helps create convenient datastructures and tracks metadata for us,
# but does no provide any solvers for dynamics or other advanced numerics.

using QuantumOpticsBase

# # Bases and States

# make a bosonic basis with a cutoff

ℱ = FockBasis(10)

# a state of five excitations

f₅ = fockstate(ℱ, 5)

# make spin basis for a two-level system

𝒮 = SpinBasis(1//2)

# and a state

sᵤ = spinup(𝒮)

#

sᵤ |> typeof

#

sᵤ |> typeof |> supertype

# # Tensor Products

Ψ = tensor(f₅, sᵤ)

# shorthand unicode notation

Ψ = f₅ ⊗ sᵤ

# has a composite basis

basis(Ψ)

# # Operators

# we can convert a ket into a density matrix

ρ = dm(f₅) ⊗ dm(sᵤ)

# have it dense or sparse

ρ = sparse(ρ)

# ## adjoints, projectors, partial traces

projector(sᵤ)

#

projector(f₅) |> sparse

#

ρ ≈ projector(f₅) ⊗ projector(sᵤ)

#

adjoint(sᵤ) * sᵤ

#

sᵤ' * sᵤ

# ## named operators

identityoperator(ℱ)

#

number(ℱ)

#

sigmax(𝒮)

#
# ## Sparsity helps with large systems, but can be detrimental with small ones

using BenchmarkTools

B = GenericBasis(100)
ψ = basisstate(B,10)+basisstate(B,1)
ρ = dense(dm(ψ))

@benchmark ρ*ψ

#

ρ = sparse(dm(ψ))

@benchmark ρ*ψ

#

B = GenericBasis(3)
ψ = basisstate(B,1)+basisstate(B,2)+basisstate(B,3)
ρ = dense(dm(ψ))

@benchmark ρ*ψ

#

ρ = sparse(dm(ψ))

@benchmark ρ*ψ

# ## In-place operations

using LinearAlgebra: mul!

ψₒ = copy(ψ)
ρ = dense(ρ)
@benchmark mul!(ψₒ, ρ, ψ)

#

ψₒ = copy(ψ)
ρ = sparse(ρ)
@benchmark mul!(ψₒ, ρ, ψ)

# ## Embedding

embed(ℱ⊗𝒮, 1, number(ℱ)) ≈ number(ℱ)⊗identityoperator(𝒮)

# ## Partial Traces

ψ₁ = fockstate(ℱ, 1)
ψ₂ = fockstate(ℱ, 2)
su = spinup(𝒮)
sd = spindown(𝒮)

Ψ = (ψ₁⊗su + ψ₂⊗sd) / √2
ρ = dm(Ψ)

purity(x) = tr(x^2)
purity(ρ)

#

ρₜ = ptrace(ρ, 2)
tr(ρₜ), purity(ρₜ)

#

basis(ρ), basis(ρₜ)
