module Glauber

export generate_random_conf, energy, magnetization, glauber_alg

"""
    generate_random_conf(N)

Return a random `Vector{Int}` of length `N` with entries in `{-1, +1}`.
"""
function generate_random_conf(N::Integer)
    return rand((-1, 1), N)
end

"""
    energy(J, σ; h=0.0)

Total energy of the 1D Ising chain `σ` with periodic boundary conditions,
coupling `J`, and external field `h`.

`E = -J Σ σ_i σ_{i+1} - h Σ σ_i`
"""
function energy(J::Real, σ::AbstractVector{<:Integer}; h::Real = 0.0)
    N = length(σ)
    E = -J * σ[N] * σ[1] - h * σ[N]
    @inbounds for i in 1:N-1
        E += -J * σ[i] * σ[i+1] - h * σ[i]
    end
    return E
end

"""
    magnetization(σ)

Return `sum(σ)` for a spin configuration.
"""
magnetization(σ::AbstractVector{<:Integer}) = sum(σ)

"""
    glauber_alg(J, T, σ; h=0.0, Nsteps=100)

Run `Nsteps` sweeps of single-spin-flip Glauber dynamics on `σ` (modified in
place) at temperature `T`, coupling `J`, and field `h`. Each sweep performs
`length(σ)` proposed flips.

Returns `(Es, Ms)`: vectors of length `Nsteps` with the energy and
magnetization recorded at the end of each sweep.
"""
function glauber_alg(J::Real, T::Real, σ::AbstractVector{<:Integer};
                     h::Real = 0.0, Nsteps::Integer = 100)
    N = length(σ)
    β = 1 / T
    Es = zeros(Nsteps)
    Ms = zeros(Nsteps)
    E = energy(J, σ; h = h)
    M = magnetization(σ)

    @inbounds for t in 1:Nsteps
        for _ in 1:N
            i  = rand(1:N)
            σl = σ[mod1(i - 1, N)]
            σr = σ[mod1(i + 1, N)]
            local_field = J * (σl + σr) + h
            if rand() < 0.5 * (1 - σ[i] * tanh(β * local_field))
                E  += 2 * σ[i] * local_field
                M  -= 2 * σ[i]
                σ[i] = -σ[i]
            end
        end
        Es[t] = E
        Ms[t] = M
    end
    return Es, Ms
end

end # module
