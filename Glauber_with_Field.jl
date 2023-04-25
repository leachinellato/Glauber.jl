using LinearAlgebra
using MKL
using Plots
using Random
using Statistics
using LaTeXStrings

function map_bool_to_integer(a::Bool)
    if a == false
        return b = -1
    else
        return b = 1
    end
end


function generate_RandomConf(N::Int64)
    aux = bitrand(N)
    σ = map_bool_to_integer.(aux)
    return σ
end


function rate_flip(J::Float64, β::Float64, h::Float64 ,σi::Int64, σpi::Int64, σmi::Int64)
    ω = 0.5 * (1.0 - σi * tanh(β*J * (σpi + σmi + h)))
    return ω
end


# function boltzmann(J::Float64, β::Float64, σ::BitVector)
#     N = length(σ)
#     sum = prod_bool(σ[N],σ[1]) #put eqaul to zero for operan boundary conditions
#     for i=1:N-1
#         sum += prod_bool(σ[i], σ[i+1])
#     end
#     P = (1/Z) * exp(β*J*sum)
#     return P
# end

function Energy(J::Float64, h::Float64, σ::Vector{Int64})
    N = length(σ)
    E = -J * σ[N] * σ[1] - h * σ[N]
    for i=1:N-1
        E += -J * σ[i] * σ[i+1] - h * σ[i]
    end
    return E
end

function Magnetization(σ::Vector{Int64})
    N = length(σ)
    M = 0
    for i=1:N
        M += σ[i]
    end
    return M
end



#######step 1##################################################
##Initialize the chain at infinite temperature               ##
##Calculate the inital energy and magnetization of the chain ##
###############################################################
#######step 2##################################################
##Choose a fix temperature T and calculate γ
##Pick a random spin along the chain, compute the energy 
## change ΔE and the rate.
###############################################################

function GlauberAlg(J::Float64, T::Float64,h::Float64, σ::Vector{Int64}; Nsteps=100)

    N = length(σ)
    Es = zeros(Nsteps)
    Ms = zeros(Nsteps)

    for jj = 1:Nsteps
    #one time step
        E = Energy(J,h, σ)
        M = Magnetization(σ)
        for _ = 1:N

            σev = deepcopy(σ)
            flip= rand(1:N)
            σev[flip] = -σev[flip]
    
            ΔE = Energy(J,h, σev) - Energy(J,h, σ)

            if flip == N 
                ω = rate_flip(J, 1/T,h, σ[flip], σ[1], σ[flip-1])    
            elseif flip == 1
                ω = rate_flip(J, 1/T, h, σ[flip], σ[flip+1], σ[N])    
            else 
                ω = rate_flip(J, 1/T, h, σ[flip], σ[flip+1], σ[flip-1])    
            end

            r = rand()
            if r < ω
                E = E + ΔE
                M = M - 2*σ[flip]
                σ = σev
            else
                nothing
            end
        end

        Es[jj] = E
        Ms[jj] = M
    end
    return Es, Ms
end




begin

    N = 500
    J = 1.0
    hs = collect(range(-1.0,1.0, step=0.1))
    T = 0.2

    σ0 = generate_RandomConf(N)
    Nsteps = 1000

    meanE = zeros(length(hs))
    meanM = zeros(length(hs))

    c = 0
    for h in hs
        c += 1
        En, Mag = GlauberAlg(J, T,h, σ0; Nsteps)
        meanE[c] = mean(En)
        meanM[c] = mean(Mag)
    end
    
    figE = plot(hs, meanE[:,1] ./ 500, xlabel =L"h[J]", ylabel =L"E[J]/N", label =L"N = 500", dpi=300)
    savefig(figE, "energy_with_h.png")
    figM = plot(hs, meanM[:,1] ./ 500, xlabel =L"h[J]", ylabel =L"M[J]/N", label =L"N = 500", dpi=300)
    savefig(figM, "magnetization_with_h.png")

    
end


   




