using LinearAlgebra
using MKL
using Plots
using Random
using Statistics

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


function rate_flip(J::Float64, β::Float64, σi::Int64, σpi::Int64, σmi::Int64)
    γ = 2.0 * tanh(2*β*J)
    ω = 0.5 * (1.0 - 0.5 * γ *(σi * (σmi + σpi)))
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

function Energy(J::Float64, σ::Vector{Int64})
    N = length(σ)
    E = σ[N] * σ[1]
    for i=1:N-1
        E += σ[i] * σ[i+1]
    end
    E = -J * E
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

function GlauberAlg(J::Float64, T::Float64, σ::Vector{Int64}; Nsteps=100)

    N = length(σ)
    Es = zeros(Nsteps)
    Ms = zeros(Nsteps)

    for jj = 1:Nsteps
    #one time step
        E = Energy(J, σ)
        M = Magnetization(σ)
        for _ = 1:N

            σev = deepcopy(σ)
            flip= rand(1:N)
            σev[flip] = -σev[flip]
    
            ΔE = Energy(J, σev) - Energy(J, σ)

            if flip == N 
                ω = rate_flip(J, 1/T, σ[flip], σ[1], σ[flip-1])    
            elseif flip == 1
                ω = rate_flip(J, 1/T, σ[flip], σ[flip+1], σ[N])    
            else 
                ω = rate_flip(J, 1/T, σ[flip], σ[flip+1], σ[flip-1])    
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

    N = 10
    J = 1.0
    T = 0.1

    σ0 = generate_RandomConf(N)

    #calculate the time dependence of the energy and magnetization for some fixed T
    Nsteps = 100
    En, Mag = GlauberAlg(J, T, σ0; Nsteps)
    
    plot(Mag)

    meanE = mean(En)
    meanM = mean(Mag)

    
end


   
begin

    N = 500
    J = 1.0
    Nsteps = 1000

    #calculate for severals temperatures
    Ts = collect(range(0.1,4.0, step=0.1))
    meanE = zeros(length(Ts))
    meanM = zeros(length(Ts))

    cc = 0
    for T in Ts
        cc += 1
        σ0 = generate_RandomConf(N)
        En, Mag = GlauberAlg(J, T, σ0; Nsteps)
        meanE[cc] = mean(En)
        meanM[cc] = mean(Mag)
    end

    plot(Ts, meanE / N)
    plot(Ts, abs.(meanM)  / N)
    
end

