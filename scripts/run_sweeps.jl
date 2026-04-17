using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

include(joinpath(@__DIR__, "..", "src", "Glauber.jl"))
using .Glauber
using Plots
using LaTeXStrings
using Random
using Statistics

Random.seed!(1)

const J = 1.0
const OUTDIR = joinpath(@__DIR__, "..")

# ---------------------------------------------------------------------------
# Sweep over temperature at h = 0 for three system sizes.
# ---------------------------------------------------------------------------
function sweep_temperature()
    Ns      = (5, 50, 500)
    Nsteps  = 1000
    Ts      = collect(range(0.1, 4.0; step = 0.1))
    meanE   = zeros(length(Ts), length(Ns))
    meanM   = zeros(length(Ts), length(Ns))

    for (ci, N) in enumerate(Ns), (ti, T) in enumerate(Ts)
        σ = generate_random_conf(N)
        Es, Ms = glauber_alg(J, T, σ; Nsteps = Nsteps)
        meanE[ti, ci] = mean(Es)
        meanM[ti, ci] = mean(Ms)
    end

    figE = plot(Ts, meanE[:, 1] ./ Ns[1]; xlabel = L"T[J]", ylabel = L"E[J]/N",
                label = L"N = 5", dpi = 300)
    plot!(figE, Ts, meanE[:, 2] ./ Ns[2]; label = L"N = 50")
    plot!(figE, Ts, meanE[:, 3] ./ Ns[3]; label = L"N = 500")
    savefig(figE, joinpath(OUTDIR, "energy.png"))

    figM = plot(Ts, meanM[:, 1] ./ Ns[1]; xlabel = L"T[J]", ylabel = L"M[J]/N",
                label = L"N = 5", dpi = 300)
    plot!(figM, Ts, meanM[:, 2] ./ Ns[2]; label = L"N = 50")
    plot!(figM, Ts, meanM[:, 3] ./ Ns[3]; label = L"N = 500")
    savefig(figM, joinpath(OUTDIR, "magnetization.png"))
end

# ---------------------------------------------------------------------------
# Sweep over external field h at fixed temperature for N = 500.
# ---------------------------------------------------------------------------
function sweep_field()
    N      = 500
    Nsteps = 1000
    T      = 0.2
    hs     = collect(range(-1.0, 1.0; step = 0.1))
    meanE  = zeros(length(hs))
    meanM  = zeros(length(hs))

    for (i, h) in enumerate(hs)
        σ = generate_random_conf(N)
        Es, Ms = glauber_alg(J, T, σ; h = h, Nsteps = Nsteps)
        meanE[i] = mean(Es)
        meanM[i] = mean(Ms)
    end

    figE = plot(hs, meanE ./ N; xlabel = L"h[J]", ylabel = L"E[J]/N",
                label = L"N = 500", dpi = 300)
    savefig(figE, joinpath(OUTDIR, "energy_with_h.png"))

    figM = plot(hs, meanM ./ N; xlabel = L"h[J]", ylabel = L"M[J]/N",
                label = L"N = 500", dpi = 300)
    savefig(figM, joinpath(OUTDIR, "magnetization_with_h.png"))
end

sweep_temperature()
sweep_field()
