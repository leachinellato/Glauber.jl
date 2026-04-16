using Test
using Random

include(joinpath(@__DIR__, "..", "src", "Glauber.jl"))
using .Glauber

@testset "Glauber" begin
    @testset "energy / magnetization on known configs" begin
        N = 10
        all_up   = fill(1, N)
        all_down = fill(-1, N)
        alt      = [isodd(i) ? 1 : -1 for i in 1:N]

        @test energy(1.0, all_up)                == -N
        @test energy(1.0, all_down)              == -N
        @test energy(1.0, alt)                   ==  N        # all bonds antiparallel
        @test energy(1.0, all_up; h = 0.5)       == -N - 0.5N
        @test magnetization(all_up)              ==  N
        @test magnetization(all_down)            == -N
        @test magnetization(alt)                 ==  0
    end

    @testset "generate_random_conf" begin
        Random.seed!(42)
        σ = generate_random_conf(100)
        @test length(σ) == 100
        @test all(x -> x == 1 || x == -1, σ)
    end

    @testset "glauber_alg output shape and energy bookkeeping" begin
        Random.seed!(0)
        N, Nsteps = 20, 50
        σ = generate_random_conf(N)
        Es, Ms = glauber_alg(1.0, 1.5, σ; Nsteps = Nsteps)
        @test length(Es) == Nsteps
        @test length(Ms) == Nsteps
        # After the run, the final recorded E/M should match a full
        # recomputation from the evolved σ — this verifies the incremental
        # bookkeeping of ΔE / ΔM inside the inner loop.
        @test isapprox(Es[end], energy(1.0, σ); atol = 1e-10)
        @test Ms[end] == magnetization(σ)
    end

    @testset "local ΔE matches full recomputation" begin
        Random.seed!(7)
        for _ in 1:20
            N = 15
            σ = generate_random_conf(N)
            J, h = 1.0, 0.3
            i = rand(1:N)
            E0 = energy(J, σ; h = h)
            σ[i] = -σ[i]
            E1 = energy(J, σ; h = h)
            σ[i] = -σ[i]  # restore
            left  = σ[mod1(i - 1, N)]
            right = σ[mod1(i + 1, N)]
            ΔE_local = 2 * σ[i] * (J * (left + right) + h)
            @test isapprox(E1 - E0, ΔE_local; atol = 1e-12)
        end
    end

    @testset "low T drives energy below random baseline" begin
        Random.seed!(123)
        N = 200
        σ = generate_random_conf(N)
        Es, _ = glauber_alg(1.0, 0.1, σ; Nsteps = 200)
        # Random start at high T has ⟨E⟩/N ≈ 0. After equilibrating at
        # T = 0.1 (deep in the ordered correlation regime for 1D Ising),
        # almost all nearest-neighbour bonds are satisfied.
        @test Es[end] / N < -0.7
    end
end
