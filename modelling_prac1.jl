cd("C:\\Users\\micho\\github\\Thesis")
pwd()

using ModelingToolkit, Plots, DifferentialEquations, LinearAlgebra
using Distributions
using Random


using Symbolics: scalarize

@variables t x(t)
@parameters τ
D = Differential(t)

sys = D(x) ~ (1 - x)/τ

@named fol_model = ODESystem(sys)

prob = ODEProblem(fol_model, [x => 0.0], (0.0, 10.0), [τ => 3.0])
sol = solve(prob)

plot(sol)

M = Matrix{Float64}(I, 4, 4)

M

Random.seed!(123)

C₁ = ones(2)

dir_dist = Dirichlet(2, 5)

C₁ = rand(dir_dist, 2)

C₁
