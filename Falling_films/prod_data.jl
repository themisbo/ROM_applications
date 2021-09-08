println("INFO: Loading packages")

using FFTW
using Parameters
using DifferentialEquations
using Plots
using Statistics
using PyCall

include("./parameters.jl")
include("./derivative.jl")
include("./initial.jl")
include("./util.jl")
include("./postprocessing.jl")

np = pyimport("numpy")

num_of_samples = 20

Re_range = range(10.0, stop = 69.0, length = num_of_samples)

struct SolverParameters
    k2cut
    N
    parameters
end

function simplified_model!(du, u, p, t)
    # unpack parameters
    k2cut = p.k2cut
    N = p.N
    parameters = p.parameters

    # alias for U
    # refactor: use lowercase `u` throughout 
    U = u

    # update derivative
    du[:] = dUdt_simplified_model(U, k2cut, N, parameters)
end

T_MAX = 25.0
radius = 0.6600
tspan = (0.0, T_MAX)
ScaleH = 1
ScaleU = 1
no_t = 50
t_range = LinRange(0.0, T_MAX, no_t)
DIMENSION = 3

Threads.@threads for sample_val in 1:num_of_samples

    global parameters = GeneralParams(Re = Re_range[sample_val], Ka = 1000.)

    N = 4 * parameters.N2 * (parameters.N1 + 1);
    global kx, kz, k2cut = aliasing_filter(radius, parameters)

    u0 = create_initial_condition(kx, kz, k2cut, N, parameters)

    p = SolverParameters(k2cut, N, parameters)

    prob = ODEProblem(simplified_model!, u0, tspan, p)

    println("INFO: Solving equation")
    @time sol = solve(prob, Tsit5(), maxiters = 10000)

    println("INFO: Postprocessing results")

    @unpack M1, M2 = parameters

    # generate plots using Nusselt scaling
    ScaleX = 2 * pi / parameters.k0x
    ScaleZ = 2 * pi / parameters.k0z

    x = (0:1 / (M1 - 1):1) * ScaleX
    z = (0:1 / (M2 - 1):1) * ScaleZ

    lll = zeros(Float64, no_t, 64, 64)
    for i in 1:no_t
        t = t_range[i]
        U = sol(t)    
        datah, dataq, datap = reconstruct_U(U, parameters)
        surf = real(datah) * ScaleH
        lll[i,:,:] = surf
    end

    nplll = np.array(lll)

    if sol.retcode == :MaxIters
        println("Unsuccesful iteration.")
    else
        println("Iteration:", sample_val)
        np.save("data/zvals_Re" * string(Re_range[sample_val]) * ".npy", nplll)
    end
end
