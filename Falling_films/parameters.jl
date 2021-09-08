using Parameters

@with_kw struct GeneralParams   

    # input parameters (Nusselt scaling)
    Re::Float64 = 40.8
    Ka::Float64 = 3923.0
    Ct::Float64 = 0.0
    k_x::Float64 = 0.0729
    k_z::Float64 = 0.0729
    h0::Float64 = 0.9147

    # derived (Shkadov scaling)
    eta::Float64 = (3.0 * Re)^(4 / 9) / Ka^(2 / 3)
    delta::Float64 = (3.0 * Re)^(11 / 9) / Ka^(1 / 3)
    zeta::Float64 = Ct * (3.0 * Re)^(2 / 9) / Ka^(1 / 3)
    k0x::Float64 = k_x / sqrt(eta)
    k0z::Float64 = k_z / sqrt(eta)

    Lx::Float64 = 2 * pi / k0x
    Lz::Float64 = 2 * pi / k0z

    # initial conditions
    stSine::Float64 = 0.2
    spSine::Float64 = 0.2
    stCosine::Int64 = 0
    spCosine::Int64 = 0
    stNwaves::Int64 = 1
    spNwaves::Int64 = 1
    stNoise::Int64 = 0
    spNoise::Int64 = 0
    
    # mesh parameters
    # 2D simulation
    # NN1::Int64 = 5
    # NN2::Int64 = 1
    # N1::Int64 = 32
    # N2::Int64 = 1
    # M1::Int64 = 64
    # M2::Int64 = 2

    # 3D simulation
    NN1::Int64 = 5
    NN2::Int64 = 5
    N1::Int64 = 32
    N2::Int64 = 32
    M1::Int64 = 64
    M2::Int64 = 64
end

parameters = GeneralParams()