
function create_initial_condition(kx, kz, k2cut, N, parameters::GeneralParams)
    datah, datap, dataq, U = initial_condition(parameters);

    datas1s2r1r2 = zeros(length(datah))'

    # remove high frequency parts by aliasing filter
    datah = datah .* k2cut;
    datap = datap .* k2cut;
    dataq = dataq .* k2cut;

    # reshape data for U vector
    datahSV = reshape(datah', 1, :)
    datapSV = reshape(datap', 1, :)
    dataqSV = reshape(dataq', 1, :)

    U = vcat(datahSV[1:N ÷ 2], dataqSV[1:N ÷ 2], datapSV[1:N ÷ 2])
    U = reshape(U, 1, :)
    return U
end

function initial_condition(parameters::GeneralParams)
    @unpack NN1, NN2 = parameters
    @unpack stSine, spSine, stCosine, spCosine, stNwaves, spNwaves, stNoise, spNoise = parameters
    @unpack k0x, k0z, h0, Lx, Lz = parameters

    Ax = stSine
    Az = spSine
    Bx = stCosine
    Bz = spCosine
    nwx = stNwaves
    nwz = spNwaves

    if (NN1 == 0) k0x = 0; Lx = 0; Ax = 0; Bx = 0 end
    if (NN2 == 0) k0z = 0; Lz = 0; Az = 0; Bz = 0 end

    N1 = 2^NN1; M1 = 2 * N1; N2 = 2^NN2; M2 = 2 * N2; N = 4 * N2 * (N1 + 1);

    datah = zeros((M1, M2))
    dataq = zeros((M1, M2))
    datap = zeros((M1, M2))

    randz = zeros(M2)
    for n = 1:M2
        randz[n] = spNoise * rand()
    end

    for m = 1:M1
        x = (m - 1) * Lx / M1
        for n = 1:M2
            z = (n - 1) * Lz / M2
            datah[m, n] = h0 * (1 + Ax * sin(nwx * k0x * x) + Bx * cos(nwx * k0x * x) + Az * sin(nwz * k0z * z) + Bz * cos(nwz * k0z * z) + stNoise * rand()) + randz[n]
            datap[m, n] = 0
            dataq[m, n] = datah[m, n].^3 / 3
        end
    end

    # calculate mass
    mass = sum(datah)
    mass = mass / (M1 * M2)

    # set mass to h0 for random function
    if spNoise != 0
        datah = datah - mass + h0
    end

    # Fourier transform of the initial condition
    datah = fft(datah)
    datap = fft(datap)
    dataq = fft(dataq)

    datah = conj(datah)
    datap = conj(datap)
    dataq = conj(dataq)

    datahSV = vec(datah')
    datapSV = vec(datap')
    dataqSV = vec(dataq')

    M1M2 = M1 * M2;

    # note ÷ for integer division
    U = vcat(datahSV[1:N ÷ 2], dataqSV[1:N ÷ 2], datapSV[1:N ÷ 2]) ./ M1M2 

    datah = datah ./ M1M2;
    dataq = dataq ./ M1M2;
    datap = datap ./ M1M2;

    return datah, datap, dataq, U

end

function aliasing_filter(radius::Float64, parameters::GeneralParams)   
    @unpack k0x, k0z, N1, N2, M1, M2 = parameters;

    m = round(radius * N1)
    n = round(radius * N2)
    k2cutCriterion = (m * k0x)^2 + ((-N2 + mod(n + N2, M2)) * k0z)^2
    kx, kz, k2, k2cut = [zeros((N1 + 1, M2)) for i = 1:4]  # initialise
    for m = 1:N1 + 1
        for n = 1:M2
            kx[m, n] = (m - 1) * k0x
            kz[m, n] = (-N2 + mod((n - 1) + N2, M2)) * k0z
            k2[m, n] = kx[m, n] * kx[m, n] + kz[m, n] * kz[m, n]
            if (k2[m, n] > k2cutCriterion)
                k2cut[m, n] = 0
            else
                k2cut[m, n] = 1
            end
        end
    end
    k2cut = [k2cut; [conj(reverse(k2cut[2:N1,1], dims=1)) conj(reverse(reverse(k2cut[2:N1,2:M2], dims=1), dims=2))]]
    return kx, kz, k2cut
end