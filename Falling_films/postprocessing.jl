function reconstruct_U(U, parameters::GeneralParams)
    @unpack N1, N2, M1, M2 = parameters
    
    reconstr = reshape(U, (N1 + 1) * M2, 3)
    datah = reshape(reconstr[:,1], M2, N1 + 1)'
    dataq = reshape(reconstr[:,2], M2, N1 + 1)'
    datap = reshape(reconstr[:,3], M2, N1 + 1)'
    
    datah   = ffti(datah * M1 * M2, parameters)
    dataq   = ffti(dataq * M1 * M2, parameters)
    datap   = ffti(datap * M1 * M2, parameters)

    return datah, dataq, datap
end

function reconstruct_velocity_fields(U, kx, kz, ScaleH, ScaleU, parameters::GeneralParams)

    @unpack M1, M2, N1, N2, h0 = parameters

    sizeU = size([1, :]);

    # assumption: simplified model
    reconstru = reshape(U[1, :], (N1 + 1) * M2, 3)
    dataS1 = zeros(size(reconstru, 1), 1)
    dataS2 = dataS1;
    dataR1 = dataS2;
    dataR2 = dataS2;

    datah = reconstru[:,1];
    dataq = reconstru[:,2];
    datap = reconstru[:,3];

    datah = reshape(datah, M2, N1 + 1)'
    dataq = reshape(dataq, M2, N1 + 1)'
    datap = reshape(datap, M2, N1 + 1)'
    dataS1 = reshape(dataS1, M2, N1 + 1)'
    dataS2 = reshape(dataS2, M2, N1 + 1)'
    dataR1 = reshape(dataR1, M2, N1 + 1)'
    dataR2 = reshape(dataR2, M2, N1 + 1)'

    datah[1,1] = h0 + imag(datah[1,1])

    kx2 = kx .* kx;
    kz2 = kz .* kz;
    kxz = kx .* kz;
    k2 = kx2 + kz2;
    M1M2 = M1 * M2;

    # refactor: derivatives are duplicated in derivative.jl
    # create a re-useable derivatives function

    # x-derivatives
    datahx = kx .* imag(datah) - kx .* real(datah) * 1im;      
    dataqx = kx .* imag(dataq) - kx .* real(dataq) * 1im;
    datapx = kx .* imag(datap) - kx .* real(datap) * 1im;
    dataS1x = kx .* imag(dataS1) - kx .* real(dataS1) * 1im;
    dataS2x = kx .* imag(dataS2) - kx .* real(dataS2) * 1im;
    dataR1x = kx .* imag(dataR1) - kx .* real(dataR1) * 1im;
    dataR2x = kx .* imag(dataR2) - kx .* real(dataR2) * 1im;

    hx  = real(ffti(datahx * M1M2, parameters))
    px  = real(ffti(datapx * M1M2, parameters));
    qx  = real(ffti(dataqx * M1M2, parameters));
    s1x  = real(ffti(dataS1x * M1M2, parameters));
    s2x  = real(ffti(dataS2x * M1M2, parameters));
    r1x  = real(ffti(dataR1x * M1M2, parameters));
    r2x  = real(ffti(dataR2x * M1M2, parameters));


    # z-derivativess
    datahz = kz .* imag(datah) - kz .* real(datah) * 1im;
    dataqz = kz .* imag(dataq) - kz .* real(dataq) * 1im;
    datapz = kz .* imag(datap) - kz .* real(datap) * 1im;
    dataS1z = kz .* imag(dataS1) - kz .* real(dataS1) * 1im;
    dataS2z = kz .* imag(dataS2) - kz .* real(dataS2) * 1im;
    dataR1z = kz .* imag(dataR1) - kz .* real(dataR1) * 1im;
    dataR2z = kz .* imag(dataR2) - kz .* real(dataR2) * 1im;

    hZ  = real(ffti(datahz * M1M2, parameters));
    pz  = real(ffti(datapz * M1M2, parameters));
    qz  = real(ffti(dataqz * M1M2, parameters));
    s1z  = real(ffti(dataS1z * M1M2, parameters));
    s2z  = real(ffti(dataS2z * M1M2, parameters));
    r1z  = real(ffti(dataR1z * M1M2, parameters));
    r2z  = real(ffti(dataR2z * M1M2, parameters));


    # double x-derivatives 
    datahxx = -kx2 .* real(datah) - kx2 .* imag(datah) * 1im;
    dataqxx = -kx2 .* real(dataq) - kx2 .* imag(dataq) * 1im;
    datapxx = -kx2 .* real(datap) - kx2 .* imag(datap) * 1im; 

    # TODO: check accuracy of ffti here versus matlab
    # have seen small deviations
    hxx = real(ffti(datahxx * M1M2, parameters));
    pxx = real(ffti(datapxx * M1M2, parameters));
    qxx = real(ffti(dataqxx * M1M2, parameters));

    # double z-derivatives
    datahzz = -kz2 .* real(datah) - kz2 .* imag(datah) * 1im;
    dataqzz = -kz2 .* real(dataq) - kz2 .* imag(dataq) * 1im;
    datapzz = -kz2 .* real(datap) - kz2 .* imag(datap) * 1im; 
    
    hzz = real(ffti(datahzz * M1M2, parameters));
    pzz = real(ffti(datapzz * M1M2, parameters));
    qzz = real(ffti(dataqzz * M1M2, parameters));

    # cross derivatives
    datahxz = -kxz .* real(datah) - kxz .* imag(datah) * 1im;
    dataqxz = -kxz .* real(dataq) - kxz .* imag(dataq) * 1im;
    datapxz = -kxz .* real(datap) - kxz .* imag(datap) * 1im;

    # transform variables back to real space
    hxz = real(ffti(datahxz * M1M2, parameters));
    pxz = real(ffti(datapxz * M1M2, parameters));
    qxz = real(ffti(dataqxz * M1M2, parameters));

    # gradient of the laplacian of h (=curvature K) 
    datahLx = -(kx .* k2 .* imag(datah) - kx .* k2 .* real(datah) * 1im);
    datahLz = -(kz .* k2 .* imag(datah) - kz .* k2 .* real(datah) * 1im); 

    dhdt = -real(dataqx) - real(datapz) - imag(dataqx) * 1im - imag(datapz) * 1im;

    h   = real(ffti(datah * M1M2, parameters));
    P   = real(ffti(datap * M1M2, parameters));
    q   = real(ffti(dataq * M1M2, parameters));
    s1   = real(ffti(dataS1 * M1M2, parameters));
    s2   = real(ffti(dataS2 * M1M2, parameters));
    r1   = real(ffti(dataR1 * M1M2, parameters));
    r2   = real(ffti(dataR2 * M1M2, parameters));



    hmax = maximum(h)

    # slightly increase upper bound compared to matlab
    # so that upperbound is included

    # cf. 
    # julia: collect(0:0 .011178042476774046:1 .3413650972128854) -> 120 elements
    # matlab: 0:0 .011178042476774046:1 .3413650972128854 -> 121 elements
    y = reshape(collect(0:0.01 * hmax:hmax * 1.21), 1, :)
    y = y * ScaleH

    # determine yBar with datah in real space
    yBar = zeros(M1, M2, size(y, 2))
    for j = 1:1:M1
        for k = 1:1:M2
            yBar[j, k, :] = (y ./ h[j, k]);
        end
    end


    # calculate functions F0-F3 and G0-G3

    F0 = yBar - 1 / 2 * yBar.^2;
    F1 = yBar - 17 / 6 * yBar.^2 + 7 / 3 * yBar.^3 - 7 / 12 * yBar.^4;
    F2 = yBar - 13 / 2 * yBar.^2 + 57 / 4 * yBar.^3 - 111 / 8 * yBar.^4 + 99 / 16 * yBar.^5 - 33 / 32 * yBar.^6;
    F3 = yBar - 531 / 62 * yBar.^2 + 2871 / 124 * yBar.^3 - 6369 / 248 * yBar.^4 + 29601 / 2480 * yBar.^5 - 9867 / 4960 * yBar.^6;

    k, l = size(q)

    uField = similar(yBar)
    wField = similar(yBar)
    vField = similar(yBar)

    for k = 1:1:l
        F0Star = F0[:, k, :]
        F1Star = F1[:, k, :]
        F2Star = F2[:, k, :]
        j, m = size(y)
        for j = 1:1:m
            uField[:, k, j] = 3 * (q[:, k] - s1[:, k] - s2[:, k]) ./ h[:, k] .* F0Star[:, j] + 45 * s1[:, k] ./ h[:, k] .* F1Star[:, j] + 210 * s2[:, k] ./ h[:, k] .* F2Star[:, j]
            wField[:, k, j] = 3 * (P[:, k] - r1[:, k] - r2[:, k]) ./ h[:, k] .* F0Star[:, j] + 45 * r1[:, k] ./ h[:, k] .* F1Star[:, j] + 210 * r2[:, k] ./ h[:, k] .* F2Star[:, j]
            vField[:, k, j] = (y[:, j].^2 .* (-24 .* h[:, k].^6 .* (pz[:, k] + qx[:, k] + 14 .* r1z[:, k] + 69 .* r2z[:, k] + 14 .* s1x[:, k] + 69 .* s2x[:, k]) - 3465 .* (hZ[:, k] .* r2[:, k] + hx[:, k] .* s2[:, k]) .* y[:, j].^5 + 495 .* h[:, k] .* y[:, j].^4 .* (42 .* hZ[:, k] .* r2[:, k] + 42 .* hx[:, k] .* s2[:, k] + r2z[:, k] .* y[:, j] + s2x[:, k] .* y[:, j]) 
        - 105 .* h[:, k].^2 .* y[:, j].^3 .* (4 .* hZ[:, k] .* r1[:, k] + 444 .* hZ[:, k] .* r2[:, k] + 4 .* hx[:, k] .* s1[:, k] + 444 .* hx[:, k] .* s2[:, k] + 33 .* r2z[:, k] .* y[:, j] + 33 .* s2x[:, k] .* y[:, j]) 
        + 84 .* h[:, k].^3 .* y[:, j].^2 .* (20 .* hZ[:, k] .* r1[:, k] + 570 .* hZ[:, k] .* r2[:, k] + 20 .* hx[:, k] .* s1[:, k] + 570 .* hx[:, k] .* s2[:, k] + r1z[:, k] .* y[:, j] + 111 .* r2z[:, k] .* y[:, j] + s1x[:, k] .* y[:, j] + 111 .* s2x[:, k] .* y[:, j]) + 8 .* h[:, k].^5 .* (6 .* hZ[:, k] .* P[:, k] + 6 .* hx[:, k] .* q[:, k] + 84 .* hZ[:, k] .* r1[:, k] + 414 .* hZ[:, k] .* r2[:, k] + 84 .* hx[:, k] .* s1[:, k] + 414 .* hx[:, k] .* s2[:, k] + pz[:, k] .* y[:, j] + qx[:, k] .* y[:, j] + 84 .* r1z[:, k] .* y[:, j] + 909 .* r2z[:, k] .* y[:, j] + 84 .* s1x[:, k] .* y[:, j] + 909 .* s2x[:, k] .* y[:, j]) 
        - 6 .* h[:, k].^4 .* y[:, j] .* (4 .* hZ[:, k] .* P[:, k] + 4 .* hx[:, k] .* q[:, k] + 336 .* hZ[:, k] .* r1[:, k] + 3636 .* hZ[:, k] .* r2[:, k] + 336 .* hx[:, k] .* s1[:, k] + 3636 .* hx[:, k] .* s2[:, k] + 70 .* r1z[:, k] .* y[:, j] + 1995 .* r2z[:, k] .* y[:, j] + 70 .* s1x[:, k] .* y[:, j] + 1995 .* s2x[:, k] .* y[:, j]))) ./ (16 .* h[:, k].^8)       
    
    
        end
    end

    uField = uField * ScaleU;
    vField = vField * ScaleU;
    wField = wField * ScaleU;
    h = h * ScaleH;

    k, l, m = size(uField)

    alpha = similar(yBar)
    for j = 1:1:M1
        for k = 1:1:M2
            for l = 1:1:m
                if y[l] > h[j, k]
                    uField[j, k, l] = NaN;
                    vField[j, k, l] = NaN;
                    wField[j, k, l] = NaN;
                    alpha[j, k, l] = 0;
                else
                    alpha[j, k, l] = 1;
                end
            end
        end
    end

    yBar = 1
    F0 = yBar - 1 / 2 * yBar.^2
    F1 = yBar - 17 / 6 .* yBar.^2 + 7 / 3 .* yBar.^3 - 7 / 12 .* yBar.^4
    F2 = yBar - 13 / 2 .* yBar.^2 + 57 / 4 * yBar.^3 - 111 / 8 .* yBar.^4 + 99 / 16 .* yBar.^5 - 33 / 32 .* yBar.^6
    F3 = yBar - 531 / 62 .* yBar.^2 + 2871 / 124 * yBar.^3 - 6369 / 248 .* yBar.^4 + 29601 / 2480 .* yBar.^5 - 9867 / 4960 .* yBar.^6


    ac, l = size(q)
    uSurf = zeros(ac, l)
    wSurf = similar(uSurf)
    vSurf = similar(uSurf)
    cWave = similar(uSurf)

    for k = 1:1:l
        F0Star = F0;
        F1Star = F1;
        F2Star = F2;
    
        j, m = size(y)
    
        j = 1

    # matlab -> julia
    # F0Star[:, j] -> F0Star
    # F1Star[:, j] -> F0Star
    # F2Star[:, j] -> F0Star
        uSurf[:,k] = 3 * (q[:, k] - s1[:, k] - s2[:, k]) ./ h[:, k] .* F0Star + 45 * s1[:, k] ./ h[:, k] .* F1Star + 210 * s2[:, k] ./ h[:, k] .* F2Star
        wSurf[:, k] = 3 * (P[:, k] - r1[:, k] - r2[:, k]) ./ h[:, k] .* F0Star + 45 * r1[:, k] ./ h[:, k] .* F1Star + 210 * r2[:, k] ./ h[:, k] .* F2Star
        vSurf[:,k] = (h[:,k].^2 .* (-24 .* h[:,k].^6 .* (pz[:,k] + qx[:,k] + 14 .* r1z[:,k] + 69 .* r2z[:,k] + 14 .* s1x[:,k] + 69 .* s2x[:,k]) - 3465 .* (hZ[:,k] .* r2[:,k] + hx[:,k] .* s2[:,k]) .* h[:,k].^5 + 495 .* h[:,k] .* h[:,k].^4 .* (42 .* hZ[:,k] .* r2[:,k] + 42 .* hx[:,k] .* s2[:,k] + r2z[:,k] .* h[:,k] + s2x[:,k] .* h[:,k])
    - 105 .* h[:,k].^2 .* h[:,k].^3 .* (4 .* hZ[:,k] .* r1[:,k] + 444 .* hZ[:,k] .* r2[:,k] + 4 .* hx[:,k] .* s1[:,k] + 444 .* hx[:,k] .* s2[:,k] + 33 .* r2z[:,k] .* h[:,k] + 33 .* s2x[:,k] .* h[:,k])
    + 84 .* h[:,k].^3 .* h[:,k].^2 .* (20 .* hZ[:,k] .* r1[:,k] + 570 .* hZ[:,k] .* r2[:,k] + 20 .* hx[:,k] .* s1[:,k] + 570 .* hx[:,k] .* s2[:,k] + r1z[:,k] .* h[:,k] + 111 .* r2z[:,k] .* h[:,k] + s1x[:,k] .* h[:,k] + 111 .* s2x[:,k] .* h[:,k])
    + 8 .* h[:,k].^5 .* (6 .* hZ[:,k] .* P[:,k] + 6 .* hx[:,k] .* q[:,k] + 84 .* hZ[:,k] .* r1[:,k] + 414 .* hZ[:,k] .* r2[:,k] + 84 .* hx[:,k] .* s1[:,k] + 414 .* hx[:,k] .* s2[:,k] + pz[:,k] .* h[:,k] + qx[:,k] .* h[:,k] + 84 .* r1z[:,k] .* h[:,k] + 909 .* r2z[:,k] .* h[:,k] + 84 .* s1x[:,k] .* h[:,k] + 909 .* s2x[:,k] .* h[:,k]) 
    - 6 .* h[:,k].^4 .* h[:,k] .* (4 .* hZ[:,k] .* P[:,k] + 4 .* hx[:,k] .* q[:,k] + 336 .* hZ[:,k] .* r1[:,k] + 3636 .* hZ[:,k] .* r2[:,k] + 336 .* hx[:,k] .* s1[:,k] + 3636 .* hx[:,k] .* s2[:,k] + 70 .* r1z[:,k] .* h[:,k] + 1995 .* r2z[:,k] .* h[:,k] + 70 .* s1x[:,k] .* h[:,k] + 1995 .* s2x[:,k] .* h[:,k]))) ./ (16 .* h[:,k].^8)

        for ls = 1:ac
            if (hx[ls, k] > 0.005)
                cWave[ls, k] = (hx[ls, k] .* uSurf[ls, k] - vSurf[ls, k]) ./ hx[ls, k];
            else
                cWave[ls, k] = NaN;
            end
        end
    end

    return uField, vField, wField, uSurf, vSurf, hx, y, cWave, h, alpha
end