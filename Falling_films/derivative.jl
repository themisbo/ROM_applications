using Parameters

export dUdt_simplified_model

function dUdt_simplified_model(yU, k2cut, N, parameters::GeneralParams)
    @unpack N1, M1, M2, h0, eta, zeta, delta = parameters;
        
    yU = reshape(yU, length(yU), :);  # assume yU is a column vector

    reconstru = reshape(yU, (N1 + 1) * M2, 3);
    datah = reconstru[:, 1];
    dataq = reconstru[:, 2];
    datap = reconstru[:, 3];

    dataq = reshape(dataq, M2, N1 + 1)';
    datah = reshape(datah, M2, N1 + 1)';
    datap = reshape(datap, M2, N1 + 1)';

    # fix the value of the average film thickness
    datah[1,1] = h0 + imag(datah[1, 1]);

    kx2 = kx .* kx;
    kz2 = kz .* kz;
    kxz = kx .* kz;
    k2 = kx2 + kz2;

    M1M2 = M1 * M2;

    # x-derivatives
    datahx = kx .* imag(datah) - kx .* real(datah) * 1im;      
    dataqx = kx .* imag(dataq) - kx .* real(dataq) * 1im;
    datapx = kx .* imag(datap) - kx .* real(datap) * 1im;


    # transform variables back to real space
    hx  = real(ffti(datahx * M1M2, parameters));
    px  = real(ffti(datapx * M1M2, parameters));
    qx  = real(ffti(dataqx * M1M2, parameters));

    # z-derivativess
    datahz = kz .* imag(datah) - kz .* real(datah) * 1im;
    dataqz = kz .* imag(dataq) - kz .* real(dataq) * 1im;
    datapz = kz .* imag(datap) - kz .* real(datap) * 1im;

    # transform variables back to real space
    hZ  = real(ffti(datahz * M1M2, parameters));
    pz  = real(ffti(datapz * M1M2, parameters));
    qz  = real(ffti(dataqz * M1M2, parameters));

    # double x-derivatives 
    datahxx = -kx2 .* real(datah) - kx2 .* imag(datah) * 1im;
    dataqxx = -kx2 .* real(dataq) - kx2 .* imag(dataq) * 1im;
    datapxx = -kx2 .* real(datap) - kx2 .* imag(datap) * 1im;  

    # transform variables back to real space  
    hxx = real(ffti(datahxx * M1M2, parameters));
    pxx = real(ffti(datapxx * M1M2, parameters));
    qxx = real(ffti(dataqxx * M1M2, parameters));


    # double z-derivatives 
    datahzz = -kz2 .* real(datah) - kz2 .* imag(datah) * 1im;
    dataqzz = -kz2 .* real(dataq) - kz2 .* imag(dataq) * 1im;
    datapzz = -kz2 .* real(datap) - kz2 .* imag(datap) * 1im; 

    # transform variables back to real space
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

    # transform variables back to real space
    h = real(ffti(datah * M1M2, parameters));
    P = real(ffti(datap * M1M2, parameters));
    q = real(ffti(dataq * M1M2, parameters));

    datahLx = ffti(datahLx * M1M2, parameters);
    datahLz = ffti(datahLz * M1M2, parameters);

    h2 = h .* h;
    h3 = h2 .* h;

    Ineq = 0;  # switch to 1 for second-order corrections
    A = 0
    hs = 0
    dqdt = (
                0.011904761904761904 .* (
                    8.0 * h .* (
                        4900.0 * A .* hs .* hx - 735.0 * eta .* (13.0 * hx .* hZ .* P + 16.0 * hx .* hx .* q + 3.0 * hZ .* hZ .* q) 
                        + 6.0 * q .* (1225.0 .+ 9.0 * delta .* (-70.0 * hZ .* P + hx .* q .* (-70.0 .+ delta .* hx .* Ineq .* q)))
                    ) 
                    + 6.0 * h2 .* (
                        245.0 * eta .* (
                            73.0 * hxz .* P + 43.0 * hZ .* px + 13.0 * hx .* pz + 96.0 * hxx .* q + 23.0 * hzz .* q + 72.0 * hx .* qx + 16.0 * hZ .* qz
                        ) + 8.0 * delta .* (560.0 * pz .* q - 17.0 * q .* (-70.0 .+ delta .* hx .* Ineq .* q) .* qx + 630.0 * P .* qz)
                    ) 
                    - 5880.0 * h3 .* (2.0 * eta .* (7.0 * pxz + 9.0 * qxx + 2.0 * qzz)) 
                    - 300.0 * hx .* (196.0 * A .* hs .* hs) 
                    - 35.0 * h2 .* h2 .* (560.0 .- 560.0 * hx .* zeta)
                    - 19600.0 * h2 .* h2 .* real(datahLx)
                )
            ) ./ (delta .* h2 .* (4.0 * h .* (-70.0 .+ delta .* hx .* Ineq .* q)));


    dpdt = (
                0.002976190476190476 .* (
                    -280.0 * A .* (2.0 * h .- 3.0 * hs) .* hs .* hZ 
                    + h .* (
                        12.0 * (
                            P .* (-70.0 .+ 7.0 * eta .* (3.0 * hx .* hx + 16.0 * hZ .* hZ) + 36.0 * delta .* hZ .* P)
                            + hx .* (91.0 * eta .* hZ + 36.0 * delta .* P) .* q
                        ) 
                        - 3.0 * h .* (
                            16.0 * delta .* (17.0 * P .* pz + 9.0 * px .* q + 8.0 * P .* qx)
                            + 7.0 * eta .* (23.0 * hxx .* P + 96.0 * hzz .* P + 16.0 * hx .* px + 72.0 * hZ .* pz + 73.0 * hxz .* q + 13.0 * hZ .* qx + 43.0 * hx .* qz)
                        )
                        + 84.0 * h2 .* (2.0 * eta .* (2.0 * pxx + 9.0 * pzz + 7.0 * qxz)) - 280.0 * h3 .* hZ .* zeta 
                    )
                    + 280.0 * h2 .* h2 .* real(datahLz)
                )
            ) ./ (delta .* h3);

        # fourier space
    dqdt = fft(dqdt);
    dpdt = fft(dpdt);

        # alaising filter removes high frequency modes
    dqdt = dqdt .* k2cut / M1M2;
    dpdt = dpdt .* k2cut / M1M2;

        # transform data back to state vector
    dhdt = conj(dhdt);
    dqdt = conj(dqdt);
    dpdt = conj(dpdt);

    dhdt = reshape(dhdt', 1, :);
    dqdt = reshape(dqdt', 1, :);
    dpdt = reshape(dpdt', 1, :);

        # required derivative 
    dUdt = reshape([dhdt[1:N ÷ 2]; dqdt[1:N ÷ 2]; dpdt[1:N ÷ 2]], 1, :)';
    return dUdt
end