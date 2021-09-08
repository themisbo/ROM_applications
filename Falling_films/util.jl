function ffti(dat, parameters::GeneralParams)
    @unpack N1, M2 = parameters
    dat = vcat(dat, [conj(reverse(dat[2:N1, 1], dims=1)) conj(rotr90(rotr90(dat[2:N1,2:M2])))])
    dat = ifft(dat)
    return dat
end