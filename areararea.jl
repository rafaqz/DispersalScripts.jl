# Area area dispersal formulation

H(x) = x > 0 ? 1 : 0

c(jx, jy, kx, ky) = begin
    r = sqrt((kx-jx)^2+(ky-jy)^2)
    θ2 = π / 2 - atan((ky - jy)/(kx - jx)) + π * H(jx - kx)
    r, θ2
end

car(jx2, jx1, jy2, jy1, kx, ky) = begin
    jx = (jx2 + jx1) / 2 
    jy = (jy2 + jy1) / 2 
    sqrt((kx-jx) ^ 2 + (ky-jy)^2)
end

caθ(jx2, jx1, jy2, jy1, kx, ky) = begin
    jx = (jx2 + jx1) / 2 
    jy = (jy2 + jy1) / 2 
    π / 2 - atan((ky - jy)/(kx - jx)) + π * H(jx - kx)
end

c((jx2 + jx1) / 2, (jy2 + jy1) / 2, kx, ky)
