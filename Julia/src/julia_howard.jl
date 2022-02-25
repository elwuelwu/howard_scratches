module julia_howard

using DSP
using Random
using LinearAlgebra
using Hadamard
using Kronecker

mm = 4
N = 2^mm

smat = rand((0,1), mm, mm)
smat = (smat + transpose(smat)) .% 2 + Diagonal(rand((0,1), mm))
bvec = rand((0,1), mm)

#display(smat)
display(bvec)

chirp = zeros(ComplexF64, N)

for i in 1:N
    xvec = reverse(digits(i-1, base=2, pad=mm))
    chirp[i] = im^(xvec'*smat*xvec + 2*bvec'*xvec)
end

# Decoder

X = [0 1; 1 0]
I = [1 0; 0 1]
seq = [X]
for i in 1:(mm-1); push!(seq, I); end

global perms = []
for i in 1:mm
    push!(perms, kron(seq...))
    global seq = circshift(seq, 1)
end

function decoderow(vecin, row)
    shiftmultiplied = conj(vecin) .* (perms[row]*vecin)
    wht = fwht_natural(shiftmultiplied)
    avec = zeros(mm); avec[row] = 1
    phasevec = [im^(avec'*digits(i-1, base=2, pad=mm)) for i in 1:N]
    display(wht .* phasevec)
    abswht = map(abs, wht)
    ix = argmax(abswht)
    row = reverse(digits(ix-1, base=2, pad=mm))
    return row
end

sest = zeros(Int8, mm, mm)
for row in 1:mm
    sest[row, :] = decoderow(chirp, row)
end




end # module