module julia_howard

using DSP
using Random
using LinearAlgebra
using Hadamard
using Kronecker

mm = 4
N = 2^mm
IdN = Diagonal(ones(N))

smat = rand((0,1), mm, mm)
smat = (smat + transpose(smat)) .% 2 + Diagonal(rand((0,1), mm))
bvec = rand((0,1), mm)

display(smat)
println("")
display(bvec)
println("")

function generate_chirp(smat, bvec)
    chirp = zeros(ComplexF64, N)
    for i in 1:N
        xvec = reverse(digits(i-1, base=2, pad=mm))
        chirp[i] = im^(xvec'*smat*xvec + 2*bvec'*xvec)
    end
    return chirp
end

chirp = generate_chirp(smat, bvec)

# Decoder

Z = [1 0; 0 -1]
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
    # Phase vector that rotates the signs to +- 1
    phasevec = [im^(-avec'*reverse(digits(i-1, base=2, pad=mm))) for i in 1:N]
    abswht = map(abs, wht)
    ix = argmax(abswht)
    row = reverse(digits(ix-1, base=2, pad=mm))
    # Determine the bit on the b-vector
    b_bit = (1 - sign((phasevec .* wht)[ix].re)) / 2
    return row, b_bit
end

function construct_pauli(avec, bvec)
    result = [1]
    for i in 1:length(avec)
        abit = avec[i]
        bbit = bvec[i]
        result = kron(result, im^(abit*bbit)*X^abit*Z^bbit)
    end
    return result
end

function project(vecin, eigenvalue, a, aS)
    # The binary chirps are common +-1 eigenvectors for E(a,aS) for all a
    # Thus an operator (1/2)*(Identity + eigenvalue*E(a,aS)) should keep the vector intact
    pauli = construct_pauli(a,aS)
    projector = (1/2)*(IdN + eigenvalue*pauli)
    return projector*vecin
end

sest = zeros(Int8, mm, mm)
best = zeros(Int8, mm)
for row in 1:mm
    sest[row, :], best[row] = decoderow(chirp, row)
end

display(sest)
display(best)
display(construct_pauli([0,0,0],[1,1,1]))

end # module
