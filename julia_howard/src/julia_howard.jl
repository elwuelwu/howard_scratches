module julia_howard

using DSP
using Random
using LinearAlgebra
using StaticArrays
using Hadamard
using Kronecker

const mm = 3
const N = 2^mm
const IdN = Diagonal(ones(N))


function generate_chirp(smat, bvec)
    xvec = Vector{Int64}(undef,mm)
    multvec1 = Vector{Float64}(undef, mm)
    multint1 = Float64
    multint2 = Float64
    chirp = Vector{ComplexF64}(undef,N)
    for i in 1:N
        digits!(xvec,i-1, base=2)
        mul!(multvec1,smat,xvec)
        multint1 = bvec' * xvec
        multint2 = multvec1' * xvec
        chirp[i] = im^(multint1 + 2*multint2)
        # chirp[i] = im^(xvec * smat * xvec + 2*bvec'*xvec)
    end
    return chirp
end


# Decoder

const Z = [1 0; 0 -1]
const X = [0 1; 1 0]
const I = [1 0; 0 1]
seq = [X]
for i in 1:(mm-1); push!(seq, I); end

global permss = []
for i in 1:mm
    push!(permss, kron(seq...))
    global seq = circshift(seq, 1)
end

global const perms = permss

function decoderow(vecin::Vector{ComplexF64}, row::Int64)
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

function decode(vecin)
    sest = zeros(Int8, mm, mm)
    best = zeros(Int8, mm)
    for row in 1:mm

        sest[row, :], best[row] = decoderow(vecin, row)
    end
    return sest, best
end

function simulate()
    drops = 1
    errors = 0
    # SNR in dB scale
    SNR = 10
    smat = Matrix{Int8}
    bvec = Vector{Int8}
    for i in 1:drops
        smat = rand((0,1), mm, mm)
        smat = (smat + transpose(smat)) .% 2 + Diagonal(rand((0,1), mm))
        bvec = rand((0,1), mm)  
        tstvec = generate_chirp(smat, bvec)
        # Add noise
        sest, best = decode(tstvec)
        if (smat != sest) | (bvec != best)
            errors += 1
        end
    end
    #println(errors)
end

simulate()

end # module
