using LinearAlgebra
using Random
using Statistics
using ProgressMeter

PauliMatrices = [[0 1; 1 0], [0 -1im; 1im 0], [1 0; 0 -1]] 
Idmatrix = [1 0; 0 1] .+ 0.0im

# kron(PauliMatrices[1], PauliMatrices[1])
⊗(x, y) = kron(x, y)


############################### hamiltonian ###########################

function Ham_Jz_ops(Jz_arr::Array{Float64, 1}, N::Int64)
    term_ops_sum = zeros(ComplexF64, 2^N, 2^N)

    index = 1
    for i in 1:N
        a = 3
        term_ops = fill(Idmatrix, N)
        term_ops[i] = PauliMatrices[a]
        term_ops[mod(i, N) + 1] = PauliMatrices[a]
        term_ops_sum[:,:] += Jz_arr[i] * foldl(⊗, term_ops)

        index += 1
    end
    
    return term_ops_sum
end


function Ham_hz_ops(hz_arr::Array{Float64, 1}, N::Int64)
    term_ops_sum = zeros(ComplexF64, 2^N, 2^N)

    index = 1
    for i in 1:N
        a = 3
        term_ops = fill(Idmatrix, N)
        term_ops[i] = PauliMatrices[a]
        term_ops_sum[:,:] += hz_arr[i] * foldl(⊗, term_ops)

        index += 1
    end
    
    return term_ops_sum
end

function Ham_Jx_ops(Jx_arr::Array{Float64, 1}, N::Int64)
    term_ops_sum = zeros(ComplexF64, 2^N, 2^N)

    index = 1
    for i in 1:N
        a = 1
        term_ops = fill(Idmatrix, N)
        term_ops[i] = PauliMatrices[a]
        term_ops_sum[:,:] += Jx_arr[i] * foldl(⊗, term_ops)

        index += 1
    end
    
    return term_ops_sum
end
