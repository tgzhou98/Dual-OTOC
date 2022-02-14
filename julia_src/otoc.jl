# using Distributed
using JLD
using DelimitedFiles
using Plots
using LinearAlgebra
# gr()
# using Arpack
ENV["GKSwstype"]="100"
gr()

# @everywhere begin
#     # using ProgressMeter
#     using Distributed
#     using SharedArrays

    include("kicked_Ising_ed.jl")

    function otoc_one_shot(Jz::Float64, Jx::Float64, hz::Float64, Lt::Int64, Lr::Int64, phi_arr::Array{Float64, 1})
        Jz_arr = Jz*ones(Lr)
        hz_arr = hz*ones(Lr)
        Jx_arr = Jx*ones(Lr)

        hamJz = Ham_Jz_ops(Jz_arr, Lr)
        hamhz = Ham_hz_ops(hz_arr, Lr)
        hamJx = Ham_Jx_ops(Jx_arr, Lr)
        
        UF = exp(1im*(hamJz+hamhz)) * exp(1im*(hamJx))

        # observable
        # at site 1

        op1 = zeros(ComplexF64, 2^Lr, 2^Lr)
        a = 3
        term_ops = fill(Idmatrix, Lr)
        term_ops[1] = PauliMatrices[a]
        op1 = foldl(⊗, term_ops)

        # boundary operator 

        # phi = philist[ii]
        # Z_op_sum = Ham_hz_ops(ones(Lr), Lr)
        # boundary_op = exp(i*phi*Z_op_sum)

        # Calculate different version of OTOC
        Utemp = Matrix{ComplexF64}(I, 2^Lr, 2^Lr)

        otoc_type_list = zeros(Float64, Lt+1, length(phi_arr), 4)
        saturate_val_list = zeros(Float64, length(phi_arr))
        # Type 1 Tr(e^{iHt} a e^{iHt} e^{iϕW}e^{iHt}a e^{iHt}e^{-iϕW})
        # Type 2 Tr(e^{iH(t+1)} a e^{iHt} e^{iϕW}e^{iHt}a e^{iHt}e^{-iϕW})
        # Type 3 Tr(e^{iHt} a e^{iHt} e^{iϕW}e^{iHt}a e^{iHt}e^{-iϕW})
        
        Z_op_sum = Ham_hz_ops(ones(Lr), Lr)

        for t in 1:Lt+1
            print("t=$(t), ")

            # Type 1
            u1 = Utemp*UF*op1*adjoint(Utemp)
            u2 = adjoint(UF)*Utemp*op1*adjoint(Utemp)
            for (phi_index, phi) in enumerate(phi_arr)
                # boundary operator 
                boundary_op = exp(1im*phi*Z_op_sum)
                otoc_type_list[t, phi_index, 1] += real(tr(u1*boundary_op*u2*adjoint(boundary_op) ) ) / (2^Lr)
            end

            # Type 2
            u1 = Utemp*UF*op1*adjoint(UF)*op1*adjoint(Utemp)*UF
            u2 = adjoint(UF)
            for (phi_index, phi) in enumerate(phi_arr)
                # boundary operator 
                boundary_op = exp(1im*phi*Z_op_sum)
                otoc_type_list[t, phi_index, 2] += real(tr(u1*boundary_op*u2*adjoint(boundary_op) ) ) / (2^Lr)
            end

            # Type 3 all four e^{iHt}
            u1 = Utemp*op1*Utemp
            u2 = Utemp*op1*Utemp
            for (phi_index, phi) in enumerate(phi_arr)
                # boundary operator 
                boundary_op = exp(1im*phi*Z_op_sum)
                otoc_type_list[t, phi_index, 3] += real(tr(u1*boundary_op*u2*adjoint(boundary_op) ) ) / (2^Lr)
            end

            # Type 4 normal OTOC
            u1 = Utemp*op1*adjoint(Utemp)
            u2 = Utemp*op1*adjoint(Utemp)
            for (phi_index, phi) in enumerate(phi_arr)
                # boundary operator 
                boundary_op = exp(1im*phi*Z_op_sum)
                otoc_type_list[t, phi_index, 4] += real(tr(u1*boundary_op*u2*adjoint(boundary_op) ) ) / (2^Lr)
            end

            Utemp = Utemp*UF


            for (phi_index, phi) in enumerate(phi_arr)
                # boundary operator 
                boundary_op = exp(1im*phi*Z_op_sum)
                saturate_val_list[phi_index] = real(tr(UF*boundary_op*adjoint(UF)*adjoint(boundary_op))/(2^Lr) * tr(UF*op1*adjoint(UF)*op1)/(2^Lr)  )
            end
        end


        print("")

        return otoc_type_list, saturate_val_list
    end

# end



Lr = 10
Lt = 15
phi_arr = [0.0, 0.2, 0.4, 0.6]

data_dirname = "dual_circuits_otoc_20211125"
if !isdir(data_dirname)
    mkdir(data_dirname)
end

# Jz_list = [2.0]
# hz_list = [0.9]
# Jx_list = [0.5]

Jz_list = []
hz_list = []
Jx_list = []

for i in 1:3
    for j in 1:3
        for z in 1:3
            append!(Jz_list, 1.5+0.5*i)
            append!(hz_list, 0.5+0.2*j)
            append!(Jx_list, 0.3+0.1*z)
        end
    end
end

otoc_type_arr_data = zeros(Lt+1, length(phi_arr), 4, length(Jz_list))
saturate_val_list_data = zeros(length(phi_arr), length(Jz_list))

# 
# 
# otoc_aver_listall = JLD.load("otoc_aver.jld", "otoc_aver")

t_list = collect(0:Lt)

label_phi = Array{String, 2}(undef, 1, length(phi_arr))
for ll in 1:length(phi_arr)
    label_phi[1,ll] = "\$ \\phi=$(phi_arr[ll]) \$"
end
using LaTeXStrings

for i in 1:length(Jz_list)
    Jz = Jz_list[i]
    hz = hz_list[i]
    Jx = Jx_list[i]
    println("index=$(i), J_z=$(Jz), h_z=$(hz), J_x=$(Jx)")

    otoc_type_arr_data[:,:,:,i], saturate_val_list_data[:,i] =  otoc_one_shot(Jz, Jx, hz, Lt, Lr, phi_arr)
    JLD.save("$(data_dirname)/otoc_Jz$(Jz)_hz$(hz)_Jx$(Jx)_Lr$(Lr)_Lt$(Lt).jld", "otoc_type_arr_data",otoc_type_arr_data)

    p1 = plot(t_list, otoc_type_arr_data[:,:,1,i], title = "\$ O(t+\\Delta t) W(\\Delta t) O(t) W(0) \$", xlabel="\$ L_t \$",ylabel=L"F(\phi, L_t)",    markershape = :auto, markersize = 3, label=label_phi, titlefont=9, legendfontsize=6, xguidefontsize=6,yguidefontsize=6)
    plot!(saturate_val_list_data[:,i], seriestype="hline")
    p2 = plot(t_list, otoc_type_arr_data[:,:,2,i], title = "\$ O(t+\\Delta t) O(t) W(\\Delta t) W(0) \$", xlabel="\$ L_t \$",ylabel=L"F(\phi, L_t)",    markershape = :auto, markersize = 3, label=label_phi, titlefont=9, legendfontsize=6, xguidefontsize=6,yguidefontsize=6)
    plot!(saturate_val_list_data[:,i], seriestype="hline")
    p3 = plot(t_list, otoc_type_arr_data[:,:,3,i], title = "\$ Tr[e^{iHt}O e^{iHt}W e^{iHt} O e^{iHt}W] \$", xlabel="\$ L_t \$",ylabel=L"F(\phi, L_t)", markershape = :auto, markersize = 3, label=label_phi, titlefont=9, legendfontsize=6, xguidefontsize=6,yguidefontsize=6)
    p4 = plot(t_list, otoc_type_arr_data[:,:,4,i], title = "\$ O(t) W(0) O(t) W(0) \$", xlabel="\$ L_t \$",ylabel=L"F(\phi, L_t)",                      markershape = :auto, markersize = 3, label=label_phi, titlefont=9, legendfontsize=6, xguidefontsize=6,yguidefontsize=6)
    # plot(p1,p2,p3,p4, layout=(2,2), title = "\$ J_z=$(Jz), h_z=$(hz), J_x=$(Jx) \$")
    plot(p1,p2,p3,p4, layout=(2,2))

    # JLD.save("$(data_dirname)/otoc_aver.jld", "otoc_aver", otoc_aver_listall)
    savefig("$(data_dirname)/otoc_Jz$(Jz)_hz$(hz)_Jx$(Jx)_Lr$(Lr)_Lt$(Lt).pdf")
    # savefig("thermalotoc.pdf")
end
