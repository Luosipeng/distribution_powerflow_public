mutable struct Sd
    z::Vector{ComplexF64}
    i::Vector{ComplexF64}
    p::Vector{ComplexF64}
end

function makeSdzip(baseMVA, bus, pw_1, pw_2, pw_3)
    (PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM,VA, 
    BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN, PER_CONSUMER) = idx_bus();
    qw_1 = pw_1
    qw_2 = pw_2
    qw_3 = pw_3
    z = (bus[:, PD] .* pw_3  + 1im * bus[:, QD] .* qw_3) / baseMVA
    i = (bus[:, PD] .* pw_2  + 1im * bus[:, QD] .* qw_2) / baseMVA
    p = (bus[:, PD] .* pw_1  + 1im * bus[:, QD] .* qw_1) / baseMVA
    sd = Sd(z[:,1], i[:,1], p[:,1])
    return sd
end

function makeSdzip(baseMVA, bus)
    (PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM,VA, 
    BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN, PER_CONSUMER) = idx_bus();
        pw = [1 0 0]
        qw = pw
    z = (bus[:, PD] * pw[3]  + 1im * bus[:, QD] * qw[3]) / baseMVA
    i = (bus[:, PD] * pw[2]  + 1im * bus[:, QD] * qw[2]) / baseMVA
    p = (bus[:, PD] * pw[1]  + 1im * bus[:, QD] * qw[1]) / baseMVA
    sd = Sd(z, i, p)
    return sd
end