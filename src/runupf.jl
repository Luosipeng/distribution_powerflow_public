"""jpc
This function is used to run a local unbalanced power flow analysis on unsymmetrical load nodes
"""
#Step1: find the unbalanced nodes in the system
#Step2: find the interface branchs of the unbalanced nodes and balanced nodes

function runupf(net ,opt)
    # 呼叫索引函数
     (GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, PC1, PC2, QC1MIN,
     QC1MAX, QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF, PW_LINEAR, POLYNOMIAL,
      MODEL, STARTUP, SHUTDOWN, NCOST, COST, MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, GEN_AREA)= PowerFlow.idx_gen()
      (PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM,VA, 
    BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN, PER_CONSUMER) = PowerFlow.idx_bus()
    # 数据处理
    jpc1 = branch_ph_2_3ph(net, opt, 1)
    jpc2, gs_eg, bs_eg = branch_ph_2_3ph(net, opt, 2)
    jpc0, _, _ = branch_ph_2_3ph(net, opt, 0)

    _, bus0, branch0, gen0, _, _, _, _, _, _ = get_pf_variables_from_jpc(jpc0)
    baseMVA, bus1, branch1, gen1, ref, pv, pq, _, _, _ = get_pf_variables_from_jpc(jpc1)
    _, bus2, branch2, gen2, _, _, _, _, _, _ = get_pf_variables_from_jpc(jpc2)

    s_abc_delta, s_abc = load_mapping(net, jpc1)
    jpc0, jpc1, jpc2, y_0_pu, y_1_pu, y_2_pu, y_0_f, y_1_f, _, y_0_t, y_1_t, _ = get_y_bus(jpc0, jpc1, jpc2)

    nb = size(bus1, 1) # number of buses in the system
    v_012_it = vcat(
    [reshape(Vector{ComplexF64}(ppc["bus"][:, VM] .* exp.(im .* deg2rad.(ppc["bus"][:, VA]))), 1, nb)
     for ppc in (jpc0, jpc1, jpc2)]...)
    #Initialize
    # For Delta transformation:
    # Voltage changed from line-earth to line-line using V_T
    # s_abc/v_abc will now give line-line currents. This is converted to line-earth
    # current using I-T
    v_del_xfmn = [1 -1 0;
                 0 1 -1;
                 -1 0 1]

    i_del_xfmn = [1 0 -1;
                 -1 1 0;
                 0 -1 1]
    v_abc_it = sequence_to_phase(v_012_it)
    
    ## Iteration using Power mismatch criterion
    outer_tolerance_mva = 3e-8
    Count = 0

    # 初始化不匹配量为布尔数组
    s_mismatch = fill(true, (2,1))

    # 记录开始时间
    t0 = time()

    # 迭代循环
    max_iteration = 30
    while any(s_mismatch .> outer_tolerance_mva) && Count < max_iteration
        # 计算标幺值功率
        s_abc_pu = -s_abc ./ jpc1["baseMVA"]
        s_abc_delta_pu = -s_abc_delta ./ jpc1["baseMVA"]
        # 计算星形连接的电流
        i_abc_it_wye = conj.(s_abc_pu ./ v_abc_it)

        # 计算三角形连接的电流
        i_abc_it_delta = i_del_xfmn * conj.(s_abc_delta_pu ./ (v_del_xfmn * v_abc_it))

        # For buses with both delta and wye loads we need to sum of their currents
        # to sum up the currents
        i_abc_it = i_abc_it_wye + i_abc_it_delta
        i012_it = phase_to_sequence(i_abc_it)
        v1_for_s1 = v_012_it[2, :]
        i1_for_s1 = -i012_it[2, :]
        v0_pu_it = transpose(v_012_it[1,:])
        v2_pu_it = transpose(v_012_it[3,:])
        i0_pu_it = transpose(i012_it[1,:])
        i2_pu_it = transpose(i012_it[3,:])
        s1 = v1_for_s1 .* conj.(i1_for_s1)

        # Current used to find S1 Positive sequence power
        jpc1["bus"][pq, PD] = real(s1[pq]) * jpc1["baseMVA"]
        jpc1["bus"][pq, QD] = imag(s1[pq]) * jpc1["baseMVA"]

        #run newton raphson power flow
        jpc1, success, iterations = run_newton_raphson_pf(jpc1,opt)
        Ybus,Yf,Yt = jpc1["internal"][:Ybus],jpc1["internal"][:Yf],jpc1["internal"][:Yt]
        V = jpc1["internal"][:V]
        bus,branch,gen = jpc1["internal"][:bus],jpc1["internal"][:branch],jpc1["internal"][:gen]
        bus, gen, branch = PowerFlow.pfsoln(baseMVA, bus, gen, branch, Ybus, Yf, Yt, V, ref, pv, pq)
        #return result
        jpc1["bus"] = bus
        jpc1["gen"] = gen
        jpc1["branch"] = branch

        # Conduct Negative and Zero sequence power flow
        v0_pu_it = V_from_I(y_0_pu, i0_pu_it)
        v2_pu_it = V_from_I(y_2_pu, i2_pu_it)

        # Evaluate Positive Sequence Power Mismatch
        i1_from_v_it = vec(I1_from_V012(v_012_it, y_1_pu))
        s_from_voltage = S_from_VI_elementwise(v1_for_s1, i1_from_v_it)
        v1_pu_it = V1_from_jpc(jpc1)
        v_012_new = combine_X012(v0_pu_it, v1_pu_it, v2_pu_it)
        s_mismatch = abs.(abs.(s1[pq]) .- abs.(s_from_voltage[pq]))
        v_012_it = v_012_new
        v_abc_it = sequence_to_phase(v_012_it)
        Count += 1
    end

    success = (Count < max_iteration)

    for ppc in [jpc0, jpc1, jpc2]
        ppc["success"] = success
    end
    # TODO: Add reference to paper to explain the following steps
    # This is required since the ext_grid power results are not correct if its
    # not done
    ref, pv, pq = PowerFlow.bustypes(jpc0["bus"], jpc0["gen"])
    jpc0["bus"][ref, GS] -= gs_eg
    jpc0["bus"][ref, BS] -= bs_eg
    y_0_pu, y_0_f, y_0_t = PowerFlow.makeYbus(jpc0["baseMVA"], jpc0["bus"], jpc0["branch"])
    # revert the change, otherwise repeated calculation with recycled elements will fail
    jpc0["bus"][ref, GS] += gs_eg
    jpc0["bus"][ref, BS] += bs_eg
    # Bus, Branch, and Gen  power values
    bus0, gen0, branch0 = PowerFlow.pfsoln(baseMVA, bus0, gen0, branch0, y_0_pu, y_0_f, y_0_t,
                                vec(v_012_it[1, :]), ref, pv, pq)
    bus1, gen1, branch1 = PowerFlow.pfsoln(baseMVA, bus1, gen1, branch1, y_1_pu, y_1_f, y_1_t,
                                vec(v_012_it[2, :]), ref, pv, pq)
    bus2, gen2, branch2 = PowerFlow.pfsoln(baseMVA, bus2, gen2, branch2, y_1_pu, y_1_f, y_1_t,
                                vec(v_012_it[3, :]) , ref, pv, pq)

    #store the results in the jpc structure
    jpc0["bus"] = bus0
    jpc0["gen"] = gen0
    jpc0["branch"] = branch0

    jpc1["bus"] = bus1
    jpc1["gen"] = gen1
    jpc1["branch"] = branch1

    jpc2["bus"] = bus2
    jpc2["gen"] = gen2
    jpc2["branch"] = branch2

    i_012_res = current_from_voltage_results(y_0_pu, y_1_pu, v_012_it)
    s_012_res = S_from_VI_elementwise(v_012_it, i_012_res) * jpc1["baseMVA"]
    eg_idx = net["ext_grid"].bus
    eg_bus_idx_ppc = Int64.(real(jpc1["gen"][eg_idx, GEN_BUS]))

    jpc0["gen"][eg_idx, PG] = real.(s_012_res[1, eg_bus_idx_ppc])
    jpc1["gen"][eg_idx, PG] = real.(s_012_res[2, eg_bus_idx_ppc])
    jpc2["gen"][eg_idx, PG] = real.(s_012_res[3, eg_bus_idx_ppc])
    jpc0["gen"][eg_idx, QG] = imag.(s_012_res[1, eg_bus_idx_ppc])
    jpc1["gen"][eg_idx, QG] = imag.(s_012_res[2, eg_bus_idx_ppc])
    jpc2["gen"][eg_idx, QG] = imag.(s_012_res[3, eg_bus_idx_ppc])

    net["jpc0"] = jpc0
    net["jpc1"] = jpc1
    net["jpc2"] = jpc2

    get_bus_v_results_3ph(net, jpc0, jpc1, jpc2)
    bus_pq = get_p_q_results_3ph(net)

    get_full_branch_zero(net,jpc0)
    get_branch_results_3ph(net, jpc0, jpc1, jpc2, bus_pq)
    get_gen_results_3ph(net, jpc0, jpc1, jpc2, bus_pq)
    get_bus_results_3ph(net, bus_pq)

    return net
end

#TODO:
# function _check_line_dc_at_b2b_buses(jpc::Dict{String, Any})

#     # 检查直流线路是否连接到背靠背VSC转换器配置的母线
#     b2b_buses = jpc["busDC"][jpc["busAC"][:, DC_BUS_TYPE] .== DC_B2B, DC_BUS_I]
#     b2b_buses = convert.(Int64, b2b_buses)
    
#     intersect_from = intersect(convert.(Int64, jpc["branch_dc"][:, DC_F_BUS]), b2b_buses)
#     intersect_to = intersect(convert.(Int64, jpc["branch_dc"][:, DC_T_BUS]), b2b_buses)
    
#     if length(intersect_from) != 0 || length(intersect_to) != 0
#         throw(NotImplementedError("发现直流线路连接到背靠背VSC转换器配置 - 未实现。直流线路只能连接到不属于背靠背配置的直流母线。"))
#     end
# end

#TODO:
# function _check_vsc_different_ac_control_modes_at_same_bus(jpc)
#     # 检查是否存在共享同一交流母线但具有不同交流控制模式的VSC转换器
#     ac_vm_pu_buses = jpc["vsc"][jpc["vsc"][:, VSC_MODE_AC] .== VSC_MODE_AC_V, VSC_BUS]
#     ac_q_mvar_buses = jpc["vsc"][jpc["vsc"][:, VSC_MODE_AC] .== VSC_MODE_AC_Q, VSC_BUS]
#     ac_slack_buses = jpc["vsc"][jpc["vsc"][:, VSC_MODE_AC] .== VSC_MODE_AC_SL, VSC_BUS]
    
#     ac_bus_intersection = vcat([intersect(a, b) for (a, b) in combinations([ac_vm_pu_buses, ac_q_mvar_buses, ac_slack_buses], 2)]...)
    
#     if length(ac_bus_intersection) != 0
#         throw(NotImplementedError("发现多个共享同一交流母线但具有不同交流控制模式的VSC转换器 - 未实现。如果VSC转换器共享同一交流母线，它们只能具有相同的交流控制模式。"))
#     end
# end

function branch_ph_2_3ph(net::Dict{Any, Any}, opt::Dict{String, Dict{String}}, sequence::Int64)
    
    # 创建新的字典，使用深拷贝避免互相影响
    jpc_new = Dict{String, Any}()
    jpc_new["baseMVA"] = deepcopy(opt["PF"]["baseMVA"])

    
    # 处理母线数据
    slack_bus = net["ext_grid"].bus
    jpc_new = PowerFlow.calculate_bus(net, jpc_new, sequence, slack_bus, opt)
    
    # 处理支路数据
    if sequence == 0
        jpc_new = PowerFlow.build_branch_Jpc_zero(net, jpc_new, opt)
    else
        jpc_new = PowerFlow.calculate_line_parameter(net, jpc_new, sequence, opt)
        jpc_new = PowerFlow.calculate_transformer_parameter(net, jpc_new, opt)
    end
    
    # 处理发电机数据
    jpc_new = PowerFlow.build_gen(net, jpc_new)
    
    # 对于负序和零序网络，添加外部电网短路阻抗
    if sequence == 2 || sequence == 0
        gs_eg, bs_eg = PowerFlow.add_grid_external_sc_impedance(jpc_new, net["ext_grid"])
        if sequence == 0
            net["_jpc0"] = deepcopy(jpc_new)
        end
        robust_process(net, jpc_new)
        
        return jpc_new, gs_eg, bs_eg
    else
        robust_process(net, jpc_new)
        return jpc_new
    end
end

function sum_by_group(bus, first_val, second_val)
    # 创建排序索引
    order = sortperm(bus)
    
    # 按排序索引重新排列数组
    sorted_bus = bus[order]
    sorted_first_val = first_val[order]
    sorted_second_val = second_val[order]
    
    # 找出唯一元素的位置
    n = length(sorted_bus)
    index = trues(n)
    for i in 1:(n-1)
        index[i] = sorted_bus[i+1] != sorted_bus[i]
    end
    
    # 计算第一组值的累积和
    cumsum_first = cumsum(sorted_first_val)
    
    # 计算第二组值的累积和
    cumsum_second = cumsum(sorted_second_val)
    
    # 提取唯一母线及对应的累积和
    unique_buses = sorted_bus[index]
    cumsum_first_unique = cumsum_first[index]
    cumsum_second_unique = cumsum_second[index]
    
    # 计算每组的和
    result_first = similar(cumsum_first_unique)
    result_second = similar(cumsum_second_unique)
    
    # 第一个元素保持不变
    result_first[1] = cumsum_first_unique[1]
    result_second[1] = cumsum_second_unique[1]
    
    # 计算差值得到每组的和
    for i in 2:length(unique_buses)
        result_first[i] = cumsum_first_unique[i] - cumsum_first_unique[i-1]
        result_second[i] = cumsum_second_unique[i] - cumsum_second_unique[i-1]
    end
    
    return unique_buses, result_first, result_second
end

function get_pf_variables_from_jpc(jpc::Dict{String,Any})
    (GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, PC1,
         PC2, QC1MIN, QC1MAX, QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, 
         RAMP_Q, APF, PW_LINEAR, POLYNOMIAL, MODEL, STARTUP, SHUTDOWN, NCOST,
          COST, MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN) = PowerFlow.idx_gen();
          (PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM,VA, 
          BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN, PER_CONSUMER) = PowerFlow.idx_bus();
    jpc_local = deepcopy(jpc)
    if isnothing(jpc_local)
        throw(ArgumentError("jpc cannot be nothing"))
    end

    #Get data from calc
    baseMVA = jpc_local["baseMVA"]
    bus = jpc_local["bus"]
    branch = jpc_local["branch"]
    gen = jpc_local["gen"]

    #get bus index list of each bus type
    ref, pv, pq = PowerFlow.bustypes(bus, gen)

    ## generator info
    on = Int.(findall(gen[:, GEN_STATUS] .> 0))  ## which generators are on?
    gbus = Int.(gen[on, GEN_BUS])  ## what buses are they at?

    ## initial state
    # V0    = ones(bus.shape[0])            ## flat start
    V0  = bus[:, VM] .* exp.(1im * pi/180 * bus[:, VA])
    V0[gbus] = gen[on, VG] ./ abs.(V0[gbus]) .* V0[gbus]

    return baseMVA, bus, branch, gen, ref, pv, pq, on, gbus, V0
end

function load_mapping(net::Dict{Any, Any}, jpc1::Dict{String, Any})

    (PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, VA,
     BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN) = PowerFlow.idx_bus()

    params = Dict{String, Any}()    
    phases = ["a", "b", "c"]
    load_types = ["wye", "delta"]
    load_elements =["load", "asymmetric_load", "gen","asymmetric_sgen"]

    for phase in phases
        for typ in load_types
            params["S$(phase)$(typ)"] = (jpc1["bus"][:, PD] .+
                                        jpc1["bus"][:, QD] .* 1im) .* 0
            params["p$(phase)$(typ)"] = Float64[]  # p values from loads/sgens
            params["q$(phase)$(typ)"] = Float64[]  # q values from loads/sgens
            params["P$(phase)$(typ)"] = Float64[]  # Aggregated Active Power
            params["Q$(phase)$(typ)"] = Float64[]  # Aggregated reactive Power
            params["b$(phase)$(typ)"] = Int64[]    # bus map for phases
            params["b$(typ)"] = Int64[]            # aggregated bus map(s_abc)
            
            for element in load_elements
                get_elements(params, net, element, phase, typ)
            end
            
            # Mapping constant power loads to buses
            if !isempty(params["b$(typ)"])
                params["b$(phase)$(typ)"] = params["b$(typ)"]
                
                params["b$(phase)$(typ)"], params["P$(phase)$(typ)"],
                temp_q = sum_by_group(params["b$(phase)$(typ)"],
                                      params["p$(phase)$(typ)"],
                                      params["q$(phase)$(typ)"] .* 1im)
                
                params["Q$(phase)$(typ)"] = imag.(temp_q)
                
                params["S$(phase)$(typ)"][params["b$(phase)$(typ)"]] .= 
                    params["P$(phase)$(typ)"] .+ params["Q$(phase)$(typ)"] .* 1im
            end
        end
    end
    Sabc_del = vcat(transpose(params["Sadelta"]), transpose(params["Sbdelta"]), transpose(params["Scdelta"]))
    Sabc_wye = vcat(transpose(params["Sawye"]), transpose(params["Sbwye"]), transpose(params["Scwye"]))
    return Sabc_del, Sabc_wye
end

function get_elements(params, net, element, phase, typ)
    """
    This function is used to get the elements of the load mapping
    Automatically skips elements that don't exist in net
    """
    # 首先检查元素是否存在于net中
    if !haskey(net, element)
        # 如果元素不存在，直接返回原参数，不做任何修改
        return params
    end

    sign = endswith(element, "sgen") ? -1 : 1
    elm = net[element]
    
    # 选择符合条件的行：active wye or active delta
    active = (elm.in_service .== true) .& (elm.type.== typ)

    asyload_p_a_mw= columnindex(elm, :p_a_mw)
    asyload_q_a_mvar= columnindex(elm, :q_a_mvar)
    asyload_p_b_mw= columnindex(elm, :p_b_mw)
    asyload_q_b_mvar= columnindex(elm, :q_b_mvar)
    asyload_p_c_mw= columnindex(elm, :p_c_mw)
    asyload_q_c_mvar= columnindex(elm, :q_c_mvar)
    
    if size(elm, 1) > 0
        if element == "load" || element == "sgen"
            vl = elm[active, asyload_scailing]
            
            # 这里需要定义p_mw_col和q_mvar_col，可能需要从asymmetric_ld_idx获取
            # 假设这些变量已经定义或从其他函数获取
            p_mw_col = asyload_p_a_mw  # 这里可能需要调整为正确的索引
            q_mvar_col = asyload_q_a_mvar  # 这里可能需要调整为正确的索引
            
            params["p$(phase)$(typ)"] = vcat(params["p$(phase)$(typ)"],
                                           elm[active, p_mw_col] ./ 3 .* vl .* sign)
            params["q$(phase)$(typ)"] = vcat(params["q$(phase)$(typ)"],
                                           (elm[active, q_mvar_col] ./ 3) .* vl .* sign)
            params["b$(typ)"] = vcat(params["b$(typ)"],
                                    convert.(Int64, elm[active, asyload_connectedbus]))
                                    
        elseif startswith(element, "asymmetric")
            vl = elm[active,:].scaling
            
            p = Dict(
                "a" => asyload_p_a_mw,
                "b" => asyload_p_b_mw,
                "c" => asyload_p_c_mw
            )
            
            q = Dict(
                "a" => asyload_q_a_mvar,
                "b" => asyload_q_b_mvar,
                "c" => asyload_q_c_mvar
            )
            
            params["p$(phase)$(typ)"] = vcat(params["p$(phase)$(typ)"],
                                           elm[active, p[phase]] .* vl .* sign)
            params["q$(phase)$(typ)"] = vcat(params["q$(phase)$(typ)"],
                                           elm[active, q[phase]] .* vl .* sign)
            params["b$(typ)"] = vcat(params["b$(typ)"],
                                    convert.(Int64, elm[active, :].bus))
        end
    end
    
    return params
end

function get_y_bus(jpc0, jpc1, jpc2)
        # build admittance matrices
        y_0_bus, y_0_f, y_0_t = PowerFlow.makeYbus(jpc0["baseMVA"], jpc0["bus"], jpc0["branch"])
        y_1_bus, y_1_f, y_1_t = PowerFlow.makeYbus(jpc1["baseMVA"], jpc1["bus"], jpc1["branch"])
        y_2_bus, y_2_f, y_2_t = PowerFlow.makeYbus(jpc2["baseMVA"], jpc2["bus"], jpc2["branch"])

    return jpc0, jpc1, jpc2, y_0_bus, y_1_bus, y_2_bus, y_0_f, y_1_f, y_2_f, y_0_t, y_1_t, y_2_t
end

function phase_shift_unit_operator(angle_deg)
    return 1 * exp(im * deg2rad(angle_deg))
end

function sequence_to_phase(X012)
    a = phase_shift_unit_operator(120)
    asq = phase_shift_unit_operator(-120)

    Tabc = [1 1 1;
            1 asq a;
            1 a asq]
    return Tabc * X012
end

function phase_to_sequence(Xabc)
    a = phase_shift_unit_operator(120)
    asq = phase_shift_unit_operator(-120)  
    T012 = [1 1 1;
        1 a asq;
        1 asq a] ./ 3 
    return T012 * Xabc
end

function store_internal(jpc, internal_storage)
    # 检查并创建internal键
    if !haskey(jpc, "internal")
        jpc["internal"] = Dict()
    end
    
    # internal storage is a dict with the variables to store in jpc["internal"]
    for (key, val) in pairs(internal_storage)
        jpc["internal"][key] = val
    end
    return jpc
end

function run_newton_raphson_pf(jpc1::Dict{String, Any}, opt)
    baseMVA, bus, branch,gen, ref, pv, pq, _, _, V0 = get_pf_variables_from_jpc(jpc1)
    # jpc1 = PowerFlow.runpf(jpc1,opt)
    Ybus, Yf, Yt = PowerFlow.makeYbus(baseMVA, bus, branch)

    Sbus = PowerFlow.makeSbus(baseMVA, bus, gen, V0)
    V, success, iterations,norm_history = PowerFlow.newtonpf(baseMVA, bus, gen, Ybus, V0, ref, pv, pq, opt["PF"]["PF_TOL"], opt["PF"]["PF_MAX_IT"], opt["PF"]["NR_ALG"])
    # keep "internal" variables in  memory / jpc["_ppc"]["internal"] -> needed for recycle.
    jpc1 = store_internal(jpc1, Dict(
    :bus => bus,
    :gen => gen,
    :branch => branch,
    :baseMVA => baseMVA,
    :V => V,
    :pv => pv,
    :pq => pq,
    :ref => ref,
    :Sbus => Sbus,
    :Ybus => Ybus,
    :Yf => Yf,
    :Yt => Yt,
))
    return jpc1, success, iterations
end

function V_from_I(Y, I)
    # 确保 I 是列向量形式
    I_vec = vec(collect(I))  # 将任何形式的 I 转换为列向量
    
    # 求解线性方程组
    V = Y \ I_vec
    
    # 返回结果
    return transpose(V)
end

function I1_from_V012(V012, Y)
    # 从对称分量中提取正序分量
    V1 = reshape(transpose(V012[2, :]), :, 1)  # 转换为列向量
    
    # 检查Y是否为稀疏矩阵
    if issparse(Y)
        # 如果是稀疏矩阵，先转为密集矩阵再计算
        i1 = Array(Matrix(Y) * V1)
        return transpose(i1)
    else
        # 如果是密集矩阵，直接计算
        i1 = Array(Y * V1)
        return transpose(i1)
    end
end

function I0_from_V012(V012, Y)
    # 从对称分量中提取零序分量
    V0 = reshape(transpose(V012[1, :]), :, 1)  # 转换为列向量
    
    # 检查Y是否为稀疏矩阵
    if issparse(Y)
        # 如果是稀疏矩阵，先转为密集矩阵再计算
        i0 = Array(Matrix(Y) * V0)
        return transpose(i0)
    else
        # 如果是密集矩阵，直接计算
        i0 = Array(Y * V0)
        return transpose(i0)
    end
end

function I2_from_V012(V012, Y)
    # 从对称分量中提取负序分量
    V2 = reshape(transpose(V012[3, :]), :, 1)  # 转换为列向量
    
    # 检查Y是否为稀疏矩阵
    if issparse(Y)
        # 如果是稀疏矩阵，先转为密集矩阵再计算
        i2 = Array(Matrix(Y) * V2)
        return transpose(i2)
    else
        # 如果是密集矩阵，直接计算
        i2 = Array(Y * V2)
        return transpose(i2)
    end
end

function S_from_VI_elementwise(v, i)
    return v .* conj(i)
end

function V1_from_jpc(jpc)
    (PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, VA,
     BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN) = PowerFlow.idx_bus()
    return transpose(
        jpc["bus"][:, VM] .* exp.(1im .* deg2rad.(jpc["bus"][:, VA]))
    )
end

function combine_X012(X0, X1, X2)
    return vcat(X0, X1, X2)
end

function current_from_voltage_results(y_0_pu, y_1_pu, v_012_pu)
    I012_pu = combine_X012(I0_from_V012(v_012_pu, y_0_pu),
                          I1_from_V012(v_012_pu, y_1_pu),
                          I2_from_V012(v_012_pu, y_1_pu))
    return I012_pu
end

function get_bus_v_results_3ph(net, jpc0, jpc1, jpc2)
    (PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, VA,
     BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN) = PowerFlow.idx_bus()

    V012_pu = zeros(ComplexF64, 3, length(jpc1["bus"][:, BUS_I]))


    V012_pu[1, :] = jpc0["bus"][:, VM] .* exp.(1im * pi/180 * jpc0["bus"][:, VA])
    V012_pu[2, :] = jpc1["bus"][:, VM] .* exp.(1im * pi/180 * jpc1["bus"][:, VA])
    V012_pu[3, :] = jpc2["bus"][:, VM] .* exp.(1im * pi/180 * jpc2["bus"][:, VA])
    # 取消注释以获得以kV为单位的结果而非标幺值
    # bus_base_kv = ppc0["bus"][:, BASE_KV] ./ sqrt(3)
    # V012_pu = V012_pu .* bus_base_kv
    
    Vabc_pu = sequence_to_phase(V012_pu)
    
    net["res_bus_3ph"] = DataFrame()
    net["res_bus_3ph"].bus = net["bus"][:, BUS_I]
    # 计算电压幅值
    net["res_bus_3ph"].vm_a_pu = abs.(Vabc_pu[1, :])
    net["res_bus_3ph"].vm_b_pu = abs.(Vabc_pu[2, :])
    net["res_bus_3ph"].vm_c_pu = abs.(Vabc_pu[3, :])
    
    # 电压角度
    net["res_bus_3ph"].va_a_degree = angle.(Vabc_pu[1, :]) .* 180 ./ π
    net["res_bus_3ph"].va_b_degree = angle.(Vabc_pu[2, :]) .* 180 ./ π
    net["res_bus_3ph"].va_c_degree = angle.(Vabc_pu[3, :]) .* 180 ./ π
    
    net["res_bus_3ph"].unbalance_percent = abs.(V012_pu[3, :] ./ V012_pu[2, :]) .* 100
    
    
    return net
end

function get_p_q_results_3ph(net)
    (PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, VA,
     BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN) = PowerFlow.idx_bus()
    # 创建结果矩阵 (bus, p in kw, q in kvar)
    bus_pq = zeros(Float64, length(net["bus"].index), 6)
    b, pA, pB, pC, qA, qB, qC = Float64[], Float64[], Float64[], Float64[], Float64[], Float64[], Float64[]

    # Todo: Voltage dependent loads
    elements = ["storage", "sgen", "load"]
    elements_3ph = ["asymmetric_load", "asymmetric_sgen"]
    
    for element in elements
        # 检查元素是否存在于jpc中
        if !haskey(net, element)
            continue  # 如果不存在，跳过此次循环
        end
        
        sign = element in ["sgen", "asymmetric_sgen"] ? -1 : 1
        if length(net[element].bus) > 0
            write_pq_results_to_element(net, net._ppc1, element, suffix="3ph")
            p_el, q_el, bus_el = get_p_q_b(net, element, suffix="3ph")
            pA = vcat(pA, sign .* p_el ./ 3)
            pB = vcat(pB, sign .* p_el ./ 3)
            pC = vcat(pC, sign .* p_el ./ 3)
            qA = vcat(qA, ac ? sign .* q_el ./ 3 : zeros(length(p_el ./ 3)))
            qB = vcat(qB, ac ? sign .* q_el ./ 3 : zeros(length(p_el ./ 3)))
            qC = vcat(qC, ac ? sign .* q_el ./ 3 : zeros(length(p_el ./ 3)))
            b = vcat(b, bus_el)
        end
    end
    
    for element in elements_3ph
        if !haskey(net, element)
            continue  # 如果不存在，跳过此次循环
        end
        sign = element in ["sgen", "asymmetric_sgen"] ? -1 : 1
        if length(net[element].bus) > 0
            write_pq_results_to_element_3ph(net, element)
            p_el_A, q_el_A, p_el_B, q_el_B, p_el_C, q_el_C, bus_el = get_p_q_b_3ph(net, element)
            pA = vcat(pA, sign .* p_el_A)
            pB = vcat(pB, sign .* p_el_B)
            pC = vcat(pC, sign .* p_el_C)
            qA = vcat(qA, sign .* q_el_A )
            qB = vcat(qB, sign .* q_el_B )
            qC = vcat(qC, sign .* q_el_C )
            b = vcat(b, bus_el)
        end
    end

    # 按组汇总每个元素的pq结果，稍后写入jpc['bus']
    b_pp, vp_A, vq_A, vp_B, vq_B, vp_C, vq_C = _sum_by_group_nvals(Int64.(b), pA, qA, pB, qB, pC, qC)
    bus_pq[b_pp, 1] = vp_A
    bus_pq[b_pp, 2] = vq_A
    bus_pq[b_pp, 3] = vp_B
    bus_pq[b_pp, 4] = vq_B
    bus_pq[b_pp, 5] = vp_C
    bus_pq[b_pp, 6] = vq_C
    
    return bus_pq
end

function write_pq_results_to_element_3ph(net, element)
    """
    get p_mw and q_mvar for a specific pq element ("load", "sgen"...).
    This function basically writes values element table to res_element table

    :param jpc: pandapower jpc
    :param element: element name (str)
    :return: jpc with updated results
    """  
    # info element
    el_data = net[element]
    res_ = "res_" * element * "_3ph"
    
    # 如果结果DataFrame不存在，创建它
    if !haskey(net, res_)
        net[res_] = DataFrame()
    end

    scaling = el_data.scaling
    element_in_service = el_data.in_service .== true

    # 假设p_mw和q_mvar的列索引
    p_mw_idx = columnindex(el_data, :p_a_mw) # 假设这是总有功功率的列索引
    q_mvar_idx = columnindex(el_data, :q_a_mvar) # 假设这是总无功功率的列索引

    # 创建结果数组
    if element in ["load", "sgen"]
        # 对于load和sgen，将总功率平均分配到三相
        p_a_mw = (el_data[:, p_mw_idx] ./ 3) .* scaling .* element_in_service
        p_b_mw = (el_data[:, p_mw_idx] ./ 3) .* scaling .* element_in_service
        p_c_mw = (el_data[:, p_mw_idx] ./ 3) .* scaling .* element_in_service
    else
        # 对于其他元素，使用各相的具体功率
        p_a_mw = el_data.p_a_mw .* scaling .* element_in_service
        p_b_mw = el_data.p_b_mw.* scaling .* element_in_service
        p_c_mw = el_data.p_c_mw .* scaling .* element_in_service
    end

    # 将结果存入jpc中
    net[res_].p_a_mw = p_a_mw
    net[res_].p_b_mw = p_b_mw
    net[res_].p_c_mw = p_c_mw

    
    # 计算无功功率结果
    if element in ["load", "sgen"]
        q_a_mvar = (el_data[:, q_mvar_idx] ./ 3) .* scaling .* element_in_service
        q_b_mvar = (el_data[:, q_mvar_idx] ./ 3) .* scaling .* element_in_service
        q_c_mvar = (el_data[:, q_mvar_idx] ./ 3) .* scaling .* element_in_service
    else
        q_a_mvar = el_data.q_a_mvar.* scaling .* element_in_service
        q_b_mvar = el_data.q_b_mvar .* scaling .* element_in_service
        q_c_mvar = el_data.q_c_mvar.* scaling .* element_in_service
    end
    
    # 将无功功率结果存入jpc中
    net[res_].q_a_mvar = q_a_mvar
    net[res_].q_b_mvar = q_b_mvar
    net[res_].q_c_mvar = q_c_mvar


    return net
end

function get_p_q_b_3ph(net, element)

    res_ = "res_" * element * "_3ph"

    # bus values are needed for stacking
    b = net[element].bus # 在 Julia 中，我们不需要 .values
    pA = net[res_].p_a_mw
    pB = net[res_].p_b_mw
    pC = net[res_].p_c_mw
    
    # 在 Julia 中使用三元运算符 ? :
    qA = net[res_].q_a_mvar 
    qB = net[res_].q_b_mvar 
    qC = net[res_].q_c_mvar
    
    return pA, qA, pB, qB, pC, qC, b
end

function _sum_by_group_nvals(bus, vals...)
    # 获取排序顺序
    order = sortperm(bus)
    sorted_bus = bus[order]
    
    # 创建索引数组，标记每个唯一bus的开始位置
    index = ones(Bool, length(sorted_bus))
    if length(sorted_bus) > 1
        index[2:end] = sorted_bus[2:end] .!= sorted_bus[1:end-1]
    end
    
    # 提取唯一的bus值
    unique_bus = sorted_bus[index]
    
    # 创建结果数组
    newvals = ntuple(i -> zeros(eltype(vals[i]), length(unique_bus)), length(vals))
    
    # 对每个值进行分组求和
    for (i, val) in enumerate(vals)
        sorted_val = val[order]
        # 计算累积和
        cumulative_sum = cumsum(sorted_val)
        
        # 提取每个组的累积和
        group_sums = cumulative_sum[index]
        
        # 计算每个组的和（通过差分）
        if length(group_sums) > 1
            group_sums[2:end] = group_sums[2:end] .- group_sums[1:end-1]
        end
        
        # 存储结果
        newvals[i] .= group_sums
    end
    
    # 返回结果元组，首个元素是唯一的bus值
    return (unique_bus, newvals...)
end

function get_branch_results_3ph(net, jpc0, jpc1, jpc2, pq_buses)
    """
    Extract the bus results and writes it in the Dataframe net.res_line and net.res_trafo.

    INPUT:

        **results** - the result of runpf loadflow calculation

        **p** - the dict to dump the "res_line" and "res_trafo" Dataframe

    """
    I012_f, S012_f, V012_f, I012_t, S012_t, V012_t = get_branch_flows_3ph(jpc0, jpc1, jpc2)
    # get_line_results_3ph(net, jpc0, jpc1, jpc2, I012_f, V012_f, I012_t, V012_t)
    # get_trafo_results_3ph(net, jpc0, jpc1, jpc2, I012_f, V012_f, I012_t, V012_t)
    # _get_trafo3w_results(net, jpc, s_ft, i_ft)
    # _get_impedance_results(net, jpc, i_ft)
    # _get_xward_branch_results(net, jpc, bus_lookup_aranged, pq_buses)
    # _get_switch_results(net, i_ft)
end

function get_branch_flows_3ph(jpc0, jpc1, jpc2)
    (F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, TAP,
     SHIFT, BR_STATUS, PF, QF, PT, QT) = PowerFlow.idx_brch()
    (PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, VA,
     BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN) = PowerFlow.idx_bus()

    br_from_idx = Int64.(real(jpc1["branch"][:, F_BUS]))
    br_to_idx = Int64.(real(jpc1["branch"][:, T_BUS]))
    
    # 计算起始端电压 - 使用 vec 代替 flatten
    V012_f = [vec(jpc["bus"][br_from_idx, VM] .* jpc["bus"][br_from_idx, BASE_KV] .*
                 exp.(1im .* deg2rad.(jpc["bus"][br_from_idx, VA]))) for jpc in [jpc0, jpc1, jpc2]]
    V012_f = transpose(hcat(V012_f...))  # 转换为与Python等效的数组形状
    
    # 计算终止端电压
    V012_t = [vec(jpc["bus"][br_to_idx, VM] .* jpc["bus"][br_to_idx, BASE_KV] .*
                 exp.(1im .* deg2rad.(jpc["bus"][br_to_idx, VA]))) for jpc in [jpc0, jpc1, jpc2]]
    V012_t = transpose(hcat(V012_t...))  # 转换为与Python等效的数组形状
    
    # 计算起始端复功率
    S012_f = [real(jpc["branch"][:, PF]) .+ 1im .* real(jpc["branch"][:, QF]) for jpc in [jpc0, jpc1, jpc2]]
    S012_f = transpose(hcat(S012_f...))  # 转换为与Python等效的数组形状
    
    # 计算终止端复功率
    S012_t = [real(jpc["branch"][:, PT]) .+ 1im .* real(jpc["branch"][:, QT]) for jpc in [jpc0, jpc1, jpc2]]
    S012_t = transpose(hcat(S012_t...))  # 转换为与Python等效的数组形状
    
    # 计算电流
    I012_f = I_from_SV_elementwise(S012_f, V012_f ./ sqrt(3))
    I012_t = I_from_SV_elementwise(S012_t, V012_t ./ sqrt(3))

    return I012_f, S012_f, V012_f, I012_t, S012_t, V012_t
end

function I_from_SV_elementwise(S, V)
    # 创建与 S 相同大小的零数组
    result = zeros(Complex{Float64}, size(S))
    
    # 对于 V 不为零的元素，计算 S/V 的共轭
    for i in eachindex(S)
        if V[i] != 0
            result[i] = conj(S[i] / V[i])
        end
    end
    
    return result
end

function get_gen_results_3ph(net, jpc0, jpc1, jpc2, pq_bus)
    (GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, PC1,
    PC2, QC1MIN, QC1MAX, QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, 
    RAMP_Q, APF, PW_LINEAR, POLYNOMIAL, MODEL, STARTUP, SHUTDOWN, NCOST,
     COST, MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN) = PowerFlow.idx_gen()
    eg_end = size(net["ext_grid"],1)
    if haskey(net, "gen")
        gen_end = eg_end + length(findall(net["gen"][:,GEN_STATUS] .== 1))
    else
        gen_end = eg_end
    end
    

    # 获取外部电网结果
    b, pA, qA, pB, qB, pC, qC = get_ext_grid_results_3ph(net, jpc0, jpc1, jpc2)

    # 获取发电机结果
    if gen_end > eg_end
        b, pA, qA, pB, qB, pC, qC = get_pp_gen_results_3ph(net, jpc0, jpc1, jpc2, b, pA, qA, pB, qB, pC, qC)
    end

    # 按组对值求和
    b_pp, pA_sum, qA_sum, pB_sum, qB_sum, pC_sum, qC_sum = _sum_by_group_nvals(
        convert.(Int64, b), pA, qA, pB, qB, pC, qC)
    
    # 更新母线功率
    b_jpc = b_pp
    pq_bus[b_jpc, 1] .-= pA_sum
    pq_bus[b_jpc, 2] .-= qA_sum
    pq_bus[b_jpc, 3] .-= pB_sum
    pq_bus[b_jpc, 4] .-= qB_sum
    pq_bus[b_jpc, 5] .-= pC_sum
    pq_bus[b_jpc, 6] .-= qC_sum
end

function get_ext_grid_results_3ph(net, jpc0, jpc1, jpc2)
    (PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM, VA,
     BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN) = PowerFlow.idx_bus()
     (GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, PC1,
     PC2, QC1MIN, QC1MAX, QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, 
     RAMP_Q, APF, PW_LINEAR, POLYNOMIAL, MODEL, STARTUP, SHUTDOWN, NCOST,
      COST, MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN) = PowerFlow.idx_gen();
    # 获取外部电网结果
    n_res_eg = size(net["ext_grid"],1)
    
    eg_idx_net = net["ext_grid"].index
    # 修复：使用发电机的母线索引而不是发电机索引
    eg_bus_idx_net = net["ext_grid"].bus
    
    # 从jpc读取这些母线的结果
    V012 = zeros(ComplexF64, (3, n_res_eg))
    V012_temp= [net["bus"][eg_bus_idx_net, VM] .* net["bus"][eg_bus_idx_net, BASE_KV] .*
                          exp.(1im .* deg2rad.(net["bus"][eg_bus_idx_net, VA]))
                          for net in [jpc0, jpc1, jpc2]]
    V012[:, eg_idx_net] = transpose(hcat(V012_temp...))  # 转换为与Python等效的数组形状
    S012 = zeros(ComplexF64, (3, n_res_eg))
    S012_temp = [(net["gen"][eg_idx_net, PG] .+ 1im .* 
                           net["gen"][eg_idx_net, QG])
                           for net in [jpc0, jpc1, jpc2]]
    S012[:, eg_idx_net] = transpose(hcat(S012_temp...))  # 转换为与Python等效的数组形状
    Sabc, Vabc = SVabc_from_SV012(S012, V012 ./ sqrt(3), n_res=n_res_eg, idx=eg_idx_net)

    pA = vec(real.(Sabc[1, :]))
    pB = vec(real.(Sabc[2, :]))
    pC = vec(real.(Sabc[3, :]))
    qA = vec(imag.(Sabc[1, :]))
    qB = vec(imag.(Sabc[2, :]))
    qC = vec(imag.(Sabc[3, :]))

    # 将结果存储在net["res"]中
    if !haskey(net, "res_ext_grid_3ph")
        net["res_ext_grid_3ph"] = DataFrame()
    end

    # 复制结果的索引
    net["res_ext_grid_3ph"].bus = net["ext_grid"].bus
    net["res_ext_grid_3ph"].p_a_mw = pA
    net["res_ext_grid_3ph"].p_b_mw = pB
    net["res_ext_grid_3ph"].p_c_mw = pC
    net["res_ext_grid_3ph"].q_a_mvar = qA
    net["res_ext_grid_3ph"].q_b_mvar = qB
    net["res_ext_grid_3ph"].q_c_mvar = qC

    # 获取pq_bus的母线值
    b = net["ext_grid"].bus
    

    return b, pA, qA, pB, qB, pC, qC
end

function SVabc_from_SV012(S012, V012; n_res=nothing, idx=nothing)
    if isnothing(n_res)
        n_res = size(S012, 2)
    end
    
    if isnothing(idx)
        idx = trues(n_res)
    end
    
    I012 = zeros(ComplexF64, (3, n_res))
    I012[:, idx] = I_from_SV_elementwise(S012[:, idx], V012[:, idx])
    
    Vabc = sequence_to_phase(V012)
    Iabc = sequence_to_phase(I012)
    Sabc = S_from_VI_elementwise(Vabc, Iabc)
    
    return Sabc, Vabc
end

function get_bus_results_3ph(net, bus_pq)
    (PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM,VA, 
    BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN, PER_CONSUMER) = PowerFlow.idx_bus();

    # update index in res bus bus
    # jpc["res_bus"].index = jpc["bus"].index
    net["res_bus_3ph"].bus = net["bus"][:, BUS_I]
    # write sum of p and q values to bus
    net["res_bus_3ph"].p_a_mw = bus_pq[:, 1]  # Julia 索引从 1 开始
    net["res_bus_3ph"].p_b_mw = bus_pq[:, 3]  # Python 的 2 对应 Julia 的 3
    net["res_bus_3ph"].p_c_mw = bus_pq[:, 5]  # Python 的 4 对应 Julia 的 5
    

    net["res_bus_3ph"].q_a_mvar = bus_pq[:, 2]  # Python 的 1 对应 Julia 的 2
    net["res_bus_3ph"].q_b_mvar = bus_pq[:, 4]  # Python 的 3 对应 Julia 的 4
    net["res_bus_3ph"].q_c_mvar = bus_pq[:, 6]  # Python 的 5 对应 Julia 的 6


    # Todo: OPF

    
    
    return nothing
end

function _clean_up(net, res=true)
    
    # 大部分代码被注释掉了，我保留了这些注释
    # mode = net.__internal_options["mode"]

    # set internal selected _is_elements to nothing. This way it is not stored (saves disk space)
    # net._is_elements = nothing

    #    mode = net["_options"]["mode"]
    #    if res
    #        res_bus = mode == "sc" ? net["res_bus_sc"] : 
    #                  mode == "pf_3ph" ? net["res_bus_3ph"] : 
    #                  net["res_bus"]
    #    end
    #    if length(net["trafo3w"]) > 0
    #        buses_3w = net["trafo3w"]["ad_bus"]
    #        deleteat!(net["bus"], buses_3w)
    #        select!(net["trafo3w"], Not(:ad_bus))
    #        if res
    #            deleteat!(res_bus, buses_3w)
    #        end
    #    end
    #
    #    if length(net["xward"]) > 0
    #        xward_buses = net["xward"]["ad_bus"]
    #        deleteat!(net["bus"], xward_buses)
    #        select!(net["xward"], Not(:ad_bus))
    #        if res
    #            deleteat!(res_bus, xward_buses)
    #        end
    #    end

    if length(net["dcline"]) > 0
        # 获取直流线路对应的发电机索引
        dc_gens = net.gen.index[(length(net.gen) - length(net.dcline) * 2):end]
        # 删除这些发电机
        net.gen = net.gen[setdiff(1:nrow(net.gen), dc_gens), :]
        if res
            net.res_gen = net.res_gen[setdiff(1:nrow(net.res_gen), dc_gens), :]
        end
    end
    
    return nothing
end

function robust_process(net, jpc)
    (PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM,VA, 
    BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN, PER_CONSUMER) = PowerFlow.idx_bus();
    (F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, TAP, SHIFT, BR_STATUS, ANGMIN,
    ANGMAX, DICTKEY, PF, QF, PT, QT, MU_SF, MU_ST, MU_ANGMIN, MU_ANGMAX, LAMBDA, SW_TIME, RP_TIME, BR_TYPE, BR_AREA) = PowerFlow.idx_brch()
         (GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, PC1,
         PC2, QC1MIN, QC1MAX, QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, 
         RAMP_Q, APF, PW_LINEAR, POLYNOMIAL, MODEL, STARTUP, SHUTDOWN, NCOST,
          COST, MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN) = PowerFlow.idx_gen();

    bus = jpc["bus"]
    branch = jpc["branch"]
    gen = jpc["gen"]

    nb = size(bus, 1)
    
    internal_index = collect(1:nb)
    # 这里的 internal_index 是一个从 1 到 nb 的整数数组，表示所有母线的索引
    # 你可以根据需要修改这个索引数组，以选择特定的母线进行处理
    bus_i = bus[:, BUS_I]

    # 创建从 internal_index 到 bus_i 的映射
    idx2bus = Dict(i => bus_i[i] for i in internal_index)

    # 创建从 bus_i 到 internal_index 的映射
    bus2idx = Dict(bus_i[i] => i for i in internal_index)

    net["bus2idx"] = bus2idx
    net["idx2bus"] = idx2bus
    #更新bus
    bus[:, BUS_I] = map(k -> bus2idx[k], bus[:, BUS_I])
    # 更新发电机和支路的母线索引
    gen[:, GEN_BUS] = map(k -> bus2idx[k], gen[:, GEN_BUS])
    branch[:, F_BUS] = map(k -> bus2idx[k], branch[:, F_BUS])
    branch[:, T_BUS] = map(k -> bus2idx[k], branch[:, T_BUS])

    # 更新外部电网的母线索引
    if haskey(net, "ext_grid")
        net["ext_grid"].bus = map(k -> bus2idx[k], net["ext_grid"].bus)
    end
    if haskey(net, "gen")
        net["gen"].bus = map(k -> bus2idx[k], net["gen"].bus)
    end
    if haskey(net, "load")
        net["load"].bus = map(k -> bus2idx[k], net["load"].bus)
    end
    if haskey(net, "sgen")
        net["sgen"].bus = map(k -> bus2idx[k], net["sgen"].bus)
    end
    if haskey(net, "storage")
        net["storage"].bus = map(k -> bus2idx[k], net["storage"].bus)
    end
    if haskey(net, "shunt")
        net["shunt"].bus = map(k -> bus2idx[k], net["shunt"].bus)
    end
    if haskey(net, "trafo")
        net["trafo"].hv_bus = map(k -> bus2idx[k], net["trafo"].hv_bus)
        net["trafo"].lv_bus = map(k -> bus2idx[k], net["trafo"].lv_bus)
    end
    if haskey(net, "trafo3w")
        net["trafo3w"].hv_bus = map(k -> bus2idx[k], net["trafo3w"].hv_bus)
        net["trafo3w"].mv_bus = map(k -> bus2idx[k], net["trafo3w"].mv_bus)
        net["trafo3w"].lv_bus = map(k -> bus2idx[k], net["trafo3w"].lv_bus)
    end
    if haskey(net, "line")
        net["line"].hv_bus = map(k -> bus2idx[k], net["line"].from_bus)
        net["line"].lv_bus = map(k -> bus2idx[k], net["line"].from_bus)
    end

    mask = branch[:, BR_STATUS] .== 1
    branch = branch[mask, :]
    jpc["bus"] = bus
    jpc["branch"] = branch
    jpc["gen"] = gen


end

function get_full_branch_zero(net,jpc0)
    (F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, TAP, SHIFT, BR_STATUS, ANGMIN,
    ANGMAX, DICTKEY, PF, QF, PT, QT, MU_SF, MU_ST, MU_ANGMIN, MU_ANGMAX, LAMBDA, SW_TIME, RP_TIME, BR_TYPE, BR_AREA) = PowerFlow.idx_brch()

    branch_not_in_service = net["_jpc0"]["branch"][findall(net["_jpc0"]["branch"][:, BR_STATUS] .== 0),:]
    branch0 = jpc0["branch"]

    lb = size(branch_not_in_service, 1)
    connected_branches = zeros(Int64, lb, 4)
    branch_not_in_service = hcat(branch_not_in_service, connected_branches)

    branch0 = vcat(branch0, branch_not_in_service)

    jpc0["branch"] = branch0
end