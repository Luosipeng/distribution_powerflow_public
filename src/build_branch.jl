function calculate_line_parameter(net, jpc, sequence, opt)
    # Call the indexing function to get the indices for the branch matrix
    (F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, TAP, SHIFT, BR_STATUS, ANGMIN,
    ANGMAX, DICTKEY, PF, QF, PT, QT, MU_SF, MU_ST, MU_ANGMIN, MU_ANGMAX, LAMBDA, SW_TIME, RP_TIME, BR_TYPE, BR_AREA) = PowerFlow.idx_brch()

    nbr = size(net["line"], 1)
    branch = zeros(nbr, 14)
    # Calculate the line parameters for the branch matrix
    if sequence ==1
        line = deepcopy(net["line"])
        branch[:, F_BUS] = line.from_bus
        branch[:, T_BUS] = line.to_bus
        # basekv = bus[line.from_bus, BASE_KV]
        basekv = net["bus"][line.from_bus,:].vn_kv
        baseR = basekv.^2 ./(3*opt["PF"]["baseMVA"])
        branch[:, BR_R] = line.length_km .* line.r_ohm_per_km ./baseR ./line.parallel
        branch[:, BR_X] = line.length_km .* line.x_ohm_per_km ./baseR ./line.parallel
        b = 2*opt["PF"]["net_hz"]*pi*line.length_km .*line.c_nf_per_km .*1e-9 .*baseR .*line.parallel
        g = line.g_us_per_km .* 1e-6 .* baseR .* line.length_km .* line.parallel
        branch[:, RATE_A] .=100
        branch[:, BR_B] = b
        branch[:,BR_STATUS] = line.in_service .== true
        branch[:,ANGMIN] .=-360
        branch[:,ANGMAX] .= 360
    elseif sequence ==2
        line = deepcopy(net["line"])
        branch[:, F_BUS] = line.from_bus
        branch[:, T_BUS] = line.to_bus
        basekv = net["bus"][line.from_bus,:].vn_kv
        baseR = basekv.^2 ./(3*opt["PF"]["baseMVA"])
        branch[:, BR_R] = line.length_km .* line.r_ohm_per_km ./baseR ./line.parallel
        branch[:, BR_X] = line.length_km .* line.x_ohm_per_km ./baseR ./line.parallel
        b = 2*opt["PF"]["net_hz"]*pi*line.length_km .*line.c_nf_per_km .*1e-9 .*baseR .*line.parallel
        g = line.g_us_per_km .* 1e-6 .* baseR .* line.length_km .* line.parallel
        branch[:, RATE_A] .=100
        branch[:, BR_B] = b
        branch[:,BR_STATUS] = line.in_service .== true
        branch[:,ANGMIN] .=-360
        branch[:,ANGMAX] .= 360
    elseif sequence == 0
        line = deepcopy(net["line"])
        branch[:, F_BUS] = line.from_bus
        branch[:, T_BUS] = line.to_bus
        basekv = net["bus"][line.from_bus,:].vn_kv
        baseR = basekv.^2 ./(3*opt["PF"]["baseMVA"])
        branch[:, BR_R] = line.length_km .* line.r0_ohm_per_km ./baseR ./line.parallel
        branch[:, BR_X] = line.length_km .* line.x0_ohm_per_km ./baseR ./line.parallel
        b = 2*opt["PF"]["net_hz"]*pi*line.length_km .*line.c0_nf_per_km .*1e-9 .*baseR .*line.parallel
        g = line.g_us_per_km .* 1e-6 .* baseR .* line.length_km .* line.parallel
        branch[:, BR_B] = b
        branch[:,BR_STATUS] = line.in_service .== true
        branch[:,ANGMIN] .=-360
        branch[:,ANGMAX] .= 360
    end

    jpc["branch"] = branch

    return jpc
end

function calculate_transformer_parameter(net, jpc, opt)
    # Call the indexing function to get the indices for the branch matrix
    (F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, TAP, SHIFT, BR_STATUS, ANGMIN,
    ANGMAX, DICTKEY, PF, QF, PT, QT, MU_SF, MU_ST, MU_ANGMIN, MU_ANGMAX, LAMBDA, SW_TIME, RP_TIME, BR_TYPE, BR_AREA) = PowerFlow.idx_brch()

    branch =jpc["branch"]
    trafo = deepcopy(net["trafo"])
    parallel = trafo.parallel
    branch_trafo = zeros(size(trafo, 1), 14)
    # Initialize the branch matrix for transformers
    branch_trafo[:, F_BUS] = trafo.hv_bus
    branch_trafo[:, T_BUS] = trafo.lv_bus
    r, x, g, b, g_asym, b_asym, ratio, shift = calc_branch_values_from_trafo_df(net, jpc, opt)
    branch_trafo[:, BR_R] = r
    branch_trafo[:, BR_X] = x
    # branch_trafo[:, BR_G] = g
    branch_trafo[:, BR_B] = b
    # branch_trafo[:, BR_G_ASYM] = g_asym
    # branch_trafo[:, BR_B_ASYM] = b_asym
    branch_trafo[:, TAP] = ratio
    branch_trafo[:, SHIFT] = shift
    branch_trafo[:, BR_STATUS] = trafo.in_service .==true

    if any(trafo.df .<= 0)
        error("Rating factor df must be positive. Transformers with false " *
            "rating factors: $(trafo[trafo.df .<= 0, :].index)")
    end

    # always set RATE_A for completeness
    # RATE_A is considered by the (PowerModels) OPF. If zero -> unlimited
    if any(trafo.max_loading_percent) !== missing
        max_load = trafo.max_loading_percent
        sn_mva = trafo.sn_mva
        df = trafo.df
        branch_trafo[:, RATE_A] = max_load ./ 100.0 .* sn_mva .* df .* parallel
    else
        # PowerModels considers "0" as "no limit"
        # todo: inf and convert only when using PowerModels to 0., pypower opf converts the zero to inf
        branch_trafo[:, RATE_A] .= 100.0
    end
        branch_trafo[:, ANGMIN] .= -360.0
        branch_trafo[:, ANGMAX] .= 360.0
        branch = [branch;branch_trafo]
        jpc["branch"] = branch

    return jpc
end

function calc_branch_values_from_trafo_df(net, jpc, opt, trafo_df=nothing, sequence=1)
    # Call the indexing function to get the indices for the bus matrix
    (PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM,VA, 
    BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN, PER_CONSUMER) = PowerFlow.idx_bus();  
    if isnothing(trafo_df)
        trafo_df = deepcopy(net["trafo"])
    end

    lv_bus = get_trafo_values(trafo_df, "lv_bus")
    vn_lv = jpc["bus"][lv_bus, BASE_KV] 
    vn_trafo_hv, vn_trafo_lv, shift = calc_tap_from_dataframe(net, trafo_df)
    ratio = calc_nominal_ratio_from_dataframe(jpc, trafo_df, vn_trafo_hv, vn_trafo_lv)
    r, x, g, b, g_asym, b_asym = calc_r_x_y_from_dataframe(net, trafo_df, vn_trafo_lv, vn_lv, jpc, opt, sequence)

    return r, x, g, b, g_asym, b_asym, ratio, shift
end


function get_trafo_values(trafo_df, par)
    # 获取原始值
    if isa(trafo_df, Dict)
        raw_values = trafo_df[par]
    elseif isa(trafo_df, DataFrame)
        # 对于DataFrame类型的处理
        raw_values = trafo_df[!, par]  # 返回整列数据
    else Symbol(par) in names(trafo_df)
        raw_values = trafo_df[:, Symbol(par)]
    end
    
    # 将字符串"NaN"转换为数值NaN，将字符串"None"转换为nothing
    if isa(raw_values, Vector) || isa(raw_values, Array)
        # 对向量或数组中的每个元素进行处理
        return map(val -> process_value(val), raw_values)
    else
        # 处理单个值
        return process_value(raw_values)
    end
end

# 辅助函数，处理单个值
function process_value(val)
    if isa(val, String)
        if val == "NaN"
            return NaN
        elseif val == "None"
            return nothing
        else
            return val
        end
    else
        return val
    end
end


# TODO:检查函数鲁棒性
function calc_tap_from_dataframe(net, trafo_df)
    vnh = get_trafo_values(trafo_df, "vn_hv_kv")
    vnl = get_trafo_values(trafo_df, "vn_lv_kv")
    trafo_shift = get_trafo_values(trafo_df, "shift_degree")

    for t in ("", "2")
        if !("tap$(t)_pos" in names(trafo_df))
            continue
        end
        
        tap_pos = get_trafo_values(trafo_df, "tap$(t)_pos")
        tap_neutral = get_trafo_values(trafo_df, "tap$(t)_neutral")
        tap_diff = tap_pos - tap_neutral
        tap_side = get_trafo_values(trafo_df, "tap$(t)_side")
        tap_step_percent = get_trafo_values(trafo_df, "tap$(t)_step_percent")
        tap_step_degree = get_trafo_values(trafo_df, "tap$(t)_step_degree")
    
        # 定义三角函数
        cos_func(x) = cos.(deg2rad.(x))
        sin_func(x) = sin.(deg2rad.(x))
        arctan_func(x) = rad2deg.(atan.(x))
    
        if "tap$(t)_changer_type" in names(trafo_df)
            # tap_changer_type is only in dataframe starting from pp Version 3.0, older version use different logic
            tap_changer_type = get_trafo_values(trafo_df, "tap$(t)_changer_type")
            
            if "tap$(t)_dependency_table" in names(trafo_df)
                tap_dependency_table = get_trafo_values(trafo_df, "tap_dependency_table")
                # 替换NaN为false
                tap_dependency_table = map(x -> (typeof(x) <: AbstractFloat && isnan(x)) ? false : x, tap_dependency_table)
            else
                tap_table = fill(false, size(tap_changer_type))
                tap_dependency_table = fill(false, size(tap_changer_type))
            end
            
            # 逻辑运算
            tap_table = tap_dependency_table .& .!(tap_changer_type .== "None")
            tap_no_table = .!tap_dependency_table .& .!(tap_changer_type .== "None")
            
            if any(tap_table)
                id_characteristic_table = get_trafo_values(trafo_df, "id_characteristic_table")
                if any(tap_dependency_table .& ismissing.(id_characteristic_table))
                    error("Trafo with tap_dependency_table True and id_characteristic_table NA detected.\n" *
                          "Please set an id_characteristic_table or set tap_dependency_table to False.")
                end
                
                for (side, vn, direction) in [("hv", vnh, 1), ("lv", vnl, -1)]
                    mask = tap_table .& (side .== tap_side)
                    
                    # 创建DataFrame进行过滤
                    filter_df = DataFrame(
                        id_characteristic = id_characteristic_table,
                        step = tap_pos,
                        mask = mask
                    )
                    
                    # 过滤并合并DataFrame
                    filtered_mask = filter_df[filter_df.mask, :]
                    filtered_df = innerjoin(net.trafo_characteristic_table, filtered_mask, on=[:id_characteristic, :step])
                    
                    cleaned_id_characteristic = id_characteristic_table[(.!ismissing.(id_characteristic_table)) .& mask]
                    
                    # 创建映射字典
                    voltage_mapping = Dict(zip(filtered_df.id_characteristic, filtered_df.voltage_ratio))
                    shift_mapping = Dict(zip(filtered_df.id_characteristic, filtered_df.angle_deg))
                    
                    if direction == 1
                        ratio = [get(voltage_mapping, id_val, 1) for id_val in cleaned_id_characteristic]
                        shift = [get(shift_mapping, id_val, 1) for id_val in cleaned_id_characteristic]
                    else
                        ratio = [get(voltage_mapping, id_val, 1) for id_val in cleaned_id_characteristic]
                        shift = [-get(shift_mapping, id_val, 1) for id_val in cleaned_id_characteristic]
                    end
                    
                    vn[mask] = vn[mask] .* ratio
                    trafo_shift[mask] .+= shift
                end
            end
            
            if any(tap_no_table)
                tap_ideal =  (tap_changer_type .== "Ideal") .& tap_no_table
                tap_complex = ((tap_changer_type .== "Ratio") .| (tap_changer_type .== "Symmetrical")) .& tap_no_table
                
                for (side, vn, direction) in [("hv", vnh, 1), ("lv", vnl, -1)]
                    mask_ideal = tap_ideal .& (tap_side .== side)
                    mask_complex = tap_complex .& (tap_side .== side)
                    
                    if any(mask_ideal)
                        degree_is_set = _replace_nan(tap_step_degree[mask_ideal]) .!= 0
                        percent_is_set = _replace_nan(tap_step_percent[mask_ideal]) .!= 0
                        
                        if any(degree_is_set .& percent_is_set)
                            error("Both tap_step_degree and tap_step_percent set for ideal phase shifter")
                        end
                        
                        trafo_shift[mask_ideal] .+= ifelse.(
                            degree_is_set,
                            (direction .* tap_diff[mask_ideal] .* tap_step_degree[mask_ideal]),
                            (direction .* 2 .* rad2deg.(asin.(tap_diff[mask_ideal] .* 
                                                        tap_step_percent[mask_ideal] ./ 100 ./ 2)))
                        )
                    end
                    
                    if any(mask_complex)
                        tap_steps = tap_step_percent[mask_complex] .* tap_diff[mask_complex] ./ 100
                        tap_angles = _replace_nan(tap_step_degree[mask_complex])
                        u1 = vn[mask_complex]
                        du = u1 .* _replace_nan(tap_steps)
                        vn[mask_complex] = sqrt.((u1 .+ du .* cos_func(tap_angles)) .^ 2 .+ (du .* sin_func(tap_angles)) .^ 2)
                        trafo_shift[mask_complex] .+= (arctan_func(direction .* du .* sin_func(tap_angles) ./
                                                      (u1 .+ du .* cos_func(tap_angles))))
                    end
                end
            end
        elseif haskey(trafo_df, "tap$(t)_phase_shifter")
            @warn("tap$(t)_phase_shifter was removed with pandapower 3.0 and replaced by " *
                  "tap$(t)_changer_type. Using old net data will still work, but usage of " *
                  "tap$(t)_phase_shifter is deprecated and will be removed in future " *
                  "releases.")
                  
            tap_phase_shifter = get_trafo_values(trafo_df, "tap$(t)_phase_shifter")
            
            for (side, vn, direction) in [("hv", vnh, 1), ("lv", vnl, -1)]
                tap_ideal = tap_phase_shifter .& (tap_side .== side)
                tap_complex = isfinite.(tap_step_percent) .& isfinite.(tap_pos) .& (tap_side .== side) .& 
                    .!tap_ideal
                    
                if any(tap_complex)
                    tap_steps = tap_step_percent[tap_complex] .* tap_diff[tap_complex] ./ 100
                    tap_angles = _replace_nan(tap_step_degree[tap_complex])
                    u1 = vn[tap_complex]
                    du = u1 .* _replace_nan(tap_steps)
                    vn[tap_complex] = sqrt.((u1 .+ du .* cos_func(tap_angles)) .^ 2 .+ (du .* sin_func(tap_angles)) .^ 2)
                    trafo_shift[tap_complex] .+= (arctan_func(direction .* du .* sin_func(tap_angles) ./
                                                (u1 .+ du .* cos_func(tap_angles))))
                end
                
                if any(tap_ideal)
                    degree_is_set = _replace_nan(tap_step_degree[tap_ideal]) .!= 0
                    percent_is_set = _replace_nan(tap_step_percent[tap_ideal]) .!= 0
                    
                    if any(degree_is_set .& percent_is_set)
                        error("Both tap_step_degree and tap_step_percent set for ideal phase shifter")
                    end
                    
                    trafo_shift[tap_ideal] .+= ifelse.(
                        degree_is_set,
                        (direction .* tap_diff[tap_ideal] .* tap_step_degree[tap_ideal]),
                        (direction .* 2 .* rad2deg.(asin.(tap_diff[tap_ideal] .*
                                                    tap_step_percent[tap_ideal] ./ 100 ./ 2)))
                    )
                end
            end
        end
    end

    return vnh, vnl, trafo_shift
end


function calc_nominal_ratio_from_dataframe(jpc, trafo_df, vn_trafo_hv, vn_trafo_lv)
    # Call indexing function
    (PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM,VA, 
    BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN, PER_CONSUMER) = PowerFlow.idx_bus()

    tap_rat = vn_trafo_hv ./ vn_trafo_lv
    hv_bus = get_trafo_values(trafo_df, "hv_bus")
    lv_bus = get_trafo_values(trafo_df, "lv_bus")
    
    # 注意：在Julia中，我们使用点运算符来进行向量化操作
    # 假设BASE_KV是一个常量，表示母线电压在jpc["bus"]矩阵中的列索引
    nom_rat = jpc["bus"][hv_bus, BASE_KV]./ jpc["bus"][lv_bus, BASE_KV]
    
    return tap_rat ./ nom_rat
end

function calc_r_x_y_from_dataframe(net, trafo_df, vn_trafo_lv, vn_lv, jpc, opt, sequence=1)
    characteristic=get(net, "characteristic", nothing)
    r, x = calc_r_x_from_dataframe(trafo_df, vn_lv, vn_trafo_lv, opt, sequence, characteristic )

    g, b = calc_y_from_dataframe( trafo_df, vn_lv, vn_trafo_lv,opt["PF"]["baseMVA"])

    r_ratio = fill(0.5, size(r)...)
        
    x_ratio = fill(0.5, size(x)...)
        
    return wye_delta(r, x, g, b, r_ratio, x_ratio)

end


function calc_r_x_from_dataframe(trafo_df, vn_lv, vn_trafo_lv, opt, sequence=1, characteristic=nothing, trafo_characteristic_table=nothing)
    """
    Calculates (Vectorized) the resistance and reactance according to the
    transformer values
    """
    parallel = get_trafo_values(trafo_df, "parallel")

    if sequence == 1
        vk_percent, vkr_percent = get_vk_values(trafo_df)
    elseif sequence == 0
        vk_percent = get_trafo_values(trafo_df, "vk0_percent")
        vkr_percent = get_trafo_values(trafo_df, "vkr0_percent")
    else
        error("Unsupported sequence")
    end

    # adjust for low voltage side voltage converter:
    sn_mva = opt["PF"]["baseMVA"]
    tap_lv = (vn_trafo_lv ./ vn_lv).^2 .* (3 * sn_mva)


    sn_trafo_mva = get_trafo_values(trafo_df, "sn_mva")
    z_sc = vk_percent ./ 100.0 ./ sn_trafo_mva .* tap_lv
    r_sc = vkr_percent / 100.0 ./ sn_trafo_mva .* tap_lv
    x_sc = sign.(z_sc) .* sqrt.(Float64.((z_sc.^2 - r_sc.^2)))

    return r_sc ./ parallel, x_sc ./ parallel
end


function get_vk_values(trafo_df, trafotype="2W")
    if trafotype == "2W"
        vk_variables = ("vk_percent", "vkr_percent")
    elseif trafotype == "3W"
        vk_variables = ("vk_hv_percent", "vkr_hv_percent", "vk_mv_percent", "vkr_mv_percent",
                        "vk_lv_percent", "vkr_lv_percent")
    else
        error("Unknown trafotype")
    end

    vals = ()

    for vk_var in vk_variables
        vk_value = get_trafo_values(trafo_df, vk_var)
        vals = (vals..., vk_value)
    end

    return vals
end

function calc_y_from_dataframe(trafo_df, vn_lv, vn_trafo_lv, net_sn_mva)

    baseZ =vn_lv.^2 ./ (3*net_sn_mva)
    vn_lv_kv = get_trafo_values(trafo_df, "vn_lv_kv")
    
    pfe_mw = (get_trafo_values(trafo_df, "pfe_kw") * 1e-3) / 3

    parallel = get_trafo_values(trafo_df, "parallel")
    trafo_sn_mva = get_trafo_values(trafo_df, "sn_mva")

    ### Calculate susceptance ###
    vnl_squared = (vn_lv_kv.^2)./3 
    g_mva = pfe_mw
    
    i0 = get_trafo_values(trafo_df, "i0_percent") / 3
  
    ym_mva = i0 ./ 100 .* trafo_sn_mva
    b_mva_squared = ym_mva.^2 - pfe_mw.^2
    b_mva_squared[b_mva_squared .< 0] .= 0
    b_mva = -sqrt.(b_mva_squared)

    g_pu = g_mva ./ vnl_squared .* baseZ .* parallel ./ ((vn_trafo_lv ./ vn_lv_kv).^2)
    b_pu = b_mva ./ vnl_squared .* baseZ .* parallel ./ ((vn_trafo_lv ./ vn_lv_kv).^2)
    
    return g_pu, b_pu
end

function wye_delta(r, x, g, b, r_ratio, x_ratio)
    tidx = (g .!= 0) .| (b .!= 0)
    # 创建临时数组来存储计算结果
    za_star = zeros(Complex{Float64}, count(tidx))
    zb_star = zeros(Complex{Float64}, count(tidx))
    zc_star = zeros(Complex{Float64}, count(tidx))
    
    # 对满足条件的元素进行计算
    za_star = r[tidx] .* r_ratio[tidx] .+ x[tidx] .* x_ratio[tidx] .* im
    zb_star = r[tidx] .* (1 .- r_ratio[tidx]) .+ x[tidx] .* (1 .- x_ratio[tidx]) .* im
    zc_star = 1 ./ (g[tidx] .+ im .* b[tidx])
    
    zSum_triangle = za_star .* zb_star .+ za_star .* zc_star .+ zb_star .* zc_star
    zab_triangle = zSum_triangle ./ zc_star
    zac_triangle = zSum_triangle ./ zb_star
    zbc_triangle = zSum_triangle ./ za_star
    
    r[tidx] = real.(zab_triangle)
    x[tidx] = imag.(zab_triangle)
    
    yf = 1 ./ zac_triangle
    yt = 1 ./ zbc_triangle
    
    # 2 because in makeYbus Bcf, Bct are divided by 2:
    g[tidx] = real.(yf) .* 2
    b[tidx] = imag.(yf) .* 2
    
    g_asym = zeros(size(g))
    b_asym = zeros(size(b))
    
    g_asym[tidx] = 2 .* real.(yt) .- g[tidx]
    b_asym[tidx] = 2 .* imag.(yt) .- b[tidx]
    
    return r, x, g, b, g_asym, b_asym
end

function build_branch_Jpc_zero(net, jpc_new, opt, k_st = nothing)
    # Call the indexing function to get the indices for the branch matrix
    (F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, TAP, SHIFT, BR_STATUS, ANGMIN,
    ANGMAX, DICTKEY, PF, QF, PT, QT, MU_SF, MU_ST, MU_ANGMIN, MU_ANGMAX, LAMBDA, SW_TIME, RP_TIME, BR_TYPE, BR_AREA) = PowerFlow.idx_brch()
    length = initialize_branch_length(net)
    jpc_new["branch"] = zeros( length, 14)
    jpc_new["branch"][:, RATE_A] .= 250
    jpc_new["branch"][:, RATE_B] .= 250
    jpc_new["branch"][:, RATE_C] .= 250
    jpc_new["branch"][:, TAP] .= 1
    jpc_new["branch"][:, SHIFT] .= 0
    jpc_new["branch"][:, BR_STATUS] .= 1
    jpc_new["branch"][:, ANGMIN] .= -360
    jpc_new["branch"][:, ANGMAX] .= 360

    add_line_sc_impedance_zero(net, jpc_new, opt)
    add_trafo_sc_impedance_zero(net, jpc_new, k_st)

    return jpc_new

end

function initialize_branch_length(net)
    Start = 0
    End = 0
    elements = ["line", "trafo", "trafo3w", "impedance", "xward"]
    for element in elements
        if element in keys(net)
            if element == "trafo3w"
                End = Start + size(net[element],1) * 3
            else
                End = Start + size(net[element],1)
            end
            Start = End
        end
    end
    return End
end

function add_line_sc_impedance_zero(net, jpc, opt)
    # Call the indexing function to get the indices for the bus matrix
    (PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM,VA, 
    BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN, PER_CONSUMER) = PowerFlow.idx_bus();
    (F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, TAP, SHIFT, BR_STATUS, ANGMIN,
    ANGMAX, DICTKEY, PF, QF, PT, QT, MU_SF, MU_ST, MU_ANGMIN, MU_ANGMAX, LAMBDA, SW_TIME, RP_TIME, BR_TYPE, BR_AREA) = PowerFlow.idx_brch()

    if !haskey(net, "line")
        return
    end
    
    line = net["line"]
    length = line.length_km
    parallel = line.parallel

    fb = line.from_bus
    tb = line.to_bus
    # baseR = jpc["bus"][fb, BASE_KV].^2 ./ opt["PF"]["baseMVA"]
   
    baseR = jpc["bus"][fb, BASE_KV].^2 ./ (3 * opt["PF"]["baseMVA"])

    
    f = 1
    t = size(line, 1)
    
    # line zero sequence impedance
    jpc["branch"][f:t, F_BUS] = fb
    jpc["branch"][f:t, T_BUS] = tb
    jpc["branch"][f:t, BR_R] = line.r0_ohm_per_km .* length ./ baseR ./ parallel
    
    jpc["branch"][f:t, BR_X] = line.x0_ohm_per_km .* length ./ baseR ./ parallel
    jpc["branch"][f:t, BR_B] = (2 * opt["PF"]["net_hz"] * π * line.c0_nf_per_km .*
                               1e-9 .* baseR .* length .* parallel)
    jpc["branch"][f:t, BR_STATUS] = line.in_service .== true
end

function add_trafo_sc_impedance_zero(net, jpc, trafo_df=nothing, k_st=nothing)
    (PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM,VA, 
    BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN, PER_CONSUMER) = PowerFlow.idx_bus();
    (F_BUS, T_BUS, BR_R, BR_X, BR_B, RATE_A, RATE_B, RATE_C, TAP, SHIFT, BR_STATUS, ANGMIN,
    ANGMAX, DICTKEY, PF, QF, PT, QT, MU_SF, MU_ST, MU_ANGMIN, MU_ANGMAX, LAMBDA, SW_TIME, RP_TIME, BR_TYPE, BR_AREA) = PowerFlow.idx_brch()
    
    if trafo_df === nothing
        trafo_df = net["trafo"]
    end
    
    if k_st === nothing
        k_st = ones(size(jpc["branch"],1))
    end
    
    if !("xn_ohm" in names(trafo_df))
        trafo_df[!, "xn_ohm"] .= 0.0
    end
    
    if !haskey(net, "trafo")
        return
    end
    
    f = size(net["line"], 1) + 1
    t = f + size(trafo_df, 1) - 1
    trafo_df[!, "_jpc_idx"] = collect(range(f, t))
    trafo_df.k_st = real.(k_st[trafo_df._jpc_idx])

    hv_bus = convert.(Int64, get_trafo_values(trafo_df, "hv_bus"))
    lv_bus = convert.(Int64, get_trafo_values(trafo_df, "lv_bus"))
    in_service = get_trafo_values(trafo_df, "in_service") .==1
    jpc["branch"][f:t, F_BUS] = hv_bus
    jpc["branch"][f:t, T_BUS] = lv_bus
    buses_all, gs_all, bs_all = Int64[], Float64[], Float64[]
    
    jpc["branch"][f:t, BR_STATUS] .= 0


    if !("vector_group" in names(trafo_df))
        error("Vector Group of transformer needs to be specified for zero " *
              "sequence modelling \n Try : net.trafo[\"vector_group\"] = 'Dyn'")
    end

    # 按vector_group分组处理
    for (vector_group, trafos) in pairs(groupby(trafo_df, :vector_group))
        (; vector_group) = vector_group
        # TODO Roman: check this/expand this
        jpc_idx = convert.(Int64, trafos._jpc_idx)

        if lowercase(vector_group) in ["yy", "yd", "dy", "dd"]
            continue
        end

        # vk_percent = convert.(Float64, trafos.vk_percent)
        # vkr_percent = convert.(Float64, trafos.vkr_percent)
        sn_trafo_mva = convert.(Float64, trafos.sn_mva)
        
        # Just put pos seq parameter if zero seq parameter is zero
        if !("vk0_percent" in names(trafos))
            error("Short circuit voltage of transformer Vk0 needs to be specified for zero " *
                  "sequence modelling \n Try : net.trafo[\"vk0_percent\"] = net.trafo[\"vk_percent\"]")
        end
        
        vk0_percent = all(convert.(Float64, trafos.vk0_percent) .!= 0.0) ? 
                      convert.(Float64, trafos.vk0_percent) : 
                      convert.(Float64, trafos.vk_percent)
        
        # Just put pos seq parameter if zero seq parameter is zero
        if !("vkr0_percent" in names(trafos))
            error("Real part of short circuit voltage Vk0(Real) needs to be specified for transformer " *
                  "modelling \n Try : net.trafo[\"vkr0_percent\"] = net.trafo[\"vkr_percent\"]")
        end
        
        vkr0_percent = all(convert.(Float64, trafos.vkr0_percent) .!= 0.0) ? 
                       convert.(Float64, trafos.vkr0_percent) : 
                       convert.(Float64, trafos.vkr_percent)
        
        lv_buses = convert.(Int64, trafos.lv_bus)
        hv_buses = convert.(Int64, trafos.hv_bus)
        lv_buses_jpc = lv_buses
        hv_buses_jpc = hv_buses
        
        if !("mag0_percent" in names(trafos))
            error("Magnetizing impedance to vk0 ratio needs to be specified for transformer " *
                  "modelling  \n Try : net.trafo[\"mag0_percent\"] = 100")
        end
        
        mag0_ratio = convert.(Float64, trafos.mag0_percent)
        
        if !("mag0_rx" in names(trafos))
            error("Magnetizing impedance R/X ratio needs to be specified for transformer " *
                  "modelling \n Try : net.trafo[\"mag0_rx\"] = 0 ")
        end
        
        mag0_rx = convert.(Float64, trafos.mag0_rx)
        
        if !("si0_hv_partial" in names(trafos))
            error("Zero sequence short circuit impedance partition towards HV side needs to be specified " *
                  "for transformer modelling \n Try : net.trafo[\"si0_hv_partial\"] = 0.9 ")
        end
        
        si0_hv_partial = convert.(Float64, trafos.si0_hv_partial)
        parallel = convert.(Float64, trafos.parallel)
        
        if "power_station_unit" in names(trafos)
            power_station_unit = convert.(Bool, coalesce.(trafos.power_station_unit, false))
        else
            power_station_unit = zeros(Bool, length(trafos))
        end
        
        in_service = trafos.in_service .== true

        jpc["branch"][jpc_idx, F_BUS] = hv_buses_jpc
        jpc["branch"][jpc_idx, T_BUS] = lv_buses_jpc

        vn_trafo_hv, vn_trafo_lv, shift = calc_tap_from_dataframe(net, trafos)
        vn_bus_lv = jpc["bus"][lv_buses_jpc, BASE_KV]
        vn_bus_hv = jpc["bus"][hv_buses_jpc, BASE_KV]
        ratio = calc_nominal_ratio_from_dataframe(jpc, trafos, vn_trafo_hv,
                                                  vn_trafo_lv)
        jpc["branch"][jpc_idx, TAP] = ratio
        jpc["branch"][jpc_idx, SHIFT] = shift

        # zero seq. transformer impedance
        tap_lv = (vn_trafo_lv ./ vn_bus_lv).^2 .* jpc["baseMVA"]
        tap_hv = (vn_trafo_hv ./ vn_bus_hv).^2 .* jpc["baseMVA"]
        
        if lowercase(vector_group) ∉ ["ynyn", "dyn", "yzn"]
            error("Calculation of 3-phase power flow is only implemented for the transformer " *
                    "vector groups 'YNyn', 'Dyn', 'Yzn'")
        end
        # =============================================================================
        #     Changing base from transformer base to Network base to get Zpu(Net)
        #     Zbase = (kV).squared/S_mva
        #     Zpu(Net)={Zpu(trafo) * Zb(trafo)} / {Zb(Net)}
        #        Note:
        #             Network base voltage is Line-Neutral voltage in each phase
        #             Line-Neutral voltage= Line-Line Voltage(vn_bus_lv) divided by sq.root(3)
        # =============================================================================
        tap_lv = (vn_trafo_lv ./ vn_bus_lv).^2 .* (3 * jpc["baseMVA"])
        tap_hv = (vn_trafo_hv ./ vn_bus_hv).^2 .* (3 * jpc["baseMVA"])

        tap_corr = lowercase(vector_group) in ("ynd", "yny") ? tap_hv : tap_lv
        z_sc = vk0_percent ./ 100.0 ./ sn_trafo_mva .* tap_corr
        r_sc = vkr0_percent ./ 100.0 ./ sn_trafo_mva .* tap_corr
        z_sc = convert.(Float64, z_sc)
        r_sc = convert.(Float64, r_sc)
        x_sc = sign.(z_sc) .* sqrt.(z_sc.^2 .- r_sc.^2)
        # TODO: This equation needs to be checked!
        # z0_k = (r_sc + x_sc * 1j) / parallel  * max(1, ratio) **2
        # z0_k = (r_sc + x_sc * 1j) / parallel * vn_trafo_hv / vn_bus_hv
        # z0_k = (r_sc + x_sc * 1j) / parallel * tap_hv
        z0_k = (r_sc .+ im .* x_sc) ./ parallel

        y0_k = 1 ./ z0_k  # adding admittance for "pi" model
        # y0_k = 1 / (z0_k * k_st_tr + 3j * z_n_ohm)  # adding admittance for "pi" model

        # =============================================================================
        #       Transformer magnetising impedance for zero sequence
        # =============================================================================
        z_m = z_sc .* mag0_ratio
        x_m = z_m ./ sqrt.(mag0_rx.^2 .+ 1)
        r_m = x_m .* mag0_rx
        r0_trafo_mag = r_m ./ parallel
        x0_trafo_mag = x_m ./ parallel
        z0_mag = r0_trafo_mag .+ im .* x0_trafo_mag
        # =============================================================================
        #         Star - Delta conversion ( T model to Pi Model)
        #      ----------- |__zc=ZAB__|-----------------
        #            _|                   _|
        #     za=ZAN|_|                  |_| zb=ZBN
        #            |                    |
        # =============================================================================
        z1 = si0_hv_partial .* z0_k
        z2 = (1 .- si0_hv_partial) .* z0_k
        z3 = z0_mag
        z_temp = z1 .* z2 .+ z2 .* z3 .+ z1 .* z3
        za = z_temp ./ z2
        # za = z_temp / (z2+z3)
        zb = z_temp ./ z1
        # zb = z_temp / (z1+z3)
        zc = z_temp ./ z3  # ZAB  Transfer impedance
        # zc = z_temp / (z1+z2)  # ZAB  Transfer impedance
        YAB = 1 ./ convert.(ComplexF64, zc)
        YAN = 1 ./ convert.(ComplexF64, za)
        YBN = 1 ./ convert.(ComplexF64, zb)

        # YAB_AN = (zc + za) /(zc * za).astype(complex)  # Series conn YAB and YAN
        # YAB_BN = (zc + zb) / (zc * zb).astype(complex)  # Series conn YAB and YBN

        YAB_AN = 1 ./ convert.(ComplexF64, zc .+ za)  # Series conn YAB and YAN
        YAB_BN = 1 ./ convert.(ComplexF64, zc .+ zb)  # Series conn YAB and YBN

        # y0_k = 1 / z0_k #adding admittance for "pi" model
        if lowercase(vector_group) == "dyn"
            buses_all = vcat(buses_all, lv_buses_jpc)

            y = convert.(ComplexF64, (YAB .+ YBN)) .* jpc["baseMVA"]  # T model

            gs_all = vcat(gs_all, real.(y) .* in_service)
            bs_all = vcat(bs_all, imag.(y) .* in_service)

        elseif lowercase(vector_group) == "ynd"
            buses_all = vcat(buses_all, hv_buses_jpc)

            y = convert.(ComplexF64, (YAB_BN .+ YAN)) .* jpc["baseMVA"]  # T model

            gs_all = vcat(gs_all, real.(y) .* in_service)
            bs_all = vcat(bs_all, imag.(y) .* in_service)

        elseif lowercase(vector_group) == "yyn"
            buses_all = vcat(buses_all, lv_buses_jpc)
            y = convert.(ComplexF64, (YAB .+ YAB_BN .+ YBN)) .* jpc["baseMVA"]  # T model
            gs_all = vcat(gs_all, real.(y) .* in_service)
            bs_all = vcat(bs_all, imag.(y) .* in_service)

        elseif lowercase(vector_group) == "ynyn"
            jpc["branch"][jpc_idx, BR_STATUS] = in_service
            # zc = ZAB
            jpc["branch"][jpc_idx, BR_R] = real.(zc)
            jpc["branch"][jpc_idx, BR_X] = imag.(zc)

            buses_all = vcat(buses_all, hv_buses_jpc)
            gs_all = vcat(gs_all, real.(YAN) .* in_service .* jpc["baseMVA"] .* tap_lv ./ tap_hv)
            bs_all = vcat(bs_all, imag.(YAN) .* in_service .* jpc["baseMVA"] .* tap_lv ./ tap_hv)

            buses_all = vcat(buses_all, lv_buses_jpc)
            gs_all = vcat(gs_all, real.(YBN) .* in_service .* jpc["baseMVA"])
            bs_all = vcat(bs_all, imag.(YBN) .* in_service .* jpc["baseMVA"])

        elseif lowercase(vector_group) == "yny"
            buses_all = vcat(buses_all, hv_buses_jpc)
            y = convert.(ComplexF64, (YAB_BN .+ YAN)) .* jpc["baseMVA"]  # T model
            gs_all = vcat(gs_all, real.(y) .* in_service)
            bs_all = vcat(bs_all, imag.(y) .* in_service)

        elseif lowercase(vector_group) == "yzn"
            buses_all = vcat(buses_all, lv_buses_jpc)
            # y = 1/(z0_mag+z0_k).astype(complex)* int(jpc["baseMVA"])#T model
            # y= (za+zb+zc)/((za+zc)*zb).astype(complex)* int(jpc["baseMVA"])#pi model
            y = convert.(ComplexF64, (YAB_AN .+ YBN)) .* jpc["baseMVA"]^2  # T model # why sn_mva squared here?
            gs_all = vcat(gs_all, (1.1547) .* real.(y) .* in_service)  # what's the 1.1547 value?
            bs_all = vcat(bs_all, (1.1547) .* imag.(y) .* in_service)

        elseif isdigit(vector_group[end])
            error("Unknown transformer vector group $vector_group - " *
                  "please specify vector group without " *
                  "phase shift number. Phase shift can be " *
                  "specified in net.trafo.shift_degree")
        else
            error("Transformer vector group $vector_group is unknown " *
                  "/ not implemented for three phase load flow")
        end
    end

    buses, gs, bs = PowerFlow.sum_by_group(buses_all, gs_all, bs_all)
    jpc["bus"][buses, GS] .+= gs
    jpc["bus"][buses, BS] .+= bs
    select!(net["trafo"], Not(:_jpc_idx))
end
