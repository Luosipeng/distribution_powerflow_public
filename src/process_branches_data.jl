function process_branches_data(cable_data::DataFrame, transline_data::DataFrame, impedance_data::DataFrame, transformer_data::DataFrame, transformer_three_data::DataFrame, HVCB_data::DataFrame, bus::Matrix{Float64}, baseMVA::Float64, baseKV::Float64, dict_bus)
    # --- Call indexing for IEEE branch data ---
    (FBUS, TBUS, R, X, B, RATEA, RATEB, RATEC, RATIO, ANGLE, STATUS, ANGMIN,
     ANGMAX, DICTKEY, PF, QF, PT, QT, MU_SF, MU_ST, MU_ANGMIN, MU_ANGMAX, LAMBDA, SW_TIME, RP_TIME, BR_TYPE, BR_AREA) = PowerFlow.idx_brch()
     (PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM,VA, 
     BASEKV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN, PER_CONSUMER) = PowerFlow.idx_bus();
    (Cable_EquipeID,Cable_inservice,Cable_Felement,Cable_Telement,Cable_length,Cable_r1,Cable_x1,Cable_y1,Cable_perunit_length_value)=PowerFlow.cable_idx()#线缆索引
    (Line_EquipmentID,Line_length,Line_inservice,Line_Felement,Line_Telement,Line_r11,Line_x1,Line_y1)=PowerFlow.xline_idx()#传输线索引
    (Trans_EquipmentID,Trans_inservice,Trans_Pelement,Trans_Selement,prisecrated_Power,Trans_Pvoltage,Trans_Svoltage,Trans_Pos_value,Trans_Pos_XRrating)=PowerFlow.XFORM2W_idx()#变压器索引
    (HVCB_ID,HVCB_FROM_ELEMENT,HVCB_TO_ELEMENT,HVCB_INSERVICE,HVCB_STATUS)=PowerFlow.hvcb_idx()#HVCB索引
    (imp_ID,imp_inservice,imp_Felement,imp_Telement,imp_R,imp_X)=PowerFlow.imp_idx()#阻抗索引
    
    # 初始化空的分支数组
    branch_cable = zeros(0, 14)
    branch_transline = zeros(0, 14)
    branch_transformer = zeros(0, 14)
    branch_impedance = zeros(0, 14)
    if !isempty(transformer_three_data)
        branch_trans_three, bus, dict_bus = process_three_winding_transformer(transformer_three_data, bus, baseMVA, baseKV, dict_bus)
    else
        branch_trans_three=[]
    end
    # --- 处理线缆数据 ---
    if size(cable_data, 1) > 0
        # --- Clean data by removing extra branches ---
        cable = filter(row -> row[Cable_Felement] !== missing, cable_data)
        cable = filter(row -> row[Cable_Telement] !== missing, cable)
        
        if size(cable, 1) > 0
            # --- Initialize matrices ---
            branch_cable = zeros(size(cable, 1), 14)
            
            # --- Transfer data ---
            branch_cable[:,STATUS] = (cable[:, Cable_inservice] .== "true").+0
            
            # --- Assign IDs ---
            branch_cable[:, FBUS] = map(k -> dict_bus[k], cable[:, Cable_Felement])
            branch_cable[:, TBUS] = map(k -> dict_bus[k], cable[:, Cable_Telement])
            
            #Transfer string to Float64
            cable[:, Cable_perunit_length_value] = parse.(Float64, cable[:, Cable_perunit_length_value])
            cable[:, Cable_length] = parse.(Float64, cable[:, Cable_length])
            cable[:, Cable_r1] = parse.(Float64, cable[:, Cable_r1])
            cable[:, Cable_x1] = parse.(Float64, cable[:, Cable_x1])
            
            # --- Assign R, X, B for each branch type ---
            cable[cable[:, Cable_perunit_length_value] .== 0, Cable_perunit_length_value] .= 1000
            branch_cable[:, R] = cable[:, Cable_length] .* cable[:, Cable_r1] ./ cable[:, Cable_perunit_length_value]
            branch_cable[:, X] = cable[:, Cable_length] .* cable[:, Cable_x1] ./ cable[:, Cable_perunit_length_value]
            branch_cable[:, B] .= 0
            
            # --- Assign Rate and Status ---
            rate_values = 100  # MVA
            branch_cable[:, RATEA] .= rate_values
            branch_cable[:, RATEB] .= rate_values
            branch_cable[:, RATEC] .= rate_values
            branch_cable[:, RATIO] .= 0
            branch_cable[:, ANGLE] .= 0
            branch_cable[:, ANGMIN] .= -180
            branch_cable[:, ANGMAX] .= 180
            
            # --- Normalize impedance to per unit ---
            cable_basekv = bus[Int.(branch_cable[:, FBUS]), BASEKV]
            branch_cable[:, R] .= branch_cable[:, R] .* baseMVA .* inv.(cable_basekv.^2)
            branch_cable[:, X] .= branch_cable[:, X] .* baseMVA .* inv.(cable_basekv.^2)
            branch_cable[:, B] .= branch_cable[:, B] ./( baseMVA .* inv.(cable_basekv.^2))
        end
    end

    # --- 处理传输线数据 ---
    if size(transline_data, 1) > 0
        # --- Clean data by removing extra branches ---
        transline = filter(row -> row[Line_Felement] !== missing, transline_data)
        transline = filter(row -> row[Line_Telement] !== missing, transline)
        
        if size(transline, 1) > 0
            # --- Initialize matrices ---
            branch_transline = zeros(size(transline, 1), 14)
            
            # --- Transfer data ---
            branch_transline[:,STATUS] = (transline[:, Line_inservice] .== "true").+0
            
            # --- Assign IDs ---
            branch_transline[:, FBUS] = map(k -> dict_bus[k], transline[:, Line_Felement])
            branch_transline[:, TBUS] = map(k -> dict_bus[k], transline[:, Line_Telement])
            
            #Transfer string to Float64
            transline[:, Line_length] = parse.(Float64, transline[:, Line_length])
            transline[:, Line_r11] = parse.(Float64, transline[:, Line_r11])
            transline[:, Line_x1] = parse.(Float64, transline[:, Line_x1])
            
            # --- Assign R, X, B for each branch type ---
            branch_transline[:, R] = transline[:, Line_length] .* transline[:, Line_r11] .* 0.000621371
            branch_transline[:, X] = transline[:, Line_length] .* transline[:, Line_x1] .* 0.000621371
            branch_transline[:, B] .= 0
            
            # --- Assign Rate and Status ---
            rate_values = 100  # MVA
            branch_transline[:, RATEA] .= rate_values
            branch_transline[:, RATEB] .= rate_values
            branch_transline[:, RATEC] .= rate_values
            branch_transline[:, RATIO] .= 0
            branch_transline[:, ANGLE] .= 0
            branch_transline[:, ANGMIN] .= -180
            branch_transline[:, ANGMAX] .= 180
            
            # --- Normalize impedance to per unit ---
            trans_basekv = bus[Int.(branch_transline[:, FBUS]), BASEKV]
            branch_transline[:, R] .= branch_transline[:, R] .* baseMVA .* inv.(trans_basekv.^2)
            branch_transline[:, X] .= branch_transline[:, X] .* baseMVA .* inv.(trans_basekv.^2)
            branch_transline[:, B] .= branch_transline[:, B] ./( baseMVA .* inv.(trans_basekv.^2))
        end
    end

    # --- 处理变压器数据 ---
    if size(transformer_data, 1) > 0
        # --- Clean data by removing extra branches ---
        transformer = filter(row -> row[Trans_Pelement] !== missing, transformer_data)
        transformer = filter(row -> row[Trans_Selement] !== missing, transformer)
        
        if size(transformer, 1) > 0
            # --- Initialize matrices ---
            branch_transformer = zeros(size(transformer, 1), 14)
            
            # --- Transfer data ---
            branch_transformer[:,STATUS] = (transformer[:, Trans_inservice] .== "true").+0
            
            # --- Assign IDs ---
            branch_transformer[:, FBUS] = map(k -> dict_bus[k], transformer[:, Trans_Pelement])
            branch_transformer[:, TBUS] = map(k -> dict_bus[k], transformer[:, Trans_Selement])
            
            #Transfer string to Float64
            transformer[:,Trans_Pos_value]=parse.(Float64,transformer[:, Trans_Pos_value])
            transformer[:, Trans_Pvoltage] = parse.(Float64,transformer[:, Trans_Pvoltage])
            transformer[:, prisecrated_Power] = parse.(Float64,transformer[:, prisecrated_Power])
            transformer[:, Trans_Pos_XRrating] = parse.(Float64,transformer[:, Trans_Pos_XRrating])
            
            # --- Assign R, X, B for each branch type ---
            Zbase = baseKV^2 / baseMVA
            branch_transformer[:, R] = 0.01 .* transformer[:, Trans_Pos_value] .* transformer[:, Trans_Pvoltage].^2 .* inv.(transformer[:, prisecrated_Power]./1000) .* inv.(sqrt.(1 .+ transformer[:, Trans_Pos_XRrating].^2)) 
            branch_transformer[:, X] = branch_transformer[:, R] .* transformer[:, Trans_Pos_XRrating]
            branch_transformer[:, B] .= 0
            
            # --- Assign Rate and Status ---
            rate_values = 100  # MVA
            branch_transformer[:, RATEA] .= rate_values
            branch_transformer[:, RATEB] .= rate_values
            branch_transformer[:, RATEC] .= rate_values
            branch_transformer[:, RATIO] .= 1.0
            branch_transformer[:, ANGLE] .= 0
            branch_transformer[:, ANGMIN] .= -180
            branch_transformer[:, ANGMAX] .= 180
            
            # --- Normalize impedance to per unit ---
            transformer_basekv = transformer[:, Trans_Pvoltage]
            branch_transformer[:, R] .= branch_transformer[:, R] .* baseMVA .* inv.(transformer_basekv.^2)
            branch_transformer[:, X] .= branch_transformer[:, X] .* baseMVA .* inv.(transformer_basekv.^2)
            branch_transformer[:, B] .= branch_transformer[:, B] ./( baseMVA .* inv.(transformer_basekv.^2))
        end
    end

    if size(impedance_data,1)>0
        impedance_data=filter(row->row[imp_Felement]!==missing,impedance_data)
        impedance_data=filter(row->row[imp_Telement]!==missing,impedance_data)
        if size(impedance_data,1)>0
            branch_impedance=zeros(size(impedance_data,1),14)
            branch_impedance[:,STATUS]=(impedance_data[:,imp_inservice].=="true").+0
            branch_impedance[:,FBUS]=map(k->dict_bus[k],impedance_data[:,imp_Felement])
            branch_impedance[:,TBUS]=map(k->dict_bus[k],impedance_data[:,imp_Telement])
            impedance_data[:,imp_R]=parse.(Float64,impedance_data[:,imp_R])./100
            impedance_data[:,imp_X]=parse.(Float64,impedance_data[:,imp_X])./100
            branch_impedance[:,R]=impedance_data[:,imp_R]
            branch_impedance[:,X]=impedance_data[:,imp_X]
            branch_impedance[:,B].=0
            rate_values=100
            branch_impedance[:,RATEA].=rate_values
            branch_impedance[:,RATEB].=rate_values
            branch_impedance[:,RATEC].=rate_values
            branch_impedance[:,RATIO].=0
            branch_impedance[:,ANGLE].=0
            branch_impedance[:,ANGMIN].=-180
            branch_impedance[:,ANGMAX].=180
            impedance_basekv=bus[Int.(branch_impedance[:,FBUS]),BASEKV]
            branch_impedance[:,R].=branch_impedance[:,R].*baseMVA.*inv.(impedance_basekv.^2)
            branch_impedance[:,X].=branch_impedance[:,X].*baseMVA.*inv.(impedance_basekv.^2)
            branch_impedance[:,B].=branch_impedance[:,B]./(baseMVA.*inv.(impedance_basekv.^2))
        end
    end
    # --- Combine all branches ---
    branch = [branch_cable; branch_transline;branch_impedance; branch_transformer]

    # --- Topology reconstruction ---
    # 处理HVCB数据
    if size(HVCB_data, 1) > 0
        # Update HVCB data
        HVCB_data = filter(row -> row[HVCB_FROM_ELEMENT] !== missing, HVCB_data)
        HVCB_data = filter(row -> row[HVCB_TO_ELEMENT] !== missing, HVCB_data)
        
        if size(HVCB_data, 1) > 0
            HVCB_data = HVCB_data[findall(HVCB_data[:, HVCB_INSERVICE] .== "true"), :]
            HVCB_data = HVCB_data[.!ismissing.(HVCB_data[:, HVCB_FROM_ELEMENT]), :]
            HVCB_data = HVCB_data[.!ismissing.(HVCB_data[:, HVCB_TO_ELEMENT]), :]
            HVCB_data[:, HVCB_FROM_ELEMENT] = map(k -> dict_bus[k], HVCB_data[:, HVCB_FROM_ELEMENT])
            HVCB_data[:, HVCB_TO_ELEMENT] = map(k -> dict_bus[k], HVCB_data[:, HVCB_TO_ELEMENT])
        end
    end

    return branch, HVCB_data, bus, dict_bus
end

function process_three_winding_transformer(transformer_three_data, bus, baseMVA, baseKV, dict_bus)
    (XFORM3W_ID,XFORM3W_INSERVICE,XFORM3W_STATE,XFORM3W_PRI_ELEMENT,XFORM3W_SEC_ELEMENT,
    XFORM3W_TERT_ELEMENT,XFORM3W_PRI_KV,XFORM3W_SEC_KV,XFORM3W_TERT_KV,XFORM3W_PRI_TAP,
    XFORM3W_SEC_TAP,XFORM3W_TERT_TAP,XFORM3W_PRI_KVA,XFORM3W_SEC_KVA,XFORM3W_TERT_KVA,
    XFORM3W_PRI_POSZ,XFORM3W_PRI_X_R,XFORM3W_SEC_POSZ,XFORM3W_SEC_X_R,
    XFORM3W_TERT_POSZ,XFORM3W_TERT_X_R)=PowerFlow.XFORM3W_idx()

    (FBUS, TBUS, R, X, B, RATEA, RATEB, RATEC, RATIO, ANGLE, STATUS, ANGMIN,
     ANGMAX, DICTKEY, PF, QF, PT, QT, MU_SF, MU_ST, MU_ANGMIN, MU_ANGMAX, LAMBDA, SW_TIME, RP_TIME, BR_TYPE, BR_AREA) = PowerFlow.idx_brch()
     (PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM,VA, 
     BASEKV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN, PER_CONSUMER) = PowerFlow.idx_bus();
    
    # --- 处理三相变压器数据 ---
    transformer_three_data = filter(row -> row[XFORM3W_INSERVICE] .== "true", transformer_three_data)
    transformer_three_data = filter(row -> row[XFORM3W_PRI_ELEMENT] !== missing, transformer_three_data)
    transformer_three_data = filter(row -> row[XFORM3W_SEC_ELEMENT] !== missing, transformer_three_data)
    transformer_three_data = filter(row -> row[XFORM3W_TERT_ELEMENT] !== missing, transformer_three_data)

    # Create center bus
    center_bus = [string("BUS_Three_wind_trans_b", i) for i in eachindex(transformer_three_data[:, XFORM3W_ID])]
    for i in eachindex(center_bus)
        dict_bus[center_bus[i]] = size(bus, 1) + i
    end

    bus = vcat(bus, zeros(length(center_bus), size(bus, 2)))
    bus[end-length(center_bus)+1:end, BUS_I] = map(k -> dict_bus[k], center_bus)
    bus[end-length(center_bus)+1:end, BUS_TYPE] .= 1  # Set all center buses to type 1
    bus[end-length(center_bus)+1:end, VM] .= 1.0  # Set voltage magnitude to 1.0
    bus[end-length(center_bus)+1:end, VA] .= 0.0  # Set voltage angle to 0.0
    bus[end-length(center_bus)+1:end, BASEKV] .= 0.0  # Set base voltage to 0.0
    bus[end-length(center_bus)+1:end, ZONE] .= 1.0 # Set zone to 1.0
    bus[end-length(center_bus)+1:end, VMIN] .= 0.8  # Set min voltage to 0.95
    bus[end-length(center_bus)+1:end, VMAX] .= 1.05  # Set max voltage to 1.0

    #Create branch matrix
    branch_transformer_three = zeros(3*size(transformer_three_data, 1), 14)

    for i in eachindex(transformer_three_data[:, XFORM3W_ID])
        branch_transformer_three[3*i-2, FBUS] = dict_bus[transformer_three_data[i, XFORM3W_PRI_ELEMENT]]
        branch_transformer_three[3*i-2, TBUS] = dict_bus[center_bus[i]]
        branch_transformer_three[3*i-1, FBUS] = dict_bus[transformer_three_data[i, XFORM3W_SEC_ELEMENT]]
        branch_transformer_three[3*i-1, TBUS] = dict_bus[center_bus[i]]
        branch_transformer_three[3*i, FBUS] = dict_bus[transformer_three_data[i, XFORM3W_TERT_ELEMENT]]
        branch_transformer_three[3*i, TBUS] = dict_bus[center_bus[i]]
       
        branch_transformer_three[3*i-2, STATUS] = 1
        branch_transformer_three[3*i-1, STATUS] = 1
        branch_transformer_three[3*i, STATUS] = 1

        #Transfer string to Float64
        transformer_three_data[i, XFORM3W_PRI_KV] = parse.(Float64, transformer_three_data[i, XFORM3W_PRI_KV])
        transformer_three_data[i, XFORM3W_SEC_KV] = parse.(Float64, transformer_three_data[i, XFORM3W_SEC_KV])
        transformer_three_data[i, XFORM3W_TERT_KV] = parse.(Float64, transformer_three_data[i, XFORM3W_TERT_KV])

        transformer_three_data[i, XFORM3W_PRI_TAP] = parse.(Float64, transformer_three_data[i, XFORM3W_PRI_TAP])
        transformer_three_data[i, XFORM3W_SEC_TAP] = parse.(Float64, transformer_three_data[i, XFORM3W_SEC_TAP])
        transformer_three_data[i, XFORM3W_TERT_TAP] = parse.(Float64, transformer_three_data[i, XFORM3W_TERT_TAP])

        transformer_three_data[i, XFORM3W_PRI_KVA] = parse.(Float64, transformer_three_data[i, XFORM3W_PRI_KVA])
        transformer_three_data[i, XFORM3W_SEC_KVA] = parse.(Float64, transformer_three_data[i, XFORM3W_SEC_KVA])
        transformer_three_data[i, XFORM3W_TERT_KVA] = parse.(Float64, transformer_three_data[i, XFORM3W_TERT_KVA])

        transformer_three_data[i, XFORM3W_PRI_POSZ] = parse.(Float64, transformer_three_data[i, XFORM3W_PRI_POSZ])
        transformer_three_data[i, XFORM3W_PRI_X_R] = parse.(Float64, transformer_three_data[i, XFORM3W_PRI_X_R])
        transformer_three_data[i, XFORM3W_SEC_POSZ] = parse.(Float64, transformer_three_data[i, XFORM3W_SEC_POSZ])
        transformer_three_data[i, XFORM3W_SEC_X_R] = parse.(Float64, transformer_three_data[i, XFORM3W_SEC_X_R])
        transformer_three_data[i, XFORM3W_TERT_POSZ] = parse.(Float64, transformer_three_data[i, XFORM3W_TERT_POSZ])
        transformer_three_data[i, XFORM3W_TERT_X_R] = parse.(Float64, transformer_three_data[i, XFORM3W_TERT_X_R])

        # --- Δ to Y transformation ---
        r_ps = 0.01*transformer_three_data[i, XFORM3W_PRI_POSZ]/sqrt(1+transformer_three_data[i, XFORM3W_PRI_X_R]^2)
        x_ps = r_ps*transformer_three_data[i, XFORM3W_PRI_X_R]

        r_pt = 0.01*transformer_three_data[i, XFORM3W_SEC_POSZ]/sqrt(1+transformer_three_data[i, XFORM3W_SEC_X_R]^2)
        x_pt = r_pt*transformer_three_data[i, XFORM3W_SEC_X_R]

        r_st = 0.01*transformer_three_data[i, XFORM3W_TERT_POSZ]/sqrt(1+transformer_three_data[i, XFORM3W_TERT_X_R]^2)
        x_st = r_st*transformer_three_data[i, XFORM3W_TERT_X_R]

        Rp = (r_ps + r_pt - r_st)/2
        Xp = (x_ps + x_pt - x_st)/2
        Rs = (r_ps + r_st - r_pt)/2
        Xs = (x_ps + x_st - x_pt)/2
        Rt = (r_pt + r_st - r_ps)/2
        Xt = (x_pt + x_st - x_ps)/2

        branch_transformer_three[3*i-2, R] = Rp
        branch_transformer_three[3*i-2, X] = Xp
        branch_transformer_three[3*i-2, B] = 0.0

        branch_transformer_three[3*i-1, R] = Rs
        branch_transformer_three[3*i-1, X] = Xs
        branch_transformer_three[3*i-1, B] = 0.0

        branch_transformer_three[3*i, R] = Rt
        branch_transformer_three[3*i, X] = Xt
        branch_transformer_three[3*i, B] = 0.0

        # --- Normalize impedance to per unit ---
        transformer_baseS = transformer_three_data[i, XFORM3W_PRI_KVA] / 1000
        branch_transformer_three[3*i-2, R] = branch_transformer_three[3*i-2, R] * baseMVA / (transformer_baseS)
        branch_transformer_three[3*i-2, X] = branch_transformer_three[3*i-2, X] * baseMVA / (transformer_baseS)
        branch_transformer_three[3*i-2, B] = branch_transformer_three[3*i-2, B] / (baseMVA / (transformer_baseS))

        branch_transformer_three[3*i-1, R] = branch_transformer_three[3*i-1, R] * baseMVA / (transformer_baseS)
        branch_transformer_three[3*i-1, X] = branch_transformer_three[3*i-1, X] * baseMVA / (transformer_baseS)
        branch_transformer_three[3*i-1, B] = branch_transformer_three[3*i-1, B] / (baseMVA / (transformer_baseS))

        branch_transformer_three[3*i, R] = branch_transformer_three[3*i, R] * baseMVA / (transformer_baseS)
        branch_transformer_three[3*i, X] = branch_transformer_three[3*i, X] * baseMVA / (transformer_baseS) 
        branch_transformer_three[3*i, B] = branch_transformer_three[3*i, B] / (baseMVA / (transformer_baseS))

        branch_transformer_three[3*i-2, RATIO] =  1 + transformer_three_data[i, XFORM3W_PRI_TAP]/100
        branch_transformer_three[3*i-1, RATIO] = 1 + transformer_three_data[i, XFORM3W_SEC_TAP]/100
        branch_transformer_three[3*i, RATIO] = 1 + transformer_three_data[i, XFORM3W_TERT_TAP]/100
    end

    # --- Assign Rate and Status ---
    rate_values = 100  # MVA
    branch_transformer_three[:, RATEA] .= rate_values
    branch_transformer_three[:, RATEB] .= rate_values
    branch_transformer_three[:, RATEC] .= rate_values
    branch_transformer_three[:, ANGLE] .= 0.0
    branch_transformer_three[:, ANGMIN] .= -180.0
    branch_transformer_three[:, ANGMAX] .= 180.0

    branch_transformer_three[:, STATUS] .= 1.0
    branch_transformer_three[:, DICTKEY] .= 0.0

    return branch_transformer_three, bus, dict_bus
end