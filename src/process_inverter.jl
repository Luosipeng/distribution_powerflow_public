function process_inverter(mpc)
    #读取索引
    (LOAD_I,LOAD_CND,LOAD_STATUS,LOAD_PD,LOAD_QD,LOADZ_PERCENT,LOADI_PERCENT,LOADP_PERCENT)=PowerFlow.idx_ld()
    (PQ, PV, REF, NONE, BUS_I, BUS_TYPE, PD, QD, GS, BS, BUS_AREA, VM,VA, 
    BASE_KV, ZONE, VMAX, VMIN, LAM_P, LAM_Q, MU_VMAX, MU_VMIN, PER_CONSUMER) = PowerFlow.idx_bus();

    #复制mpc
    mpc_result = deepcopy(mpc)
    #读取数据
    Cr = mpc_result["Cr"]
    Ci = mpc_result["Ci"]
    P_inv = mpc_result["P_inv"]
    Q_inv = mpc_result["Q_inv"]
    P_inv_dc = mpc_result["P_inv_dc"]
    busAC = mpc_result["busAC"]
    busDC = mpc_result["busDC"]
    if haskey(mpc_result, "loadAC") == false
        loadAC = zeros(size(busAC, 1), 8)
        loadAC[:, LOAD_I] = collect(1:size(busAC, 1))
        loadAC[:, LOAD_CND] = busAC[:, BUS_I]
        loadAC[:, LOAD_STATUS] .= 1
        loadAC[:, LOAD_PD] .= busAC[:, PD]
        loadAC[:, LOAD_QD] .= busAC[:, QD]
        loadAC[:, LOADZ_PERCENT] .= 0.0
        loadAC[:, LOADI_PERCENT] .= 0.0
        loadAC[:, LOADP_PERCENT] .= 1.0
        loadDC = zeros(size(busDC, 1), 8)
        loadDC[:, LOAD_I] = collect(1:size(busDC, 1))
        loadDC[:, LOAD_CND] = busDC[:, BUS_I]
        loadDC[:, LOAD_STATUS] .= 1
        loadDC[:, LOAD_PD] .= busDC[:, PD]
        loadDC[:, LOAD_QD] .= busDC[:, QD]
        loadDC[:, LOADZ_PERCENT] .= 0.0
        loadDC[:, LOADI_PERCENT] .= 0.0
        loadDC[:, LOADP_PERCENT] .= 1.0
    else
        loadAC = mpc_result["loadAC"]
        loadDC = mpc_result["loadDC"]
    end

    #处理交流逆变器
    # #处理逆变器传输功率
    if P_inv!==nothing
        inverter_connected_bus=Ci
        for i in inverter_connected_bus
            if(i in loadAC[:,LOAD_CND])
                idx=findall(x->x==i,loadAC[:,LOAD_CND])
                loadAC[idx,LOADP_PERCENT]=(loadAC[idx,LOAD_PD].*loadAC[idx,LOADP_PERCENT].-P_inv[inverter_connected_bus.==i]./1000)./(loadAC[idx,LOAD_PD].-P_inv[inverter_connected_bus.==i]./1000)
                # load[idx,LOADI_PERCENT]=(load[idx,LOAD_PD].*load[idx,LOADI_PERCENT].+P_inv[inverter_connected_bus.==i]./1000)./(load[idx,LOAD_PD].+P_inv[inverter_connected_bus.==i]./1000)
                loadAC[idx,LOADI_PERCENT].=0
                loadAC[idx,LOADZ_PERCENT]=1 .-loadAC[idx,LOADP_PERCENT].-loadAC[idx,LOADI_PERCENT]

                loadAC[idx,LOAD_PD].-=P_inv[inverter_connected_bus.==i]./1000
                loadAC[idx,LOAD_QD].-=Q_inv[inverter_connected_bus.==i]./1000
                
                busAC[busAC[:,BUS_I].==i,PD] .= loadAC[idx,LOAD_PD]
                busAC[busAC[:,BUS_I].==i,QD] .= loadAC[idx,LOAD_QD]

            else
                load_inv=zeros(1,8)
                load_inv[LOAD_I]=size(loadAC,1)+1
                load_inv[LOAD_CND]=i
                load_inv[LOAD_STATUS]=1
                load_inv[:,LOAD_PD].=-P_inv[inverter_connected_bus.==i]./1000
                load_inv[:,LOAD_QD].=-Q_inv[inverter_connected_bus.==i]./1000
                load_inv[LOADZ_PERCENT]=0
                load_inv[LOADI_PERCENT]=0
                load_inv[LOADP_PERCENT]=1.0
                loadAC=vcat(loadAC,load_inv)
                busAC[busAC[:,BUS_I].==i,PD] .= load_inv[LOAD_PD]
                busAC[busAC[:,BUS_I].==i,QD] .= load_inv[LOAD_QD]
            end
        end
    end

    #处理直流整流器
    if(P_inv_dc!==nothing)
        #Find the connected bus for the inverter
        inverter_index=Cr
        for i in inverter_index
            if(i in loadDC[:,LOAD_CND])
                idx=findall(x->x==i,loadDC[:,LOAD_CND])
                loadDC[idx,LOADP_PERCENT].=(loadDC[idx,LOAD_PD].*loadDC[idx,LOADP_PERCENT].+P_inv_dc[inverter_index.==i]./1000)./(loadDC[idx,LOAD_PD].+P_inv_dc[inverter_index.==i]./1000)
                loadDC[idx,LOADI_PERCENT].=0
                loadDC[idx,LOADZ_PERCENT].=1 .-loadDC[idx,LOADP_PERCENT].-loadDC[idx,LOADI_PERCENT]

                loadDC[idx,LOAD_PD].+=P_inv_dc[inverter_index.==i]./1000
                loadDC[idx,LOAD_QD].+=0

                busDC[busDC[:,BUS_I].==i,PD] .= loadDC[idx,LOAD_PD]
                busDC[busDC[:,BUS_I].==i,QD] .= loadDC[idx,LOAD_QD]

            else
                P_add = P_inv_dc[inverter_index.==i][1]./1000
                busDC[busDC[:,BUS_I].==i,PD] .= P_add
                busDC[busDC[:,BUS_I].==i,QD] .= 0
                
                # 创建新行并添加到 dcload
                new_row = [size(loadDC,1)+1, i, 1, P_add, 0, 0, 0, 1]'
                loadDC = vcat(loadDC, new_row)
            end
        end
    end
    mpc_result["busAC"] = busAC
    mpc_result["busDC"] = busDC
    mpc_result["loadAC"] = loadAC
    mpc_result["loadDC"] = loadDC
    # mpc_result["genAC"][1,6] = busAC[1,8]
    # mpc_result["genDC"][1,6] = busDC[7,8]
    return mpc_result
end