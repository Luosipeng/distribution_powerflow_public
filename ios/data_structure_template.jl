
function create_distribution_system_empty_template()
    # 创建一个新的Excel工作簿
    xlsx_file = "配电系统数据录入模板.xlsx"
    XLSX.openxlsx(xlsx_file, mode="w") do xf
        # 创建指定的工作表，移除auto_trafo
        sheet_names = [
            "bus", "line", "line_dc", 
            "trafo", "trafo3w", "load", "asymmetric_load", "sgen", "asymmetric_sgen", "gen", "ext_grid", 
            "storage", "converter", "switch", 
            "equipment_carbon", "carbon_time_series", "carbon_scenario", 
            "vpp", 
            "mobile_storage", "charging_station", "charger", 
            "ev_aggregator", "v2g_service"
        ]
        
        # 创建说明工作表
        XLSX.addsheet!(xf, "说明")
        
        # 添加说明信息
        instructions = [
            "配电系统数据录入模板使用说明:",
            "",
            "1. 本模板包含多个工作表，每个工作表对应配电系统的一个组件类型",
            "2. 每个工作表的第1行是表头，包含字段名称",
            "3. 从第2行开始，按照表头要求填入数据",
            "4. 数据类型说明:",
            "   - 整数字段: index, bus编号等",
            "   - 浮点数字段: 电压, 功率, 阻抗等",
            "   - 文本字段: 名称, 类型等",
            "   - 布尔值字段: in_service等，填写TRUE/FALSE",
            "",
            "5. 每个组件必须有唯一的index作为标识",
            "6. 对于引用其他组件的字段(如bus引用)，请确保引用的组件已存在",
            "7. bus表中，is_dc字段表示是否为直流节点，TRUE为直流节点，FALSE为交流节点",
            "8. vpp表包含了虚拟电厂基本信息、资源和聚合负荷的数据",
            "9. gen表包含常规发电机组的数据",
            "10. trafo表包含两卷变压器数据，trafo3w表包含三卷变压器的数据",
            "",
            "11. 完成数据填写后，可以使用相应的软件导入此Excel文件进行配电系统分析"
        ]
        
        for (i, line) in enumerate(instructions)
            XLSX.setdata!(xf["说明"], "A$i", line)
        end
        
        # 创建其他工作表
        for name in sheet_names
            XLSX.addsheet!(xf, name)
        end

        # 1. 合并后的节点表(Bus)，包含AC和DC节点
        bus_headers = ["index", "name", "type", "zone", "vn_kv", "is_dc", "max_vm_pu", "min_vm_pu", "in_service"]
        
        # 写入表头（字段名称）
        for (col, header) in enumerate(bus_headers)
            XLSX.setdata!(xf["bus"], XLSX.CellRef(1, col), header)
        end
        
        # 2. 交流线路(Line)
        line_headers = [
            "index", "name", "from_bus", "to_bus", "length_km", 
            "r_ohm_per_km", "x_ohm_per_km", "c_nf_per_km",
            "r0_ohm_per_km", "x0_ohm_per_km", "c0_nf_per_km", 
            "g_us_per_km", "max_i_ka", "type", "max_loading_percent", 
            "parallel", "df", "in_service", "mtbf_hours", "mttr_hours", 
            "failure_rate_per_year", "planned_outage_hours_per_year", 
            "forced_outage_rate", "permanent_fault_rate_per_km_year", 
            "temporary_fault_rate_per_km_year", "repair_time_permanent_hours", 
            "auto_reclosing_success_rate"
        ]
        
        # 写入表头
        for (col, header) in enumerate(line_headers)
            XLSX.setdata!(xf["line"], XLSX.CellRef(1, col), header)
        end
        
        # 3. 直流线路(Line DC)
        line_dc_headers = [
            "index", "name", "from_bus", "to_bus", "length_km", 
            "r_ohm_per_km", "g_us_per_km", "max_i_ka", "type", 
            "max_loading_percent", "parallel", "df", "in_service"
        ]
        
        # 写入表头
        for (col, header) in enumerate(line_dc_headers)
            XLSX.setdata!(xf["line_dc"], XLSX.CellRef(1, col), header)
        end
        
        # 4. 两卷变压器(Transformer 2W)
        trafo_headers = [
            "index", "name", "std_type", "hv_bus", "lv_bus", "sn_mva", 
            "vn_hv_kv", "vn_lv_kv", "vk_percent", "vkr_percent", "pfe_kw", 
            "i0_percent", "vk0_percent", "vkr0_percent", "mag0_percent", "mag0_rx", "si0_hv_partial", "shift_degree", "tap_side", "tap_neutral", 
            "tap_min", "tap_max", "tap_step_percent", "tap_step_degree", 
            "tap_pos", "tap_phase_shifter", "parallel", "max_loading_percent", 
            "df", "in_service", "tap_changer_type", "oltc", "power_station_unit", "vector_group", 
            "hv_connection", "lv_connection", "thermal_capacity_mw", 
            "cooling_type", "oil_volume_liters", "winding_material"
        ]
        
        # 写入表头
        for (col, header) in enumerate(trafo_headers)
            XLSX.setdata!(xf["trafo"], XLSX.CellRef(1, col), header)
        end

        # 5. 三卷变压器(Transformer 3W)
        trafo3w_headers = [
            "index", "name", "std_type", "hv_bus", "mv_bus", "lv_bus", "sn_hv_mva", "sn_mv_mva", "sn_lv_mva", 
            "vn_hv_kv", "vn_mv_kv", "vn_lv_kv", "vk_hv_percent", "vk_mv_percent", "vk_lv_percent", 
            "vkr_hv_percent", "vkr_mv_percent", "vkr_lv_percent", "pfe_kw", "i0_percent", "shift_mv_degree", 
            "shift_lv_degree", "tap_side", "tap_neutral", "tap_min", "tap_max", "tap_step_percent", 
            "tap_step_degree", "tap_at_star_point", "tap_pos", "in_service", "vector_group_hv_mv", 
            "vector_group_hv_lv", "vector_group_mv_lv", "hv_connection", "mv_connection", "lv_connection", 
            "thermal_capacity_mw", "cooling_type", "oil_volume_liters", "winding_material"
        ]
        
        # 写入表头
        for (col, header) in enumerate(trafo3w_headers)
            XLSX.setdata!(xf["trafo3w"], XLSX.CellRef(1, col), header)
        end
        
        # 6. 负荷(Load)
        load_headers = [
            "index", "name", "bus", "p_mw", "q_mvar", "const_z_percent", 
            "const_i_percent", "const_p_percent", "scaling", "in_service", "type"
        ]
        
        
        # 写入表头
        for (col, header) in enumerate(load_headers)
            XLSX.setdata!(xf["load"], XLSX.CellRef(1, col), header)
        end
        # 7. 不对称负荷(asymmetric_load)
        asymmetric_load_headers = [
            "index", "name", "bus", "p_a_mw", "p_b_mw", "p_c_mw", 
            "q_a_mvar", "q_b_mvar", "q_c_mvar", "sn_mva", "scaling", 
            "in_service", "type"
        ]
        # 负荷表头包含两种类型：对称负荷和不对称负荷

        for (col, header) in enumerate(asymmetric_load_headers)
            XLSX.setdata!(xf["asymmetric_load"], XLSX.CellRef(1, col), header)
        end
        # 8. 静态发电机(Static Generator)
        sgen_headers = [
            "index", "name", "bus", "p_mw", "q_mvar", "scaling", 
            "max_p_mw", "min_p_mw", "max_q_mvar", "min_q_mvar", 
            "k", "rx", "in_service", "type", "controllable"
        ]
        
        # 写入表头
        for (col, header) in enumerate(sgen_headers)
            XLSX.setdata!(xf["sgen"], XLSX.CellRef(1, col), header)
        end
        
        # 9. 不对称静态发电机(asymmetric_sgen)
        asymmetric_sgen_headers = [
            "index", "name", "bus", "p_a_mw", "p_b_mw", "p_c_mw", 
            "q_a_mvar", "q_b_mvar", "q_c_mvar", "sn_mva", 
            "scaling", "in_service", "type", "current_source"
        ]

        # 写入表头
        for (col, header) in enumerate(asymmetric_sgen_headers)
            XLSX.setdata!(xf["asymmetric_sgen"], XLSX.CellRef(1, col), header)
        end
        # 8. 常规发电机组(Generator)
        gen_headers = [
            "index", "name", "bus", "p_mw", "vm_pu", "sn_mva", "scaling", 
            "max_p_mw", "min_p_mw", "max_q_mvar", "min_q_mvar", "vn_kv", 
            "xdss_pu", "rdss_pu", "cos_phi", "controllable", "in_service", 
            "type", "generator_type", "fuel_type", "startup_time_cold_h", 
            "startup_time_warm_h", "startup_time_hot_h", "min_up_time_h", 
            "min_down_time_h", "ramp_up_rate_mw_per_min", "ramp_down_rate_mw_per_min", 
            "efficiency_percent", "heat_rate_mmbtu_per_mwh", "co2_emission_rate_kg_per_mwh"
        ]
        
        # 写入表头
        for (col, header) in enumerate(gen_headers)
            XLSX.setdata!(xf["gen"], XLSX.CellRef(1, col), header)
        end
        
        # 9. 外部电网(External Grid)
        ext_grid_headers = [
            "index", "name", "bus", "vm_pu", "va_degree", "in_service", 
            "s_sc_max_mva", "s_sc_min_mva", "rx_max", "rx_min", "r0x0_max", 
            "x0x_max", "controllable"
        ]
        
        # 写入表头
        for (col, header) in enumerate(ext_grid_headers)
            XLSX.setdata!(xf["ext_grid"], XLSX.CellRef(1, col), header)
        end
        
        # 10. 储能(Storage)
        storage_headers = [
            "index", "name", "bus", "p_mw", "q_mvar", "sn_mva", 
            "soc_percent", "min_e_mwh", "max_e_mwh", "charge_efficiency_percent", 
            "discharge_efficiency_percent", "max_p_mw", "min_p_mw", "max_q_mvar", 
            "min_q_mvar", "in_service", "type", "controllable"
        ]
        
        # 写入表头
        for (col, header) in enumerate(storage_headers)
            XLSX.setdata!(xf["storage"], XLSX.CellRef(1, col), header)
        end
        
        # 11. 换流器(Converter) - 原vsc
        converter_headers = [
            "index", "name", "bus_ac", "bus_dc", "p_mw", "q_mvar", 
            "vm_ac_pu", "vm_dc_pu", "loss_percent", "loss_mw", "max_p_mw", 
            "min_p_mw", "max_q_mvar", "min_q_mvar", "control_mode", 
            "droop_kv", "in_service", "controllable"
        ]
        
        # 写入表头
        for (col, header) in enumerate(converter_headers)
            XLSX.setdata!(xf["converter"], XLSX.CellRef(1, col), header)
        end
        
        # 12. 开关(Switch)
        switch_headers = [
            "index", "name", "bus_from", "bus_to", "element_type", "element_id", 
            "closed", "type", "z_ohm", "in_service"
        ]
        
        # 写入表头
        for (col, header) in enumerate(switch_headers)
            XLSX.setdata!(xf["switch"], XLSX.CellRef(1, col), header)
        end
        
        # 13. 设备碳排放(Equipment Carbon)
        equipment_carbon_headers = [
            "index", "element_type", "element_id", "carbon_embodied_kgCO2e", 
            "carbon_operational_kgCO2e_per_year", "lifetime_years", 
            "manufacturing_date", "installation_date", "recycling_rate_percent"
        ]
        
        # 写入表头
        for (col, header) in enumerate(equipment_carbon_headers)
            XLSX.setdata!(xf["equipment_carbon"], XLSX.CellRef(1, col), header)
        end
        
        # 14. 碳排放时间序列(Carbon Time Series)
        carbon_time_series_headers = [
            "index", "timestamp", "grid_carbon_intensity_kgCO2e_per_MWh", 
            "renewable_generation_carbon_intensity_kgCO2e_per_MWh", 
            "storage_carbon_intensity_kgCO2e_per_MWh"
        ]
        
        # 写入表头
        for (col, header) in enumerate(carbon_time_series_headers)
            XLSX.setdata!(xf["carbon_time_series"], XLSX.CellRef(1, col), header)
        end
        
        # 15. 碳排放情景(Carbon Scenario)
        carbon_scenario_headers = [
            "index", "name", "description", "year", 
            "grid_carbon_intensity_kgCO2e_per_MWh", 
            "renewable_penetration_percent", "ev_adoption_percent", 
            "storage_adoption_percent"
        ]
        
        # 写入表头
        for (col, header) in enumerate(carbon_scenario_headers)
            XLSX.setdata!(xf["carbon_scenario"], XLSX.CellRef(1, col), header)
        end
        
        # 16. 虚拟电厂(VPP) - 合并了vpp, vpp_resources, vpp_aggregated_load
        vpp_headers = [
            # 基本信息 - 来自原vpp表
            "index", "name", "description", "control_area", "capacity_mw", 
            "energy_mwh", "response_time_s", "ramp_rate_mw_per_min", 
            "availability_percent", "operator", "in_service",
            
            # 资源信息 - 来自原vpp_resources表
            "resource_type", "resource_id", "capacity_share_percent", 
            "control_priority", "resource_response_time_s", "max_duration_h",
            
            # 负荷信息 - 来自原vpp_aggregated_load表
            "timestamp", "p_mw", "q_mvar", "flexibility_up_mw", 
            "flexibility_down_mw", "flexibility_duration_h"
        ]
        
        # 写入表头
        for (col, header) in enumerate(vpp_headers)
            XLSX.setdata!(xf["vpp"], XLSX.CellRef(1, col), header)
        end
        
        # 17. 移动储能(Mobile Storage)
        mobile_storage_headers = [
            "index", "name", "type", "capacity_mwh", "power_mw", 
            "soc_percent", "charge_efficiency_percent", 
            "discharge_efficiency_percent", "owner", "availability", 
            "current_location", "status"
        ]
        
        # 写入表头
        for (col, header) in enumerate(mobile_storage_headers)
            XLSX.setdata!(xf["mobile_storage"], XLSX.CellRef(1, col), header)
        end
        
        # 18. 充电站(Charging Station)
        charging_station_headers = [
            "index", "name", "bus", "location", "operator", 
            "num_chargers", "max_power_kw", "in_service"
        ]
        
        # 写入表头
        for (col, header) in enumerate(charging_station_headers)
            XLSX.setdata!(xf["charging_station"], XLSX.CellRef(1, col), header)
        end
        
        # 19. 充电桩(Charger)
        charger_headers = [
            "index", "name", "station_id", "type", "max_power_kw", 
            "min_power_kw", "connector_type", "v2g_capable", 
            "in_service", "availability"
        ]
        
        # 写入表头
        for (col, header) in enumerate(charger_headers)
            XLSX.setdata!(xf["charger"], XLSX.CellRef(1, col), header)
        end
        
        # 20. 电动汽车聚合商(EV Aggregator)
        ev_aggregator_headers = [
            "index", "name", "num_evs", "total_capacity_mwh", 
            "max_power_mw", "control_strategy", "service_area", "operator"
        ]
        
        # 写入表头
        for (col, header) in enumerate(ev_aggregator_headers)
            XLSX.setdata!(xf["ev_aggregator"], XLSX.CellRef(1, col), header)
        end
        
        # 21. 车网互动服务(V2G Service)
        v2g_service_headers = [
            "index", "aggregator_id", "service_type", "start_time", 
            "end_time", "capacity_mw", "energy_mwh", "price_per_mwh", 
            "grid_area", "reliability_percent"
        ]
        
        # 写入表头
        for (col, header) in enumerate(v2g_service_headers)
            XLSX.setdata!(xf["v2g_service"], XLSX.CellRef(1, col), header)
        end
    end
    
    println("Excel模板已创建: 配电系统数据录入模板.xlsx")
    println("请打开Excel文件，按照表头要求在每个工作表中从第2行开始填写数据")
    return xlsx_file
end

# 添加一个辅助函数，用于导入数据
function import_distribution_system_data(xlsx_file)
    data = Dict()
    
    XLSX.openxlsx(xlsx_file) do xf
        # 获取所有工作表名称
        sheet_names = XLSX.sheetnames(xf)
        
        # 排除说明工作表
        filter!(name -> name != "Sheet1", sheet_names)
        filter!(name -> name != "说明", sheet_names)
        
        # 导入每个工作表的数据
        for sheet_name in sheet_names
            # 读取工作表数据，从第2行开始（跳过表头行）
            # df = DataFrame(XLSX.readtable(xf[sheet_name], first_row=2))
            df = DataFrame(XLSX.readtable(xlsx_file, sheet_name, first_row=1))
            
            if !isempty(df)
                # 获取字段名称（第1行）
                headers = []
                for col in 1:XLSX.get_dimension(xf[sheet_name]).stop.column_number
                    header_value = XLSX.getdata(xf[sheet_name], XLSX.CellRef(1, col))
                    if header_value !== nothing
                        push!(headers, string(header_value))
                    else
                        break
                    end
                end
                
                # 重命名列名为实际的字段名
                # rename_dict = Dict(Symbol("x$i") => Symbol(headers[i]) for i in 1:length(headers))
                # df = rename(df, rename_dict)
                
                # 存储处理后的数据
                data[sheet_name] = df
            end
        end
    end
    
    return data
end

# 执行函数创建Excel模板
# template_file = create_distribution_system_empty_template()
# println("您可以在 $template_file 中按行填写数据")