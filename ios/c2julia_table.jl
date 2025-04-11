using Libdl

# 加载DLL
if Sys.iswindows()
    const lib_path = joinpath(@__DIR__, "JuliaCall.dll")
else
    const lib_path = joinpath(@__DIR__, "libJuliaCall.so")
end

# 尝试加载库文件
let
    global lib
    try
        lib = Libdl.dlopen(lib_path)
        println("成功加载库: $lib_path")
    catch e
        error("无法加载库 $lib_path: $e")
    end
end

# 定义结构体
const RTDB_KEY_STRU = UInt64
primitive type SWITCH_STRU 112*8 end
const SWITCH_FIELD = "id,name,first_node_id,second_node_id,point,p,q,i_a_value,tpcolor"

# 从SWITCH_STRU中提取字段的函数
function extract_id(sw::SWITCH_STRU)
    bytes = reinterpret(UInt8, [sw])
    return reinterpret(UInt64, bytes[1:8])[1]
end

function extract_name(sw::SWITCH_STRU)
    bytes = reinterpret(UInt8, [sw])
    null_pos = findfirst(==(0x00), bytes[9:72])
    if null_pos === nothing
        return String(bytes[9:72])
    else
        return String(bytes[9:(8+null_pos-1)])
    end
end

function extract_ind(sw::SWITCH_STRU)
    bytes = reinterpret(UInt8, [sw])
    return reinterpret(Int64, bytes[73:80])[1]
end

function extract_jnd(sw::SWITCH_STRU)
    bytes = reinterpret(UInt8, [sw])
    return reinterpret(Int64, bytes[81:88])[1]
end

function extract_point(sw::SWITCH_STRU)
    bytes = reinterpret(UInt8, [sw])
    return bytes[89]
end

function extract_p(sw::SWITCH_STRU)
    bytes = reinterpret(UInt8, [sw])
    return reinterpret(Float32, bytes[90:93])[1]
end

function extract_q(sw::SWITCH_STRU)
    bytes = reinterpret(UInt8, [sw])
    return reinterpret(Float32, bytes[94:97])[1]
end

function extract_i(sw::SWITCH_STRU)
    bytes = reinterpret(UInt8, [sw])
    return reinterpret(Float32, bytes[98:101])[1]
end

function extract_tpcolor(sw::SWITCH_STRU)
    bytes = reinterpret(UInt8, [sw])
    return bytes[102]
end

# 只打印前3条记录的详细信息
function print_switch_record(sw::SWITCH_STRU, index::Int)
    id = extract_id(sw)
    name_str = extract_name(sw)
    ind = extract_ind(sw)
    jnd = extract_jnd(sw)
    point = extract_point(sw)
    p = extract_p(sw)
    q = extract_q(sw)
    i_val = extract_i(sw)
    tpcolor = extract_tpcolor(sw)
    
    println("记录 #$index: ID=$id, 名称=\"$name_str\", 首端=$ind, 末端=$jnd")
end

function read_switch_data_direct(app_no::Int, table_no::Int)
    println("尝试从C++接口获取数据...")
    
    # 确保readRTDB函数存在
    try
        sym = Libdl.dlsym(lib, :readRTDB)
    catch e
        error("无法找到readRTDB函数: $e")
    end
    
    # 获取记录数
    rec_num = Ref{Cint}(0)
    GC.@preserve SWITCH_FIELD begin
        ccall(Libdl.dlsym(lib, :readRTDB), Cvoid, 
              (Cint, Cint, Cstring, Ptr{Nothing}, Cint, Ptr{Cint}),
              app_no, table_no, SWITCH_FIELD, C_NULL, 0, rec_num)
    end
    
    println("C++函数返回的记录数: $(rec_num[])")
    expected_records = rec_num[] > 0 ? rec_num[] : 291147
    
    # 方法1：使用UInt8数组接收原始字节
    println("==== 方法1: 使用UInt8数组 ====")
    raw_buffer = Vector{UInt8}(undef, expected_records * sizeof(SWITCH_STRU))
    actual_count = Ref{Cint}(0)
    
    GC.@preserve SWITCH_FIELD raw_buffer begin
        ccall(Libdl.dlsym(lib, :readRTDB), Cvoid, 
              (Cint, Cint, Cstring, Ptr{UInt8}, Cint, Ptr{Cint}),
              app_no, table_no, SWITCH_FIELD, raw_buffer, expected_records, actual_count)
    end
    
    println("方法1返回的记录数: $(actual_count[])")
    nonzero_count = count(!=(0), raw_buffer)
    println("原始数据中非零字节数: $nonzero_count / $(length(raw_buffer))")
    
    if nonzero_count > 0
        # 尝试从原始字节中提取SWITCH_STRU结构
        switch_array = reinterpret(SWITCH_STRU, raw_buffer)
        
        # 只检查前几条记录
        valid_count = 0
        for i in 1:min(100, length(switch_array))
            id = extract_id(switch_array[i])
            if id != 0
                valid_count += 1
                if valid_count <= 3  # 只打印前3条记录
                    print_switch_record(switch_array[i], valid_count)
                end
                
                if valid_count >= 10
                    println("...")  # 表示还有更多记录
                    break
                end
            end
        end
        
        if valid_count > 0
            println("从原始字节中找到 $valid_count 条有效记录")
            return (method="原始字节", valid_records=switch_array, count=actual_count[])
        end
    end
    
    # 方法2：使用readSwitchData函数（如果存在）
    println("==== 方法2: 尝试使用readSwitchData函数 ====")
    try
        sym = Libdl.dlsym(lib, :readSwitchData)
        
        # 分配缓冲区
        result_buffer = Vector{SWITCH_STRU}(undef, expected_records)
        actual_count = Ref{Cint}(0)
        
        GC.@preserve result_buffer begin
            ccall(sym, Cvoid, 
                  (Cint, Cint, Ptr{SWITCH_STRU}, Cint, Ptr{Cint}),
                  app_no, table_no, result_buffer, expected_records, actual_count)
        end
        
        println("方法2返回的记录数: $(actual_count[])")
        
        # 检查数据有效性
        valid_count = 0
        valid_records = SWITCH_STRU[]
        
        for i in 1:expected_records
            id = extract_id(result_buffer[i])
            if id != 0
                valid_count += 1
                push!(valid_records, result_buffer[i])
                
                if valid_count <= 3  # 只打印前3条记录
                    print_switch_record(result_buffer[i], valid_count)
                end
                
                if valid_count >= 10
                    println("...")  # 表示还有更多记录
                    break
                end
            end
        end
        
        println("方法2找到有效记录数: $valid_count")
        
        if valid_count > 0
            return (method="readSwitchData", valid_records=valid_records, count=valid_count)
        end
    catch e
        println("readSwitchData函数不存在或调用失败")
    end
    
    # 方法3：使用不同的字段名
    println("==== 方法3: 使用不同的字段名 ====")
    simple_field = "id"
    result_buffer = Vector{UInt64}(undef, expected_records)
    actual_count = Ref{Cint}(0)
    
    GC.@preserve simple_field result_buffer begin
        ccall(Libdl.dlsym(lib, :readRTDB), Cvoid, 
              (Cint, Cint, Cstring, Ptr{UInt64}, Cint, Ptr{Cint}),
              app_no, table_no, simple_field, result_buffer, expected_records, actual_count)
    end
    
    println("方法3返回的记录数: $(actual_count[])")
    
    # 检查ID是否有效
    valid_count = count(id -> id != 0, result_buffer[1:min(expected_records, actual_count[])])
    
    if valid_count > 0
        println("前3个有效ID:")
        valid_ids = filter(id -> id != 0, result_buffer[1:min(expected_records, actual_count[])])
        for (i, id) in enumerate(valid_ids[1:min(3, length(valid_ids))])
            println("  ID #$i: $id")
        end
        if length(valid_ids) > 3
            println("  ...")
        end
        
        return (method="仅ID字段", valid_ids=valid_ids, count=valid_count)
    end
    
    # 如果所有方法都失败，返回失败结果
    return (method="全部失败", valid_records=SWITCH_STRU[], count=0)
end

function save_to_csv(switches::Vector{SWITCH_STRU}, filename::String)
    open(filename, "w") do io
        # 写入表头
        println(io, "ID,名称,首端,末端,点号,P,Q,I,拓扑颜色")
        
        # 写入数据
        for sw in switches
            id = extract_id(sw)
            name_str = extract_name(sw)
            ind = extract_ind(sw)
            jnd = extract_jnd(sw)
            point = extract_point(sw)
            p = extract_p(sw)
            q = extract_q(sw)
            i_val = extract_i(sw)
            tpcolor = extract_tpcolor(sw)
            
            println(io, "$id,\"$name_str\",$ind,$jnd,$point,$p,$q,$i_val,$tpcolor")
        end
    end
    
    println("数据已保存到文件: $filename")
end

function save_ids_to_csv(ids::Vector{UInt64}, filename::String)
    open(filename, "w") do io
        # 写入表头
        println(io, "ID")
        
        # 写入数据
        for id in ids
            println(io, "$id")
        end
    end
    
    println("ID数据已保存到文件: $filename")
end

function cleanup()
    Libdl.dlclose(lib)
    println("已关闭库")
end

# 主函数
function run_main()
    try
        # 测试参数
        app_no = 6500000
        table_no = 1706
        
        println("尝试直接读取数据")
        result = read_switch_data_direct(app_no, table_no)
        
        # 根据结果执行不同的操作
        if result.method == "原始字节" || result.method == "readSwitchData"
            println("成功: 使用【$(result.method)】读取到有效数据")
            println("找到 $(result.count) 条有效记录")
            
            # 保存数据到CSV文件
            save_to_csv(result.valid_records, "switch_data.csv")
        elseif result.method == "仅ID字段"
            println("部分成功: 使用【$(result.method)】读取到有效ID")
            println("找到 $(result.count) 条有效ID")
            
            # 保存ID到CSV文件
            save_ids_to_csv(result.valid_ids, "switch_ids.csv")
        else
            println("失败: 所有方法都未能读取到有效数据")
        end
        
    catch e
        println("运行过程中出错: $e")
    finally
        # 清理资源
        cleanup()
    end
end

# 运行主函数
run_main()
