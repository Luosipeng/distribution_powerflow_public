using Libdl
using Printf
using Dates

# 加载C共享库
const libpath = "C:/Users/DELL/source/repos/data_import/data_import/libpowersystem.dll"  # 替换为实际的DLL路径
const lib = Libdl.dlopen(libpath)

# 定义数据结构
struct Bus
    bus_id::Int32
    voltage::Float64
    angle::Float64
    p_load::Float64
    q_load::Float64
    last_update::NTuple{20, UInt8}
end

struct Generator
    gen_id::Int32
    bus_id::Int32
    p_gen::Float64
    q_gen::Float64
    v_setpoint::Float64
    last_update::NTuple{20, UInt8}
end

struct Branch
    branch_id::Int32
    from_bus::Int32
    to_bus::Int32
    resistance::Float64
    reactance::Float64
    charging::Float64
    last_update::NTuple{20, UInt8}
end

# 连接到数据库
function connect_to_db()
    f = Libdl.dlsym(lib, :connect_to_db)
    result = ccall(f, Int32, ())
    if result == 0
        error("无法连接到数据库")
    end
    return result
end

# 断开数据库连接
function disconnect_db()
    f = Libdl.dlsym(lib, :disconnect_db)
    ccall(f, Cvoid, ())
    println("数据库连接已关闭")
end

# 获取母线数量
function get_buses_count()
    f = Libdl.dlsym(lib, :get_buses_count)
    return ccall(f, Int32, ())
end

# 获取发电机数量
function get_generators_count()
    f = Libdl.dlsym(lib, :get_generators_count)
    return ccall(f, Int32, ())
end

# 获取支路数量
function get_branches_count()
    f = Libdl.dlsym(lib, :get_branches_count)
    return ccall(f, Int32, ())
end

# 获取母线数据
function get_buses()
    count = get_buses_count()
    if count == 0
        return Bus[]
    end
    
    buses = Vector{Bus}(undef, count)
    f = Libdl.dlsym(lib, :get_buses)
    actual_count = ccall(f, Int32, (Ptr{Bus}, Int32), buses, count)
    
    return buses[1:actual_count]
end

# 获取发电机数据
function get_generators()
    count = get_generators_count()
    if count == 0
        return Generator[]
    end
    
    generators = Vector{Generator}(undef, count)
    f = Libdl.dlsym(lib, :get_generators)
    actual_count = ccall(f, Int32, (Ptr{Generator}, Int32), generators, count)
    
    return generators[1:actual_count]
end

# 获取支路数据
function get_branches()
    count = get_branches_count()
    if count == 0
        return Branch[]
    end
    
    branches = Vector{Branch}(undef, count)
    f = Libdl.dlsym(lib, :get_branches)
    actual_count = ccall(f, Int32, (Ptr{Branch}, Int32), branches, count)
    
    return branches[1:actual_count]
end

# 获取总负荷
function get_total_load()
    f = Libdl.dlsym(lib, :get_total_load)
    return ccall(f, Float64, ())
end

# 获取总发电量
function get_total_generation()
    f = Libdl.dlsym(lib, :get_total_generation)
    return ccall(f, Float64, ())
end

# 将字节数组转换为字符串
function bytes_to_string(bytes)
    idx = findfirst(isequal(0), bytes)
    if idx === nothing
        return String(UInt8[bytes...])
    else
        return String(UInt8[bytes[1:idx-1]...])
    end
end

# 打印母线数据
function print_buses(buses)
    println("\n--- 母线数据 ---")
    println("母线ID  电压(pu)  相角(rad)  有功负荷  无功负荷  最后更新时间")
    println("------------------------------------------------------------------")
    for bus in buses
        last_update = bytes_to_string(bus.last_update)
        @printf("%6d  %8.2f  %9.2f  %8.2f  %8.2f  %s\n", 
                bus.bus_id, bus.voltage, bus.angle, bus.p_load, bus.q_load, last_update)
    end
end

# 打印发电机数据
function print_generators(generators)
    println("\n--- 发电机数据 ---")
    println("发电机ID  母线ID  有功出力  无功出力  电压设定  最后更新时间")
    println("------------------------------------------------------------------")
    for gen in generators
        last_update = bytes_to_string(gen.last_update)
        @printf("%9d  %6d  %8.2f  %8.2f  %8.2f  %s\n", 
                gen.gen_id, gen.bus_id, gen.p_gen, gen.q_gen, gen.v_setpoint, last_update)
    end
end

# 打印支路数据
function print_branches(branches)
    println("\n--- 支路数据 ---")
    println("支路ID  起始点  终止点  电阻    电抗    充电    最后更新时间")
    println("------------------------------------------------------------------")
    for branch in branches
        last_update = bytes_to_string(branch.last_update)
        @printf("%7d  %6d  %6d  %6.2f  %6.2f  %6.2f  %s\n", 
                branch.branch_id, branch.from_bus, branch.to_bus, 
                branch.resistance, branch.reactance, branch.charging, last_update)
    end
end

# 打印系统摘要
function print_system_summary()
    total_load = get_total_load()
    total_generation = get_total_generation()
    
    println("\n--- 系统摘要 ---")
    println("母线数量: $(get_buses_count())")
    println("发电机数量: $(get_generators_count())")
    println("支路数量: $(get_branches_count())")
    println("总负荷: $(round(total_load, digits=2)) MW")
    println("总发电量: $(round(total_generation, digits=2)) MW")
    
    if total_generation > 0
        loss = total_generation - total_load
        loss_percent = (loss / total_generation) * 100
        println("系统损耗: $(round(loss, digits=2)) MW ($(round(loss_percent, digits=2))%)")
    end
end

# 读取电力系统数据
function read_power_system_data()
    try
        # 连接数据库
        connect_to_db()
        println("成功连接到电力系统数据库")
        
        # 获取数据
        println("正在读取数据...")
        buses = get_buses()
        generators = get_generators()
        branches = get_branches()
        
        # 打印数据
        print_buses(buses)
        print_generators(generators)
        print_branches(branches)
        print_system_summary()
        
        # 返回数据供进一步处理
        return buses, generators, branches
        
    catch e
        println("错误: $e")
        return nothing, nothing, nothing
    finally
        # 确保断开数据库连接
        disconnect_db()
    end
end

# 直接运行主函数
buses, generators, branches=read_power_system_data()
