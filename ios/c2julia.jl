using Libdl
using JLD2

# 加载DLL
# 根据操作系统选择正确的库文件名
if Sys.iswindows()
    const lib_path = joinpath(@__DIR__, "JuliaCall.dll")
else
    const lib_path = joinpath(@__DIR__, "libJuliaCall.so")
end

# 尝试加载库文件
try
    global lib = Libdl.dlopen(lib_path)
    println("成功加载库: $lib_path")
catch e
    error("无法加载库 $lib_path: $e")
end

"""
    saveCharData(data::Vector{UInt8}, filename::String)

将二进制数据保存到文件。
"""
# function saveCharData(data::Vector{UInt8}, filename::String)
#     ccall(Libdl.dlsym(lib, :saveCharData), Cvoid, 
#           (Ptr{UInt8}, Csize_t, Cstring), 
#           data, length(data), filename)
# end

# """
#     loadCharData(filename::String)

# 从文件加载二进制数据。
# """
# function loadCharData(filename::String)
#     loaded_length = Ref{Csize_t}(0)
#     data_ptr = ccall(Libdl.dlsym(lib, :loadCharData), Ptr{UInt8}, 
#                     (Cstring, Ptr{Csize_t}), 
#                     filename, loaded_length)
    
#     if data_ptr == C_NULL
#         return UInt8[]
#     end
    
#     # 将数据复制到Julia数组
#     data = unsafe_wrap(Array, data_ptr, loaded_length[], own=false)
#     result = copy(data)  # 创建一个副本
    
#     # 释放C++分配的内存
#     ccall(Libdl.dlsym(lib, :freeCharData), Cvoid, (Ptr{UInt8},), data_ptr)
    
#     return result
# end

# """
#     getFormattedDate()

# 获取格式化的当前日期（yyyyMMdd格式）。
# """
# function getFormattedDate()
#     date_ptr = ccall(Libdl.dlsym(lib, :getFormattedDate), Cstring, ())
#     return unsafe_string(date_ptr)
# end

"""
    createDirectoryIfNotExists(path::String)

如果目录不存在，则创建目录。
"""
# function createDirectoryIfNotExists(path::String)
#     result = ccall(Libdl.dlsym(lib, :createDirectoryIfNotExists), Cint, 
#                   (Cstring,), path)
#     return result != 0
# end

"""
    readRTDB(app_no::Int, table_no::Int, column_name::String, buffer_size::Int=0)

从RTDB读取数据。
如果buffer_size为0，则自动确定所需的缓冲区大小。
返回UInt64数组。
"""
function readRTDB(app_no::Int, table_no::Int, column_name::String, buffer_size::Int=0)
    if buffer_size <= 0
        # 第一次调用，获取记录总数
        rec_num = Ref{Cint}(0)
        ccall(Libdl.dlsym(lib, :readRTDB), Cvoid, 
              (Cint, Cint, Cstring, Ptr{UInt64}, Cint, Ptr{Cint}),
              app_no, table_no, column_name, C_NULL, 0, rec_num)
        
        total_records = rec_num[]
        if total_records <= 0
            println("没有找到记录")
            return UInt64[]
        end
        
        println("获取设备信息成功，记录数:$(total_records)")
        
        # 第二次调用，分配足够大的缓冲区
        buffer_size = total_records
        result_buffer = zeros(UInt64, buffer_size)
        rec_num[] = 0  # 重置记录数
        
        ccall(Libdl.dlsym(lib, :readRTDB), Cvoid, 
              (Cint, Cint, Cstring, Ptr{UInt64}, Cint, Ptr{Cint}),
              app_no, table_no, column_name, result_buffer, buffer_size, rec_num)
        
        actual_count = rec_num[]
        println("设备向量 size:$(actual_count)")
        
        if actual_count <= 0
            return UInt64[]
        end
        
        if actual_count > buffer_size
            @warn "返回的记录数($actual_count)超出了缓冲区大小($buffer_size)，将被截断"
            actual_count = buffer_size
        end
        
        return result_buffer[1:actual_count]
    else
        # 使用指定的缓冲区大小
        rec_num = Ref{Cint}(0)
        result_buffer = zeros(UInt64, buffer_size)
        
        ccall(Libdl.dlsym(lib, :readRTDB), Cvoid, 
              (Cint, Cint, Cstring, Ptr{UInt64}, Cint, Ptr{Cint}),
              app_no, table_no, column_name, result_buffer, buffer_size, rec_num)
        
        actual_count = rec_num[]
        println("获取设备信息成功，记录数:$(actual_count)")
        println("设备向量 size:$(actual_count)")
        
        if actual_count <= 0
            return UInt64[]
        end
        
        if actual_count > buffer_size
            @warn "返回的记录数($actual_count)超出了缓冲区大小($buffer_size)，将被截断"
            actual_count = buffer_size
        end
        
        return result_buffer[1:actual_count]
    end
end


"""
    cleanup()

关闭动态库，释放资源。
"""
function cleanup()
    Libdl.dlclose(lib)
    println("已关闭库")
end

# 使用示例
if abspath(PROGRAM_FILE) == @__FILE__
    try
        # # 创建目录
        # dir_created = createDirectoryIfNotExists("data_output")
        # println("目录创建状态: $dir_created")
        
        # # 获取当前日期
        # current_date = getFormattedDate()
        # println("当前日期: $current_date")
        
        # # 保存一些测试数据
        # test_data = Vector{UInt8}("这是一个测试数据")
        # saveCharData(test_data, "data_output/test.bin")
        
        # # 读取测试数据
        # loaded_data = loadCharData("data_output/test.bin")
        # println("读取的数据: ", String(loaded_data))
        
        # 从RTDB读取数据
        app_no = 6500000    # 根据实际情况调整
        table_no = 1706  # 根据实际情况调整
        column_name = "id"  # 根据实际情况调整
        
        rtdb_data = readRTDB(app_no, table_no, column_name)
        @save "rtdb_data.jld2" rtdb_data
        println("RTDB数据 ($(length(rtdb_data))条记录):")
        if !isempty(rtdb_data)
            for (i, value) in enumerate(rtdb_data)
                if i <= 10  # 只显示前10条
                    println("  $i: $value")
                end
            end
            if length(rtdb_data) > 10
                println("  ... (共$(length(rtdb_data))条记录)")
            end
        else
            println("  无数据")
        end
    catch e
        println("错误: $e")
    finally
        # 清理资源
        cleanup()
    end
end
