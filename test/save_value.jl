using JLD2

# 方法1：使用 @load 宏（推荐）
@load "rtdb_data.jld2"  # 这会将文件中的变量加载到当前作用域

# 此时应该有一个名为 rtdb_data 的变量可用
# 检查它的类型
println(typeof(rtdb_data))

# 打印数组的长度
println("rtdb_data 的长度: ", length(rtdb_data))

# 查看前几个元素
if !isempty(rtdb_data)
    println("前10个元素:")
    for (i, value) in enumerate(rtdb_data)
        if i <= 10
            println("  $i: $value")
        else
            break
        end
    end
    
    # 如果数组长度超过10，显示省略信息
    if length(rtdb_data) > 10
        println("  ... (共$(length(rtdb_data))个元素)")
    end
else
    println("数组为空")
end

# 如果需要，可以进行一些基本的统计分析
if !isempty(rtdb_data)
    println("\n基本统计信息:")
    println("  最小值: ", minimum(rtdb_data))
    println("  最大值: ", maximum(rtdb_data))
    println("  平均值: ", sum(rtdb_data) / length(rtdb_data))
    
    # 检查唯一值的数量
    unique_count = length(unique(rtdb_data))
    println("  唯一值数量: ", unique_count)
    println("  重复率: ", 1.0 - unique_count/length(rtdb_data))
end
