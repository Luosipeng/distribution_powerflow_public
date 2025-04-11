"""
    Main function for the AC power flow
"""

push!(LOAD_PATH, pwd()*"/src/")
include(pwd()*"/data/case118.jl")
include(pwd()*"/data/case33bw.jl")
include(pwd()*"/data/case_reliability.jl")

using PowerFlow
# using MATLAB
using Base.Threads
using Plots

opt = PowerFlow.options() # The initial settings 
opt["PF"]["NR_ALG"] = "bicgstab";
opt["PF"]["ENFORCE_Q_LIMS"] = 0;
opt["PF"]["DC_PREPROCESS"] = 1;
#test find_islands and delete island
# mat"addpath('C:/Users/DELL/Desktop/matpower8.0/data')"
# mpc = mat"case1888rte"
# file_path = "C:/Users/13733/Desktop/onelinedata.xlsx"
# DCfile_path = "C:/Users/13733/Desktop/etap-main/parameters.xlsx"
file_path = joinpath(pwd(), "data", "acparameters.xlsx")
# file_path = "C:/Users/13733/Desktop/testdata.xlsx"

mpc, dict_bus, node_mapping, pv_curves = PowerFlow.excel2jpc(file_path)

# result_path = "C:/Users/13733/Desktop/etap-main/dcresult.xlsx"
# PowerFlow.dc_element_validate("DCLUMPLOAD","RatedKW",Dict_busdc,result_path,opt,mpc,connectedname=2)

# mpc=PowerFlow.runpf(mpc,opt)
mpc_list, isolated = PowerFlow.extract_islands(mpc)
n_islands = length(mpc_list)
println("共提取出 $(n_islands) 个孤岛")

# 创建结果数组
results_array = Vector{Any}(undef, n_islands)

println("开始多线程计算...")
t_start = time()

# 使用多线程计算每个孤岛的潮流
@threads for i in 1:n_islands
    results_array[i] = PowerFlow.runpf(mpc_list[i], opt)
end

t_end = time()
elapsed = t_end - t_start

# 构造类似@timed返回的结果
results = (value=results_array, time=elapsed)

println("计算完成，耗时: $(results.time) 秒")
PowerFlow.process_result(results, isolated, "powerflow_report.txt")

# """
#     Compare the results with ETAP
# """

result_path = joinpath(pwd(), "data", "result.xlsx")

PowerFlow.resultcompare(results[1][1],dict_bus,node_mapping)
misv,misa=PowerFlow.loadflow_result_ETAP(result_path,"output.xlsx")

plot(misv, 
    xlabel = "Bus",  # x轴标签
    ylabel = "Relative Error", # y轴标签
    title = "Relative Error Figure",
    legend = false)     # 不显示图例

