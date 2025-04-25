"""

This function is designed to test the unbalanced power flow analysis on unsymmetrical load nodes.
It includes the following steps:
"""
# push!(LOAD_PATH, pwd()*"/src/")

using PowerFlow

# Load the test data
# xlsx_file = "data/bus_2_test.xlsx"
# xlsx_file = "data/test_trafo_YNyn.xlsx"
# xlsx_file = "data/test_trafo_Dyn.xlsx"
xlsx_file = "data/test_trafo_Yzn.xlsx"
net = PowerFlow.import_distribution_system_data(xlsx_file)
# jpc = test_2bus_test()

# Set the options for the power flow analysis
opt = PowerFlow.options()
opt["PF"]["NR_ALG"] = "bicgstab";
opt["PF"]["ENFORCE_Q_LIMS"] = 0;
opt["PF"]["baseMVA"] = 100.0;
opt["PF"]["net_hz"] = 50.0;
# # Run the unbalanced power flow analysis
net = PowerFlow.runupf(net, opt)

# # Check the results
# PowerFlow.check_it(net)