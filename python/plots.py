import csv
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

csv_file_SGDVP = "runtimes_in_minutes_and_costs_for_SGDVP_data.csv"
csv_file_1000G = "runtimes_in_minutes_and_costs_for_1000G_data.csv"
    
df_SGDVP = pd.read_csv(csv_file_SGDVP)
df_1000G = pd.read_csv(csv_file_1000G)  
   
dragen_runtime_values_SGDVP = df_SGDVP['dragen_runtime'].values
workflow_runtime_values_SGDVP = df_SGDVP['workflow_runtime'].values
approximate_cost_values_SGDVP = df_SGDVP['approximate_cost'].values
actual_cost_values_SGDVP = df_SGDVP['actual_cost'].values
sample_size_values_SGDVP = df_SGDVP['sample_size'].values
approximate_cost_values_SGDVP = df_SGDVP['approximate_cost'].values

dragen_runtime_values_1000G = df_1000G['dragen_runtime'].values
workflow_runtime_values_1000G = df_1000G['workflow_runtime'].values
actual_cost_values_1000G = df_1000G['actual_cost'].values
sample_size_values_1000G = df_1000G['sample_size'].values

# use polyfit to get line of best fit for 1000G data
x_1000G = np.sort(dragen_runtime_values_1000G)
y_1000G = actual_cost_values_1000G
grad_1000G, c_1000G = np.polyfit(x_1000G, y_1000G, 1)
plt.scatter(x_1000G, y_1000G, c='red')
plt.title('Cost vs DRAGEN Runtime for 1000G Samples')
plt.xlabel('DRAGEN Runtime (minutes)')
plt.ylabel('Cost (USD)')
plt.plot(x_1000G, grad_1000G * x_1000G + c_1000G, '-')
plt.show()

# # use polyfit to get line of best fit for SGDVP data
# x_SGDVP = np.sort(dragen_runtime_values_SGDVP)
# y_SGDVP = actual_cost_values_SGDVP
# grad_SGDVP, c_SGDVP = np.polyfit(x_SGDVP, y_SGDVP, 1)

# # corr_matrix = np.corrcoef(x_SGDVP, y_SGDVP)
# # corr = corr_matrix[0,1]
# # R_sq = corr**2

# plt.scatter(x_SGDVP, y_SGDVP, c='red')
# plt.plot(x_SGDVP, grad_SGDVP * x_SGDVP + c_SGDVP, '-')
# plt.title('Cost vs DRAGEN Runtime for SGDVP Samples')
# plt.xlabel('DRAGEN Runtime (minutes)')
# plt.ylabel('Cost (USD)')
# plt.show()
