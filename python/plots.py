import csv
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

csv_file_SGDVP = "python/runtimes_and_costs_for_SGDVP_data.csv"
csv_file_1000G = "python/runtimes_and_costs_for_1000G_data.csv"
    
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

plt.scatter(sample_size_values_SGDVP, actual_cost_values_SGDVP, label='cost for SGDVP samples')
plt.title('Cost vs Sample Size for SGVDP Samples')
plt.xlabel('Sample Size (GB)')
plt.ylabel('Cost (USD)')
plt.show()

plt.scatter(sample_size_values_1000G, actual_cost_values_1000G, label='cost for 1000G samples')
plt.title('Cost vs Sample Size for 1000G Samples')
plt.xlabel('Sample Size (GB)')
plt.ylabel('Cost (USD)')
plt.show()

plt.scatter(dragen_runtime_values_SGDVP, actual_cost_values_SGDVP)
plt.title('DRAGEN Runtime vs Sample Size for SGDVP Samples')
plt.xlabel('DRAGEN Runtime (seconds)')
plt.ylabel('Cost (USD)')
plt.show()

plt.scatter(dragen_runtime_values_1000G, actual_cost_values_1000G)
plt.title('Cost vs DRAGEN Runtime for 1000G Samples')
plt.xlabel('DRAGEN Runtime (seconds)')
plt.ylabel('Cost (USD)')
plt.show()

plt.scatter(dragen_runtime_values_SGDVP, actual_cost_values_SGDVP)
plt.title('Cost vs DRAGEN Runtime for SGDVP Samples')
plt.xlabel('DRAGEN Runtime (seconds)')
plt.ylabel('Cost (USD)')
plt.show()

plt.scatter(sample_size_values_SGDVP, dragen_runtime_values_SGDVP)
plt.title('DRAGEN Runtime vs Sample Size for SGDVP Samples')
plt.xlabel('Sample Size (GB)')
plt.ylabel('DRAGEN Runtime (s)')
plt.show()

plt.scatter(sample_size_values_1000G, dragen_runtime_values_1000G)
plt.title('DRAGEN Runtime vs Sample Size for 1000G Samples')
plt.xlabel('Sample Size (GB)')
plt.ylabel('DRAGEN Runtime (s)')
plt.show()
