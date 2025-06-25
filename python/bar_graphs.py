import numpy as np
import matplotlib.pyplot as plt

batch_of_two_total_individual_runtimes = 588.75 # minutes
batch_of_two_total_wall_clock_time = 317.65 # minutes

batch_of_three_total_individual_runtimes = 472.87 # minutes
batch_of_three_total_wall_clock_time = 299.25 # minutes

batch_of_four_total_individual_runtimes = 708.58 # minutes
batch_of_four_total_wall_clock_time = 391.33 # minutes

batch_of_six_total_individual_runtimes = 994 # minutes
batch_of_six_total_wall_clock_time = 372 # minutes

batch_of_eight_total_individual_runtimes = 1317 # minutes
batch_of_eight_total_wall_clock_time = 271 # minutes

batches = ("Batch of Two", "Batch of Three", "Batch of Four", "Batch of Six", "Batch of Eight")
batch_times = {
    'Wall-Clock Time': (batch_of_two_total_wall_clock_time, batch_of_three_total_wall_clock_time, batch_of_four_total_wall_clock_time, batch_of_six_total_wall_clock_time, batch_of_eight_total_wall_clock_time),
    'Sum Total of Workflow Runtimes': (batch_of_two_total_individual_runtimes, batch_of_three_total_individual_runtimes, batch_of_four_total_individual_runtimes, batch_of_six_total_individual_runtimes, batch_of_eight_total_individual_runtimes)
}

x = np.arange(len(batches))
width_of_bars = 0.25
multiplier = 0

fig, ax = plt.subplots(layout='constrained')

for attribute, measurement in batch_times.items():
    offset = width_of_bars*multiplier
    rects = ax.bar(x + offset, measurement, width=width_of_bars, label=attribute)
    ax.bar_label(rects, padding=3)
    multiplier += 1
    
ax.set_ylabel('Time (min)')
ax.set_title('Wall-Clock Runtimes & Sum Total Workflow Runtimes')
ax.set_xticks(x + width_of_bars, batches)
ax.legend(loc='upper left', ncols=5)
ax.set_ylim(0, 1500)

plt.show()