import math
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


### inputs
# total flow rate (kbd)
qt_kbd = 15
# reservoir radius (m)
re = 1500
# wellbore radius (m)
rw = 0.05
# density (kg/m3)
density = 717
# viscosity (cP)
viscosity = 0.6
# formation volume factor
fvf = 1.414
# oil gradient (bar/m)
oil_grad = 0.07

# path for csv input files
kh_path = "kh_inputs.csv"
zone_path = "zone_inputs.csv"

# kh matrix
kh_init = pd.read_csv(kh_path)
kh_pre_flip = kh_init.values
kh_pre = np.flip(kh_pre_flip, axis=0)
calc_mat = np.array([1, 9.869233e-16]).reshape(1, -1)
num_rows = kh_pre.shape[0]
calc_mat_full = np.tile(calc_mat, (num_rows, 1))
kh_mat_neg = kh_pre * calc_mat_full
kh_mat = np.absolute(kh_mat_neg)

# zone matrix
zone_init = pd.read_csv(zone_path)
zone_mat = zone_init.values
num_zones = zone_mat.shape[0]
pres_mat_kpa = zone_mat[:, 1].reshape(-1, 1)
skin_mat = zone_mat[:, 2].reshape(-1, 1)
top_mat = zone_mat[:, 3].reshape(-1, 1)
bot_mat = zone_mat[:, 4].reshape(-1, 1)

# calculations
qt = qt_kbd * 1000 / 24 / 60 / 60 / 6.29287
g = 100000 * oil_grad / density
x_num = 2 * math.pi / (viscosity / 1000 * fvf)
x_ln = math.log(re / rw)

# creating perf, k, and h matrix
h_mat = kh_mat[:, 0].reshape(-1, 1)
perf_mat = np.array(h_mat)
k_pre = kh_mat[:, 1].reshape(-1, 1)
k_pre_tile = np.tile(k_pre, (1, num_zones))
k_mat = np.hstack((h_mat, k_pre_tile))

# create skin and pressure lists
for i, (zone_top, zone_bot) in enumerate(zip(list(top_mat), list(bot_mat))):
    perf_mat[np.logical_and(perf_mat >= zone_top, perf_mat <= zone_bot)] = -15
    perf_mat[np.logical_and(perf_mat != -15, perf_mat <= zone_bot)] = -999
    k_mat[:, i + 1][np.logical_or(k_mat[:, 0] < zone_top, k_mat[:, 0] > zone_bot)] = 0

# finalize perf & k matrix
perf_mat[perf_mat != -15] = -999
k_mat = np.delete(k_mat, 0, 1)

# constant matrix (+ skin)
x_mat = x_num / (x_ln + skin_mat)

# creating kh matrix
# kh = k_av * h_delta
k_mat_m0 = np.delete(k_mat, 0, 0)
k_mat_mn = np.delete(k_mat, num_rows - 1, 0)
k_mat_avg = (k_mat_m0 + k_mat_mn) / 2
h_mat_m0 = np.delete(h_mat, 0, 0)
h_mat_mn = np.delete(h_mat, num_rows - 1, 0)
h_mat_diff = h_mat_mn - h_mat_m0
h_mat_diff_tile = np.tile(h_mat_diff, (1, num_zones))

# final kh matrices
kh_mat_avg = k_mat_avg * h_mat_diff_tile
kh_mat_sum = np.sum(kh_mat_avg, axis=0).reshape(-1, 1)

# reservoir middle matrix & pres(kPa) matrix
mid_mat = (top_mat + bot_mat) / 2
pres_mat = pres_mat_kpa * 1000

# perf mid difference matrix
first_perf = np.ones((num_zones, 1)) * mid_mat[0,0]
z_mat = mid_mat - first_perf

# calculating pwf1 terms
pwf1_n1 = kh_mat_sum * x_mat * pres_mat
pwf1_n2 = density * g * kh_mat_sum * x_mat * z_mat
pwf1_d = kh_mat_sum * x_mat

# calculating pwf1
pwf1 = (np.sum(pwf1_n1) - np.sum(pwf1_n2) - qt) / np.sum(pwf1_d)

# hydrostatics
hydro_add_mat = density * g * z_mat

# calculating pwf for each zone
pwf_mat = pwf1 + hydro_add_mat

# calculating flow for each zone
dd_mat = pres_mat - pwf_mat
q_mat = x_mat * kh_mat_sum * dd_mat

# kh percentage matrix
cum_kh_mat = np.cumsum(kh_mat_avg, axis=0)
cum_kh_max_mat = np.tile(kh_mat_sum.T, (num_rows - 1, 1))
cum_kh_perc_mat = cum_kh_mat / cum_kh_max_mat

# no dd matrix for graphing
cum_kh_mat_all = np.sum(cum_kh_mat, axis=1)
cum_kh_max_mat_all = np.sum(cum_kh_max_mat, axis=1)
cum_kh_perc_mat_graph_pre = cum_kh_mat_all / cum_kh_max_mat_all
cum_kh_perc_mat_graph = np.hstack((0, cum_kh_perc_mat_graph_pre)).reshape(-1, 1)
cum_kh_perc_mat_graph_perc = cum_kh_perc_mat_graph * 100

# flow percentage matrix
q_mat_tile = np.tile(q_mat.T, (num_rows - 1, 1))
q_perc_mat = cum_kh_perc_mat * q_mat_tile

# total flow matrix
q_total_perc_mat = np.sum(q_perc_mat, axis=1) / np.sum(q_mat)

# create matrix for graphing
q_final_mat = np.hstack((0, q_total_perc_mat)).reshape(-1, 1)
q_final_mat_perc = q_final_mat * 100

# plotting
plt.plot(q_final_mat_perc, h_mat, label='Flow')
plt.plot(perf_mat, h_mat, label='Perforations')
plt.plot(cum_kh_perc_mat_graph_perc, h_mat, '--', label='No DD Flow')
plt.ylabel("Depth (mTVD)")
plt.xlim([-20,120])
plt.xlabel("Flow Percentage (%)")
plt.legend()
plt.title("PLT", y=1.09)
plt.twiny()
plt.plot(pres_mat_kpa, mid_mat, 'ro', label = 'Reservoir Pressures')
plt.xlabel("Reservoir Pressure (kPa)")
plt.gca().invert_yaxis()
plt.legend()
plt.savefig("plt_output.png")
