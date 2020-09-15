import numpy as np
import re
import os, sys
import matplotlib.pyplot as plt
import mylib as my


"""
This is test command
"""
print("Hello there...")
"""
"""

start_str = r'Step Density TotEng Press Temp'
end_str = r'Loop time of'
pic_ext = 'eps'
"""
This is test comment.
"""
args = sys.argv[1:]
if(len(args) == 0):
	print('format:\n' + sys.argv[0] + ' model_name')
	exit(1)
model_name = args[0]
if(model_name == '.'):
	model_name = ''
	filename = 'wat.log'
else:
	filename = 'wat.log.' + model_name

# ================= parse log file ===================

with open(filename, 'r') as log_file:
	log_data = log_file.read()
mass_H = np.float(re.search(r'mass 1 ' + my.float_digital, log_data).group(0)[7:])
mass_O = np.float(re.search(r'mass 2 ' + my.float_digital, log_data).group(0)[7:])
N_molecules = np.int(re.search(r'[1-9]\d+ atoms', log_data).group(0).split()[0]) / 3
mu = (mass_O + 2 * mass_H) / 1000   # g/mol -> kg/mol

start_stats_pos = log_data.find(start_str)
end_stats_pos = log_data.find(end_str)
log_data = log_data[:end_stats_pos] + log_data[end_stats_pos + 2:]
end_stats_pos = log_data.find(end_str)
if(end_stats_pos == -1):
	end_stats_pos = len(log_data)

stats_data = log_data[start_stats_pos + len(start_str) + 2 : end_stats_pos].split('\n')
del stats_data[-1]

Npoints = len(stats_data)
Nfields = len(stats_data[0].split())

stats = np.empty([Npoints, Nfields])
for i in range(Npoints):
	stats[i, :] = np.array(stats_data[i].split()).astype(np.float)

# ===================== proc =====================
# stats: time, dens, Etot, P, T

k_B = 1.38e-23
Na = 6.02e23
t_mem = max(stats[1, 0] - stats[0, 0], 20)  # if dummp step is > 20 we can only take it
dt = 1
stab_time = -1  # from Etot(time), fs
t_tot = Npoints * dt
N_stat_points = t_tot / t_mem
s2ps = 1e12
atm2Pa = 101325
fs2s = 1e-15
Kcal2Jmol = 4184 / Na
gcm3_2_kgm3 = 1e3

stab_ind = (stats[:, 0] > stab_time)
t = stats[stab_ind, 0] * dt * fs2s                   #  fs       -> s
dens = stats[stab_ind, 1]*gcm3_2_kgm3         #  g/cm^3   -> kg/m^3
Etot = stats[stab_ind, 2]*Kcal2Jmol           #  Kcal/mol -> J
P = stats[stab_ind, 3]*atm2Pa                 #  atm      -> Pa
T = stats[stab_ind, 4]                        #  K        -> K

stat_Q = 1 / np.sqrt(N_stat_points)
dens_av = np.mean(dens)
d_dens_av = np.std(dens) * stat_Q
nu = N_molecules / Na
mass_tot = mu * nu    #  kg
V = mass_tot / dens   #  m^3
H = Etot + P*V
T_av = np.mean(T)
d_T_av = np.std(T) * stat_Q
P_av = np.mean(P)
d_P_av = np.std(P) * stat_Q
s_H = np.std(H)
Cp = s_H**2 / (k_B * T_av**2)
s_E = np.std(Etot)
Cv = s_E**2 / (k_B * T_av**2)

fig,ax = my.get_fig('$time (ps)$', r'$E_{tot} (Mcal/mol)$', r'$E_{tot}(t); C_v = ' + my.f2str(Cv / mass_tot / 1e3) + ' (KJ/K/kg)$')
ax.plot(t * s2ps, Etot / Kcal2Jmol / 1e3)
plt.savefig(model_name + 'E.' + pic_ext)

fig,ax = my.get_fig('$time (ps)$', r'$\rho (g/cm^3)$', r'$\rho(t); \rho = ' + my.f2str(dens_av / gcm3_2_kgm3) + ' \pm ' + my.f2str(d_dens_av / gcm3_2_kgm3) + '$')
ax.plot(t * s2ps, dens / gcm3_2_kgm3)
plt.savefig(model_name + 'rho.' + pic_ext)

fig,ax = my.get_fig('$time (ps)$', r'$P (atm)$', r'$P(t); P = ' + my.f2str(P_av / atm2Pa) + ' \pm ' + my.f2str(d_P_av / atm2Pa) + '$')
ax.plot(t * s2ps, P / atm2Pa)
plt.savefig(model_name + 'P.' + pic_ext)

#fig,ax = my.get_fig('$time (ps)$', r'$V (nm^3)$', r'$V(t)$')
#ax.plot(t * fs2ps, V*1e27)
#plt.savefig(model_name + 'V.' + pic_ext)

fig,ax = my.get_fig('$time (ps)$', r'$T (K)$', r'$T(t); T = ' + my.f2str(T_av) + ' \pm ' + my.f2str(d_T_av) + '$')
ax.plot(t * s2ps, T)
plt.savefig(model_name + 'T.' + pic_ext)

fig,ax = my.get_fig('$time (ps)$', r'$H (Mcal/mol)$', r'$H(t); C_p = ' + my.f2str(Cp / mass_tot / 1e3) + ' (KJ/K/kg)$')
ax.plot(t * s2ps, H / Kcal2Jmol / 1e3)
plt.savefig(model_name + 'H.' + pic_ext)

plt.show()

