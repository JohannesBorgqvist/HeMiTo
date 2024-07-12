
# =========================================================
# Script: generate_read_save_data
# Written by: Johannes Borgqvist
# Date: 2024-07-12
# Description: This script generates the plot presented in
# the HeMiTo manuscript written by Borgqvist and Gretarsson
# Alexandersen.
# =========================================================
# Import our beloved libraries
from load_libraries import * 
# Make the plotting fancy with LaTeX
plt.rcParams['text.usetex'] = True
# =========================================================
# =========================================================
# Define dimensionless parameters
# =========================================================
# =========================================================
# Define dimensionless parameters
c1 = 1.75
c2 = 0.7
vareps = 0.4
# =========================================================
# =========================================================
# DEFINE FUNCTIONS
# =========================================================
# =========================================================
# ----------------------------------------------------------
# Function 1: ODE system with f=1
def aB(y, t, p):
    # Unpack and re-name our parameters so that
    # they are understandable
    c1 = p[0]
    c2 = p[1]
    vareps = p[2]
    # Define the system of ODEs
    dydt = [c1-c2*y[0] - vareps*y[0]*y[1], vareps*(y[0]*y[1]-y[1])]
    return dydt
# ----------------------------------------------------------   
# Function 2: ODE system with Hill-type conversion function
def aB_Hill(y, t, p):
    # Unpack and re-name our parameters so that
    # they are understandable
    c1 = p[0]
    c2 = p[1]
    vareps = p[2]
    # Define our non-linear Hill function
    f = 0.57+((y[1]**4)/(1+(y[1]**5)))
    #f = 0.65+((y[1]**4)/(1+(y[1]**5)))    
    # Define the system of ODEs
    dydt = [c1-c2*y[0] - vareps*f*y[0]*y[1], vareps*(f*y[0]*y[1]-y[1])]
    return dydt
# ----------------------------------------------------------
# Function 3: SS_calculator
def SS_calculator(f,g_func,v0):
    u1_star = (1/(c1*f(v0)))
    v2_star = findIntersection(f,g_func,1.5)
    u2_star = (1/(c1*f(v2_star)))
    return u1_star, u2_star, v2_star
# Define lambda functions
f1 = lambda x : 1 
f2 = lambda x : 0.57 + ((x**4)/(1+x**5))
g = lambda x: ((c2)/(c1-vareps*x))
# Function 4: Define intersection between f and g
def findIntersection(fun1,fun2,x0):
 return fsolve(lambda x : fun1(x) - fun2(x),x0)[0]
# ----------------------------------------------------------
# Approximations in the Mi-phase
# u_Mi
def u_Mi(tau,C1,C2,C3,C4):
 return (C1+C2*tau)*np.exp(-C3*tau)+C4
# v_Mi 
def v_Mi(tau,C1,C2):
 return C1*np.exp(C2*tau)
# =========================================================
# =========================================================
# Calculate TSSs and plot them
# =========================================================
# =========================================================
# Create numeric vector
v = np.linspace(0,2.75,200)
# Calculate function f=1
f1_vec = np.array([f1(v_temp) for v_temp in v])
# Calculate Hill type conversion function f
f2_vec = np.array([f2(v_temp) for v_temp in v])
# Define the function g as well
g_vec = np.array([g(v_temp) for v_temp in v])
# Print the steady states yeah
print("TSS for f1")
print("(u_2^star,v_2^star)\t=\t(%1.2f,%1.2f)\n"%(1/f1(SS_calculator(f1,g,1)[2]),SS_calculator(f1,g,1)[2]))
print(SS_calculator(f1,g,1)[2])
print("TSS for f2")
print("(u_2^star,v_2^star)\t=\t(%1.2f,%1.2f)\n"%(1/f2(SS_calculator(f2,g,1)[2]),SS_calculator(f2,g,1)[2]))
print(f2(SS_calculator(f2,g,1)[2]))
#Define the first figure
f, ax = plt.subplots(1, 1, constrained_layout=True, figsize=(20, 8))
#---------------------------------------------------------------------------------
#---------------------------------------------------------------------------------
# Plot the simulated data and the underlying model
ax.plot(v,f1_vec,color=(0/256,68/256,27/256),label="$f_{1}(v)=1$",linewidth=3.0)
ax.plot(v,f2_vec,color=(35/256,139/256,69/256),label="$f_{2}(v)=0.57+\\frac{v^{4}}{1+v^{5}}$",linewidth=3.0)
ax.plot(v,g_vec,color=(153/256,216/256,201/256),label="$g(v)=\\frac{c_{2}}{c_{1}-\\varepsilon{v}}$",linewidth=3.0)
ax.plot(np.array([SS_calculator(f1,g,1)[2]*(i+2)/(i+2) for i in range(100)]),np.linspace(0,f1(SS_calculator(f1,g,1)[2]),100),"--",color=(0,0,0),label="$v_{2}^{\\star}$",linewidth=3)
ax.plot(np.array([SS_calculator(f2,g,1)[2]*(i+2)/(i+2) for i in range(100)]),np.linspace(0,f2(SS_calculator(f2,g,1)[2]),100),"--",color=(0,0,0),label="$v_{2}^{\\star}$",linewidth=3)
# Set the y-limit
ax.set_ylim([0, 1.5])
# Set a grid and define a legend
ax.grid()
ax.legend(loc='best',prop={"size":40})
# Set the x-labels and y-labels
ax.set_xlabel(xlabel="$v$",fontsize=40)
ax.set_ylabel(ylabel="$f(v)$",fontsize=40)
ax.set_title(label="Illustration of TSSs",fontsize=50)    
ax.xaxis.set_tick_params(labelsize=35)
ax.yaxis.set_tick_params(labelsize=35)
#plt.show()
# Title and saving the figure
f.savefig('../Figures/TSS_visualisation.png')
# =========================================================
# Plot solutions, yeah?
# =========================================================
# Initial conditions
y0 = [0.2, 0.05]
#y0 = [0.19, 0.06]
#y0 = [0.15, 0.10]
# Parameters yeah
p = [c1, c2, vareps]
# Define the age or time if you will
t = np.linspace(0,35,200)
#-------------------------------------------------------------------
# Simulate patient with simple conversion term
prions_1 = odeint(aB, y0, t, args=(p,))
# Extract each species for patient
u_1, v_1 = prions_1.T
# Simulate patient with simple conversion term
prions_2 = odeint(aB_Hill, y0, t, args=(p,))
# Extract each species for patient
u_2, v_2 = prions_2.T
#---------------------------------------------------------------------------------
#---------------------------------------------------------------------------------
# Plot the trajectories in the phase plane
#---------------------------------------------------------------------------------
#---------------------------------------------------------------------------------
#Define the first figure
f_prions, ax_prions = plt.subplots(1, 2, constrained_layout=True, figsize=(20, 8))
#---------------------------------------------------------------------------------
#---------------------------------------------------------------------------------
# Subplot A
#---------------------------------------------------------------------------------
# Plot the simulated data and the underlying model
ax_prions[0].plot(t,u_1,color=(0/256,68/256,27/256),label="$u(\\tau)$",linewidth=3.0)
ax_prions[0].plot(t,v_1,color=(35/256,139/256,69/256),label="$v(\\tau)$",linewidth=3.0)
# Set a grid and define a legend
ax_prions[0].grid()
ax_prions[0].legend(loc='best',prop={"size":40})
# Set the x-labels and y-labels
ax_prions[0].set_xlabel(xlabel="Time, $\\tau$",fontsize=40)
ax_prions[0].set_ylabel(ylabel="Prion conc.",fontsize=40)
ax_prions[0].set_title(label="$f(v)=1$",fontsize=50)    
ax_prions[0].xaxis.set_tick_params(labelsize=35)
ax_prions[0].yaxis.set_tick_params(labelsize=35)
#---------------------------------------------------------------------------------
# Subplot B
#---------------------------------------------------------------------------------    
# Plot the simulated data and the underlying model
ax_prions[1].plot(t,u_2,color=(2/256,56/256,88/256),label="$u(\\tau)$",linewidth=3.0)
ax_prions[1].plot(t,v_2,color=(5/256,112/256,176/256),label="$v(\\tau)$",linewidth=3.0)
# Set a grid and define a legend
ax_prions[1].grid()
ax_prions[1].legend(loc='best',prop={"size":40})
# Set the x-labels and y-labels
ax_prions[1].set_xlabel(xlabel="Time, $\\tau$",fontsize=40)
ax_prions[1].set_ylabel(ylabel="Prion conc.",fontsize=40)
ax_prions[1].set_title(label="$f(v)=0.57+\\frac{v^{4}}{1+v^{5}}$",fontsize=50)    
ax_prions[1].xaxis.set_tick_params(labelsize=35)
ax_prions[1].yaxis.set_tick_params(labelsize=35)    
#plt.show()
f_prions.savefig('../Figures/HeMiTo_dynamics.png')
#---------------------------------------------------------------------------------
#---------------------------------------------------------------------------------
# Calculate approximations in the He-phase
#---------------------------------------------------------------------------------
#---------------------------------------------------------------------------------
# Extract an end_index
He_end_index = 14
He_end_index_2 = 25
# Define the initial conditions
u0 = y0[0]
v0 = y0[1]
# Define our approximation
v_He = np.array([v0*((i+1)/(i+1)) for i in range(He_end_index)])
u_He = np.array([(c1/c2)-((c1/c2)-u0)*np.exp(-c2*t[i]) for i in range(He_end_index)])
v_He_2 = np.array([v0*((i+1)/(i+1)) for i in range(He_end_index_2)])
u_He_2 = np.array([(c1/c2)-((c1/c2)-u0)*np.exp(-c2*t[i]) for i in range(He_end_index_2)])
#---------------------------------------------------------------------------------
#---------------------------------------------------------------------------------
# Plot approximations in the He phase
#---------------------------------------------------------------------------------
#---------------------------------------------------------------------------------
#Define the first figure
f_He, ax_He = plt.subplots(1, 2, constrained_layout=True, figsize=(20, 8))
#---------------------------------------------------------------------------------
#---------------------------------------------------------------------------------
# Subplot A
#---------------------------------------------------------------------------------
# Plot the simulated data and the underlying model
ax_He[0].plot(t[0:He_end_index],u_1[0:He_end_index],color=(0/256,68/256,27/256),label="$u(\\tau)$",linewidth=3.0)
ax_He[0].plot(t[0:He_end_index],v_1[0:He_end_index],color=(35/256,139/256,69/256),label="$v(\\tau)$",linewidth=3.0)
ax_He[0].plot(t[0:He_end_index],u_He,color=(0/256,68/256,27/256),linestyle = "dashed",label="$u_{\\mathrm{He}}(\\tau)$",linewidth=3.0)
ax_He[0].plot(t[0:He_end_index],v_He,color=(35/256,139/256,69/256),linestyle = "dashed",label="$v_{\\mathrm{He}}(\\tau)$",linewidth=3.0)
# Set the x- and y-limits
ax_He[0].set_xlim([0, round(t[He_end_index])])
# Set a grid and define a legend
ax_He[0].grid()
ax_He[0].legend(loc='best',prop={"size":40})
# Set the x-labels and y-labels
ax_He[0].set_xlabel(xlabel="Time, $\\tau$",fontsize=40)
ax_He[0].set_ylabel(ylabel="Prion conc.",fontsize=40)
ax_He[0].set_title(label="$f(v)=1$",fontsize=50)    
ax_He[0].xaxis.set_tick_params(labelsize=35)
ax_He[0].yaxis.set_tick_params(labelsize=35)
#---------------------------------------------------------------------------------
# Subplot B
#---------------------------------------------------------------------------------    
# Plot the simulated data and the underlying model
ax_He[1].plot(t[0:He_end_index_2],u_2[0:He_end_index_2],color=(2/256,56/256,88/256),label="$u(\\tau)$",linewidth=3.0)
ax_He[1].plot(t[0:He_end_index_2],v_2[0:He_end_index_2],color=(5/256,112/256,176/256),label="$v(\\tau)$",linewidth=3.0)
ax_He[1].plot(t[0:He_end_index_2],v_He_2,color=(2/256,56/256,88/256),linestyle = "dashed",label="$u_{\\mathrm{He}}(\\tau)$",linewidth=3.0)
ax_He[1].plot(t[0:He_end_index_2],u_He_2,color=(5/256,112/256,176/256),linestyle = "dashed",label="$v_{\\mathrm{He}}(\\tau)$",linewidth=3.0)
# Set the y-limit
ax_He[1].set_xlim([0, round(t[He_end_index_2])])
# Set a grid and define a legend
ax_He[1].grid()
ax_He[1].legend(loc='best',prop={"size":40})
# Set the x-labels and y-labels
ax_He[1].set_xlabel(xlabel="Time, $\\tau$",fontsize=40)
ax_He[1].set_ylabel(ylabel="Prion conc.",fontsize=40)
ax_He[1].set_title(label="$f(v)=0.57+\\frac{v^{4}}{1+v^{5}}$",fontsize=50)    
ax_He[1].xaxis.set_tick_params(labelsize=35)
ax_He[1].yaxis.set_tick_params(labelsize=35)    
# Save figure
f_He.savefig('../Figures/He_phase.png')
#---------------------------------------------------------------------------------
#---------------------------------------------------------------------------------
# Calculate approximations in the Mi-phase
#---------------------------------------------------------------------------------
#---------------------------------------------------------------------------------
# End index
Mi_end_index = He_end_index + 25
# Time series for the first curve?
t_data_original = t[He_end_index:round(1.15*Mi_end_index)]
t_data = np.array([t[i]-t[He_end_index] for i in range(He_end_index,round(1.15*Mi_end_index))])
t_data_original_2 = t[He_end_index:round(1.15*Mi_end_index)]
t_data_2 = np.array([t[i]-t[He_end_index] for i in range(He_end_index,round(1.15*Mi_end_index))])
u_data_1 = u_1[He_end_index:round(1.15*Mi_end_index)]
v_data_1 = v_1[He_end_index:round(1.15*Mi_end_index)]
u_data_2 = u_2[He_end_index:round(1.15*Mi_end_index)]
v_data_2 = v_2[He_end_index:round(1.15*Mi_end_index)]
# Fit parameters for u curve
popt_u_1, pcov = curve_fit(u_Mi, t_data, u_data_1)
# Extract parameters yeah?
C3_1, C4_1, C5_1, C6_1 = popt_u_1
u_max_1 = ((C4_1)/(C5_1))*np.exp(-(1-((C3_1*C5_1)/(C4_1))))+C6_1
# Fit parameters for u curve again
popt_u_2, pcov = curve_fit(u_Mi, t_data_2, u_data_2)
# Extract parameters
C3_2, C4_2, C5_2, C6_2 = popt_u_2
u_max_2 = ((C4_2)/(C5_2))*np.exp(-(1-((C3_2*C5_2)/(C4_2))))+C6_2
print("popt_u_1")
print(popt_u_1)
print("popt_u_2")
print(popt_u_2)
print("u_max_1\t=\t%0.3f"%(np.max([u_max_1,C6_1])))
print("u_max_2\t=\t%0.3f"%(np.max([u_max_2,C6_2])))
# Fit parameters for v curve f(v)=1
popt_v_1, pcov = curve_fit(v_Mi, t_data, v_data_1)
# Fit parameters for v curve f(v)=1
popt_v_2, pcov = curve_fit(v_Mi, t_data_2, v_data_2)
print("popt_v_1")
print(popt_v_1)
print("popt_v_2")
print(popt_v_2)
#---------------------------------------------------------------------------------
#---------------------------------------------------------------------------------
# Plot approximations in the He phase
#---------------------------------------------------------------------------------
#---------------------------------------------------------------------------------
#Define the first figure
f_Mi, ax_Mi = plt.subplots(1, 2, constrained_layout=True, figsize=(20, 8))
#---------------------------------------------------------------------------------
#---------------------------------------------------------------------------------
# Subplot A
#---------------------------------------------------------------------------------
# Plot the simulated data and the underlying model
ax_Mi[0].plot(t_data_original,u_data_1,color=(0/256,68/256,27/256),label="$u(\\tau)$",linewidth=3.0)
ax_Mi[0].plot(t_data_original,v_data_1,color=(35/256,139/256,69/256),label="$v(\\tau)$",linewidth=3.0)
ax_Mi[0].plot(t_data_original,u_Mi(t_data,*popt_u_1),color=(0/256,68/256,27/256),linestyle = "dashed",label="$u_{\\mathrm{Mi}}(\\tau)$",linewidth=3.0)
ax_Mi[0].plot(t_data_original,v_Mi(t_data,*popt_v_1),color=(35/256,139/256,69/256),linestyle = "dashed",label="$v_{\\mathrm{Mi}}(\\tau)$",linewidth=3.0)
# Set a grid and define a legend
ax_Mi[0].grid()
ax_Mi[0].legend(loc='best',prop={"size":40})
# Set the x-labels and y-labels
ax_Mi[0].set_xlabel(xlabel="Time, $\\tau$",fontsize=40)
ax_Mi[0].set_ylabel(ylabel="Prion conc.",fontsize=40)
ax_Mi[0].set_title(label="$f(v)=1$",fontsize=50)    
ax_Mi[0].xaxis.set_tick_params(labelsize=35)
ax_Mi[0].yaxis.set_tick_params(labelsize=35)
#---------------------------------------------------------------------------------
# Subplot B
#---------------------------------------------------------------------------------    
# Plot the simulated data and the underlying model
ax_Mi[1].plot(t_data_original_2,u_data_2,color=(2/256,56/256,88/256),label="$u(\\tau)$",linewidth=3.0)
ax_Mi[1].plot(t_data_original_2,v_data_2,color=(5/256,112/256,176/256),label="$v(\\tau)$",linewidth=3.0)
ax_Mi[1].plot(t_data_original_2,u_Mi(t_data_2,*popt_u_2),linestyle = "dashed",label="$u_{\\mathrm{Mi}}(\\tau)$",linewidth=3.0)
ax_Mi[1].plot(t_data_original_2,v_Mi(t_data_2,*popt_v_2),color=(5/256,112/256,176/256),linestyle = "dashed",label="$v_{\\mathrm{Mi}}(\\tau)$",linewidth=3.0)
# Set a grid and define a legend
ax_Mi[1].grid()
ax_Mi[1].legend(loc='best',prop={"size":40})
# Set the x-labels and y-labels
ax_Mi[1].set_xlabel(xlabel="Time, $\\tau$",fontsize=40)
ax_Mi[1].set_ylabel(ylabel="Prion conc.",fontsize=40)
ax_Mi[1].set_title(label="$f(v)=0.57+\\frac{v^{4}}{1+v^{5}}$",fontsize=50)    
ax_Mi[1].xaxis.set_tick_params(labelsize=35)
ax_Mi[1].yaxis.set_tick_params(labelsize=35)    
plt.show()
# Save figure
f_Mi.savefig('../Figures/Mi_phase.png')



