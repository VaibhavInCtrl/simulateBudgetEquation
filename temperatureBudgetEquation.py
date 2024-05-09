import numpy as np

def get_gradient(arr1, arr2):
    coefficients = np.polyfit(arr1, arr2, 1)
    return coefficients[0]

# Function to calculate penetrative shortwave radiation (Qpen)
h = 50  # Mixed layer depth (m)
QSWR = 300  # Shortwave radiation flux (W/m^2)
R = 0.58  # Reflection coefficient for shortwave radiation
zeta1 = 0.35  # Attenuation coefficient for the first layer (1/m)
zeta2 = 23  # Attenuation coefficient for the second layer (1/m)
def Qpen(QSWR, h, R, zeta1, zeta2):
    return QSWR * (R * np.exp(-h/zeta1) + (1 - R) * np.exp(-h/zeta2))


tau_x = 0.1  # Wind stress in x-direction (N/m^2)
tau_y = 0.1  # Wind stress in y-direction (N/m^2)
f = 1e-4  # Coriolis parameter (1/s), typical for mid-latitudes
h = 50  # Mixed layer depth (m)
rho = 1025  # Density of seawater (kg/m^3)
additional_currents = (0.05, 0.05)  # Other currents (u_other, v_other) in m/s
def calculate_temperature_advection(dTdx, dTdy):
    # Calculate Ekman velocities from wind stress
    u_e = tau_y / (rho * f * h)
    v_e = -tau_x / (rho * f * h)
    # Include additional current components if available
    u_total = u_e + (additional_currents[0] if additional_currents else 0)
    v_total = v_e + (additional_currents[1] if additional_currents else 0)
    # Calculate horizontal advection
    horizontal_advection = u_total * dTdx + v_total * dTdy
    return horizontal_advection


Qnet = 100  # Net surface heat flux (W/m^2)
cp = 3985  # Specific heat capacity at constant pressure (J/(kg*K))
h_array = [50,55,60]
t_array = [10,13,17]
dhdt = get_gradient(h_array, t_array) # Mixed Layer Depth change with time, entrainment velocity we
Kz = 1e-5  # Vertical diffusivity (m2/s)
T_array = [35,37,40]
z_array = [2,3,5]
dTdz = get_gradient(T_array, z_array)  # Temperature gradient (K/m)
# Budget Equations
def temperature_residual(T, dTdx, dTdy, wh, Th, QSWR, dTdt):
    # Compute Qpen
    penetration = Qpen(QSWR, h, R, zeta1, zeta2)
    # Heaviside function H(wh) - equals 1 if wh > 0 (entrainment), 0 otherwise (detrainment)
    H = 1 if wh > 0 else 0

    # Calculate tendencies
    heating = (Qnet - penetration) / (rho * cp * h)
    advection = -calculate_temperature_advection(dTdx, dTdy)
    entrainment = -(H * (wh + dhdt) * ((T- Th) / h))
    diffusion = -(Kz / h) * dTdz
    
    # Assuming a perfect model, residual is zero if all terms are included
    residual = dTdt - heating - advection - entrainment - diffusion
    return residual


T = 20  # Surface temperature (C)
Th = 19  # Temperature at the base of the mixed layer (C)
Sh = 34.5  # Salinity at the base of the mixed layer (psu)
x_array = [3,5,8]
dTdx = get_gradient(T_array, x_array)
y_array = [5,6,7]
dTdy = get_gradient(T_array, y_array)  # Horizontal temperature gradients (K/m)
wh = 0.01  # Vertical velocity at the base of the mixed layer (m/s)
dTdt = get_gradient(T_array, t_array) # temperature

residual = temperature_residual(T, dTdx, dTdy, wh, Th, QSWR, dTdt)
print(f"Temperature Budget Residual: {residual:.6f} K/s")