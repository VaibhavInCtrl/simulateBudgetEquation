import numpy as np

def get_gradient(arr1, arr2):
    coefficients = np.polyfit(arr1, arr2, 1)
    return coefficients[0]

SST = 1000
QLHE = 100
def calculate_evaporation(QLHE, SST):
    Le = (2.501 - 0.00237 * SST) * 1000000
    E = QLHE / (rho * Le)
    return E


tau_x = 0.1  # Wind stress in x-direction (N/m^2)
tau_y = 0.1  # Wind stress in y-direction (N/m^2)
f = 1e-4  # Coriolis parameter (1/s), typical for mid-latitudes
h = 50  # Mixed layer depth (m)
rho = 1025  # Density of seawater (kg/m^3)
additional_currents = (0.05, 0.05)  # Other currents (u_other, v_other) in m/s
def calculate_salinity_advection(dSdx, dSdy):
    # Calculate Ekman velocities from wind stress
    u_e = tau_y / (rho * f * h)
    v_e = -tau_x / (rho * f * h)
    # Include additional current components if available
    u_total = u_e + (additional_currents[0] if additional_currents else 0)
    v_total = v_e + (additional_currents[1] if additional_currents else 0)
    # Calculate horizontal advection
    horizontal_advection = u_total * dSdx + v_total * dSdy
    return horizontal_advection


h = 50  # Mixed layer depth (m)
P = 0.0005  # Precipitation rate (m/s)
S = 35  # Surface salinity (psu)
t_array = [10,13,17]
h_array = [50,55,60]
Kz = 1e-5  # Placeholder for vertical diffusivity calculation
dhdt = get_gradient(h_array, t_array)
def calculate_salinity_residual(S, dSdx, dSdy, wh, Sh, dSdt, dSdz):
    # Net surface freshwater flux
    E = calculate_evaporation(QLHE, SST)
    net_freshwater_flux = (E - P) * S / h
    # Advective flux
    advective_flux = calculate_salinity_advection(dSdx, dSdy)
    H = 1 if wh > 0 else 0
    # Entrainment flux
    entrainment_flux = (H * (wh + dhdt) * ((S- Sh) / h))
    # Diffusive flux
    diffusive_flux = Kz / h * dSdz
    # Total tendency of salinity
    total_tendency = net_freshwater_flux - advective_flux - entrainment_flux - diffusive_flux
    # Residual calculation
    residual = dSdt - total_tendency
    return residual


Sh = 34.5  # Salinity at the base of the mixed layer (psu)
wh = 0.01  # Vertical velocity at the base of the mixed layer (m/s)
S_array = [35,37,40]
z_array = [2,3,5]
x_array = [3,5,8]
y_array = [5,6,7]
dSdt = get_gradient(S_array, t_array)
dSdx = get_gradient(S_array, x_array)
dSdy = get_gradient(S_array, y_array)
dSdz = get_gradient(S_array, z_array)
# Calculating the salinity tendency and residual
residual = calculate_salinity_residual(S, dSdx, dSdy, wh, Sh, dSdt, dSdz)

print(f"Residual: {residual:.6f} psu/s")