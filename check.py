import numpy as np

# Constants
rho = 1025  # Density of seawater (kg/m^3)
cp = 3985  # Specific heat capacity at constant pressure (J/(kg*K))
g = 9.81  # Gravitational acceleration (m/s^2)
alpha = 2.07e-4  # Thermal expansion coefficient (1/K)
beta = 7.86e-4  # Haline contraction coefficient (1/psu)
R = 0.58  # Reflection coefficient for shortwave radiation
zeta1 = 0.35  # Attenuation coefficient for the first layer (1/m)
zeta2 = 23  # Attenuation coefficient for the second layer (1/m)

# Parameters for the mixed layer
h = 50  # Mixed layer depth (m)
Qnet = 100  # Net surface heat flux (W/m^2)
QSWR = 300  # Shortwave radiation flux (W/m^2)
P = 0.0005  # Precipitation rate (m/s)
E = 0.0003  # Evaporation rate (m/s)
S = 35  # Surface salinity (psu)
MLD = 20 # Mixed Layer Depth

# Function to calculate penetrative shortwave radiation (Qpen)
def Qpen(QSWR, h, R, zeta1, zeta2):
    return QSWR * (R * np.exp(-h/zeta1) + (1 - R) * np.exp(-h/zeta2))

# Buoyancy flux
B0 = g * (alpha * Qnet / (rho * cp) + beta * (P - E) * S / (1 - S / 1000))

# Function to calculate the density
def density(S, T):
    return rho * (1 - alpha * (T - 10) + beta * (S - 35))

# Vertical diffusivity (m2/s)
Kz = 1e-5  # Placeholder for vertical diffusivity calculation

# Temperature and salinity gradients (assumed values)
dTdz = -0.02  # Temperature gradient (K/m)
dSdz = 0.01  # Salinity gradient (psu/m)

# Budget Equations
def temperature_residual(u, v, T, dTdx, dTdy, wh, Th, QSWR):
    # Compute Qpen
    penetration = Qpen(QSWR, h, R, zeta1, zeta2)
    # Heaviside function H(wh) - equals 1 if wh > 0 (entrainment), 0 otherwise (detrainment)
    H = 1 if wh > 0 else 0
    print(H)
    # Calculate tendencies
    heating = (Qnet - penetration) / (rho * cp * h)
    print(heating)
    advection = -(u * dTdx + v * dTdy)
    print(advection)
    entrainment = -H (wh + MLD) * ((T- Th) / h)
    diffusion = -(Kz / h) * dTdz
    
    # Temperature tendency
    dTdt = heating + advection + entrainment + diffusion
    
    # Assuming a perfect model, residual is zero if all terms are included
    residual = dTdt - heating - advection - entrainment - diffusion
    return residual

# Example usage
u, v = 0.1, 0.1  # Current velocities (m/s)
T = 20  # Surface temperature (C)
Th = 19  # Temperature at the base of the mixed layer (C)
S = 35  # Surface salinity (psu)
Sh = 34.5  # Salinity at the base of the mixed layer (psu)
dTdx, dTdy = 0.001, -0.001  # Horizontal temperature gradients (K/m)
wh = 0.01  # Vertical velocity at the base of the mixed layer (m/s)

H = 1 if wh > 0 else 0
entrainment = -H (wh + MLD) * ((T- Th) / h)
print(entrainment)