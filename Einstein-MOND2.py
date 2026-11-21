import numpy as np
from scipy.special import i0e, i1e, k0e, k1e
from scipy.integrate import quad
from scipy.optimize import brentq

# ==============================================================================
# 1. FUNDAMENTAL CONSTANTS
# ==============================================================================
G = 6.6743e-11                 # Gravitational constant (m^3 kg^-1 s^-2)
c = 2.99792458e8               # Speed of light (m/s)
M_sun = 1.98847e30             # Solar mass (kg)
R_sun = 6.957e8                # Solar radius (m)
kpc = 3.08567758e19            # Kiloparsec (m)
Mpc = 1000.0 * kpc             # Megaparsec (m)
pc = kpc / 1000.0              # Parsec (m)
AU = 1.495978707e11            # Astronomical Unit (m)
a0 = 1.2e-10                   # Cosmic vacuum acceleration floor (m/s^2)
H_DISK_UNIVERSAL = 0.45 * kpc  # Universal vertical confinement scale (450 pc)
RAD_TO_ARCSEC = 206264.806     # Radians to arcseconds

# ==============================================================================
# 2. ASTROPHYSICAL CATALOG (10 SYSTEMS ACROSS 10 ORDERS OF MAGNITUDE)
# ==============================================================================
catalog = {
    'Solar_System_Earth': {
        'M_total': 1.0, 'R_G': 1.0 * AU, 'r_au': 1.0, 'obs_v': 29.78
    },
    'Solar_System_Mercury': {
        'M_total': 1.0, 'R_G': 0.387098 * AU, 'r_au': 0.387098,
        'eccentricity': 0.205630, 'obs_v': 47.870, 'obs_precession': 42.98
    },
    'Milky_Way': {
        'M_bulge': 0.0, 'r_bulge': 0.8,
        'M_disk': 5.60e10, 'R_d': 2.6, 'R_G': 19.5 * kpc,
        'R_opt': 8.3, 'V_opt': 235.0,
        'radii_kpc': np.array([1.0, 2.5, 5.0, 6.0, 8.2, 15.0, 19.5, 25.0]),
        'observed_v': np.array([115.0, 200.0, 233.0, 236.0, 233.0, 212.0, 200.0, 185.0])
    },
    'Andromeda_M31': {
        'M_bulge': 1.20e10, 'r_bulge': 1.0,
        'M_disk': 7.80e10, 'R_d': 5.3, 'R_G': 35.0 * kpc,
        'R_opt': 16.0, 'V_opt': 250.0,
        'radii_kpc': np.array([1.0, 2.5, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0]),
        'observed_v': np.array([165.0, 200.0, 235.0, 250.0, 248.0, 242.0, 236.0, 230.0, 225.0])
    },
    'M33': {
        'M_total': 6.20e9, 'R_d': 2.6, 'R_G': 25.0 * kpc,
        'R_opt': 8.3, 'V_opt': 118.0,
        'radii_kpc': np.array([0.6, 1.5, 3.0, 5.0, 8.0, 11.7, 15.0]),
        'observed_v': np.array([48.0, 73.0, 92.0, 102.0, 112.0, 120.0, 125.0])
    },
    'NGC_3198': {
        'M_total': 1.15e10, 'R_d': 2.8, 'R_G': 45.0 * kpc,
        'R_opt': 9.0, 'V_opt': 152.0,
        'radii_kpc': np.array([1.0, 2.5, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0]),
        'observed_v': np.array([75.0, 120.0, 144.0, 151.0, 150.0, 149.0, 149.0, 148.0])
    },
    'NGC_1560_Dwarf': {
        'M_total': 9.50e8, 'R_d': 2.6, 'R_G': 15.0 * kpc,
        'R_opt': 8.3, 'V_opt': 78.0,
        'radii_kpc': np.array([1.0, 2.5, 4.0, 6.0, 8.3]),
        'observed_v': np.array([31.0, 52.0, 64.0, 74.0, 78.0])
    },
    'NGC_1052_DF2_UDG': {
        'M_total': 2.00e8, 'r_core': 2.2, 'R_G': 4.5 * kpc,
        'radii_kpc': np.array([1.0, 2.2, 4.5, 7.6])
    },
    'Omega_Centauri_GC': {
        'M_total': 3.50e6, 'r_core': 4.5, 'R_G': 25.0 * pc,
        'radii_pc': np.array([2.0, 5.0, 10.0, 15.0, 25.0])
    },
    'M87_Edge_Dropoff': {
        'M_total': 8.80e11, 'R_G': 80.0 * kpc,
        'radii_kpc': np.array([20.0, 50.0, 80.0, 110.0, 140.0])
    }
}

# ==============================================================================
# 3. EINSTEIN-MOND2 RELATIVISTIC FIELD & GEODESIC EQUATIONS (C^INFINITY)
# ==============================================================================
def compute_einstein_mond2(r_m, a_N, R_G_m, q_axis=0.0):
    """
    Dressed Einstein Spacetime Acceleration & Relativistic Geodesic Velocity:
    - a(r) = 0.5*a_N + sqrt[ (0.5*a_N)^2 + a_N * (a_0 * eta) * D(r) ]
    - v_rel = sqrt( (r * a) / (1 - 3 * r * a / c^2) )
    """
    q_disk = 1.0 / np.sqrt(1.0 + (R_G_m / H_DISK_UNIVERSAL) ** 2)
    q = q_axis + (1.0 - q_axis) * q_disk
    eta = 1.0 / np.sqrt(15.0 * q)
    dilution = 1.0 / (1.0 + (r_m / R_G_m) ** 2)

    a = 0.5 * a_N + np.sqrt((0.5 * a_N) ** 2 + a_N * (a0 * eta) * dilution)
    rel_factor = 1.0 - (3.0 * r_m * a) / (c ** 2)
    return np.sqrt((r_m * a) / rel_factor) / 1000.0, a

def baryonic_aN_disk(r_kpc, cfg):
    """Exact bare acceleration a_N for a 2D Freeman exponential disk + optional bulge."""
    r_m = r_kpc * kpc
    a_bulge = 0.0
    if 'M_bulge' in cfg and cfg['M_bulge'] > 0:
        r_b = cfg['r_bulge']
        m_b = cfg['M_bulge'] * (r_kpc ** 2) / ((r_kpc + r_b) ** 2) * M_sun
        a_bulge = G * m_b / (r_m ** 2)

    M_d = cfg.get('M_disk', cfg.get('M_total', 0.0)) * M_sun
    R_d_m = cfg['R_d'] * kpc
    y = r_m / (2.0 * R_d_m)
    a_disk = (G * M_d / (R_d_m ** 2)) * y * (i0e(y) * k0e(y) - i1e(y) * k1e(y))
    return a_bulge + a_disk

def enclosed_plummer_mass(r, M_total, r_core):
    return M_total * (r ** 3) / ((r ** 2 + r_core ** 2) ** 1.5)

def deflection_angle_arcsec(b_m, M_msun, R_G_m, r_core_m=0.0, q_axis=0.0):
    """Relativistic photon deflection angle alpha_hat(b) in arcseconds."""
    def integrand(u):
        z = b_m * u
        r = np.sqrt(b_m ** 2 + z ** 2 + r_core_m ** 2)
        a_N = G * (M_msun * M_sun) / (r ** 2)
        _, acc = compute_einstein_mond2(r, a_N, R_G_m, q_axis=q_axis)
        return (acc * (b_m ** 2)) / r

    integral, _ = quad(integrand, 0.0, np.inf, epsabs=1e-12, epsrel=1e-8)
    return ((4.0 / (c ** 2)) * integral) * RAD_TO_ARCSEC

def solve_einstein_radius(M_msun, R_G_m, D_L_m, D_S_m, D_LS_m, r_core_m=3.8*kpc, q_axis=0.0):
    """Solves strong-lensing Einstein ring radius theta_E in arcseconds."""
    def lens_eq(theta_arcsec):
        b_m = D_L_m * (theta_arcsec / RAD_TO_ARCSEC)
        alpha = deflection_angle_arcsec(b_m, M_msun, R_G_m, r_core_m=r_core_m, q_axis=q_axis)
        return theta_arcsec - (D_LS_m / D_S_m) * alpha

    return brentq(lens_eq, 0.05, 10.0) if lens_eq(0.05) * lens_eq(10.0) < 0 else 0.0

def compute_pure_urc(r_kpc, R_opt, V_opt):
    x = r_kpc / R_opt
    shape = np.sqrt((1.97 * (x ** 1.22)) / ((x ** 2 + 0.78 ** 2) ** 1.43) + 0.36 * (x ** 2) / (x ** 2 + 0.6 ** 2))
    norm = np.sqrt((1.97 * (1.0 ** 1.22)) / ((1.0 ** 2 + 0.78 ** 2) ** 1.43) + 0.36 * (1.0 ** 2) / (1.0 ** 2 + 0.6 ** 2))
    return V_opt * (shape / norm)

def enclosed_plummer_mass(r, M_total, r_core):
    return M_total * (r ** 3) / ((r ** 2 + r_core ** 2) ** 1.5)

# ==============================================================================
# 4. EINSTEIN-MOND2 MULTI-SCALE BENCHMARK EXECUTION
# ==============================================================================
if __name__ == "__main__":
    print("=== 1. SOLAR SYSTEM (ORBITS, PRECESSION & LIGHT REFRACTION) ===")
    for name in ['Solar_System_Earth', 'Solar_System_Mercury']:
        p = catalog[name]
        r_p = p['r_au'] * AU
        aN_p = G * (p['M_total'] * M_sun) / (r_p ** 2)
        v_ein, a_ein = compute_einstein_mond2(r_p, aN_p, p['R_G'], q_axis=0.0)
        out = f"{name:22s} | Obs = {p['obs_v']:6.3f} | Einstein-MOND2 = {v_ein:6.3f} km/s"
        if 'eccentricity' in p:
            prec = ((6.0 * np.pi * r_p * a_ein) / ((c ** 2) * (1.0 - p['eccentricity'] ** 2))) * RAD_TO_ARCSEC * 415.205
            out += f" | Precession = {prec:.2f}\"/century (GR: {p['obs_precession']}\")"
        print(out)

    sun_m2 = deflection_angle_arcsec(R_sun, 1.0, R_sun, r_core_m=0.0, q_axis=0.0)
    print(f"Solar Limb Light Bend  | Einstein-MOND2 = {sun_m2:.5f}\" (Obs/GR: 1.75000\")\n")

    for idx, key in enumerate(['Milky_Way', 'Andromeda_M31', 'M33', 'NGC_3198', 'NGC_1560_Dwarf'], start=2):
        cfg = catalog[key]
        print(f"=== {idx}. {key.upper()} (2D FREEMAN DISK) ===")
        for r_kpc, v_obs in zip(cfg['radii_kpc'], cfg['observed_v']):
            aN = baryonic_aN_disk(r_kpc, cfg)
            v_ein, _ = compute_einstein_mond2(r_kpc * kpc, aN, cfg['R_G'], q_axis=0.0)
            v_urc = compute_pure_urc(r_kpc, cfg['R_opt'], cfg['V_opt'])
            print(f"r = {r_kpc:4.1f} kpc | Obs = {v_obs:5.1f} | Einstein-MOND2 = {v_ein:6.2f} | Pure URC = {v_urc:6.2f} km/s")
        print()

    print("=== 7. 3D SPHERICAL SYSTEMS (VELOCITY DISPERSION sigma) ===")
    for key, unit, scale_unit in [('NGC_1052_DF2_UDG', 'kpc', kpc), ('Omega_Centauri_GC', 'pc', pc)]:
        cfg = catalog[key]
        r_arr = cfg['radii_kpc'] if unit == 'kpc' else cfg['radii_pc']
        for r_val in r_arr:
            r_m = r_val * scale_unit
            aN = G * (enclosed_plummer_mass(r_val, cfg['M_total'], cfg['r_core']) * M_sun) / (r_m ** 2)
            v_ein, _ = compute_einstein_mond2(r_m, aN, cfg['R_G'], q_axis=1.0)
            print(f"{key:17s} | r = {r_val:4.1f} {unit:3s} | Einstein sigma = {v_ein/np.sqrt(3):5.2f} km/s")

    m87 = catalog['M87_Edge_Dropoff']
    for r_kpc in m87['radii_kpc']:
        r_m = r_kpc * kpc
        aN = G * (m87['M_total'] * M_sun) / (r_m ** 2)
        v_ein, _ = compute_einstein_mond2(r_m, aN, m87['R_G'], q_axis=1.0)
        print(f"M87_Edge_Dropoff  | r = {r_kpc:5.1f} kpc| Einstein sigma = {v_ein/np.sqrt(3):6.2f} km/s")

    print("\n=== 8. STRONG GRAVITATIONAL LENSING (SWELLS / SLACS DISK LENS) ===")
    theta_m2 = solve_einstein_radius(2.5e11, 30.0*kpc, 700.0*Mpc, 1400.0*Mpc, 870.0*Mpc, r_core_m=3.8*kpc, q_axis=0.0)
    print(f"Einstein-MOND2 (No DM) -> theta_E = {theta_m2:.3f}\" (Obs: ~1.1-1.4\", Bare GR: 0.744\")")
