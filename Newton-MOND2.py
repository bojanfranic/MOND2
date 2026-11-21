import numpy as np
from scipy.special import i0e, i1e, k0e, k1e

# Fundamental Constants
G = 6.6743e-11          # Gravitational constant (m^3 kg^-1 s^-2)
M_sun = 1.98847e30      # Solar mass (kg)
kpc = 3.08567758e19     # Kiloparsec (m)
pc = kpc / 1000.0       # Parsec (m)
AU = 1.495978707e11     # Astronomical Unit (m)
a0 = 1.2e-10            # Cosmic vacuum acceleration floor (m/s^2)
H_DISK_UNIVERSAL = 0.45 * kpc  # Single universal disk thickness (450 pc)

catalog = {
    'Solar_System_Earth': {
        'is_disk': False, 'M_total': 1.0, 'r_scale_m': 1.0 * AU,
        'R_G': 1.0 * AU, 'h': 1.0 * AU, 'r_au': 1.0
    },
    'Milky_Way': {
        'is_disk': True, 'M_bulge': 0.0, 'r_bulge': 0.8,
        'M_disk': 5.60e10, 'R_d': 2.6,          # Realistic MW baryonic disk mass (McMillan 2017)
        'R_G': 19.5 * kpc, 'h': 0.45 * kpc,
        'R_opt': 8.3, 'V_opt': 235.0,
        'radii_kpc': np.array([1.0, 2.5, 5.0, 6.0, 8.2, 15.0, 19.5, 25.0]),
        'observed_v': np.array([115.0, 200.0, 233.0, 236.0, 233.0, 212.0, 200.0, 185.0])
    },
    'Andromeda_M31': {
        'is_disk': True, 'M_bulge': 1.20e10, 'r_bulge': 1.0,
        'M_disk': 7.80e10, 'R_d': 5.3,          # Realistic M31 baryonic mass (9.0e10 M_sun total)
        'R_G': 35.0 * kpc, 'h': 0.50 * kpc,
        'R_opt': 16.0, 'V_opt': 250.0,
        'radii_kpc': np.array([1.0, 2.5, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0]),
        'observed_v': np.array([165.0, 200.0, 235.0, 250.0, 248.0, 242.0, 236.0, 230.0, 225.0])
    },
    'M33': {
        'is_disk': True, 'M_total': 6.20e9, 'R_d': 2.6,  # Realistic M33 stars + HI gas mass
        'R_G': 25.0 * kpc, 'h': 0.45 * kpc,
        'R_opt': 8.3, 'V_opt': 118.0,
        'radii_kpc': np.array([0.6, 1.5, 3.0, 5.0, 8.0, 11.7, 15.0]),
        'observed_v': np.array([48.0, 73.0, 92.0, 102.0, 112.0, 120.0, 125.0])
    },
    'NGC_3198': {
        'is_disk': True, 'M_total': 1.15e10, 'R_d': 2.8, # Realistic NGC 3198 baryonic mass
        'R_G': 45.0 * kpc, 'h': 0.45 * kpc,
        'R_opt': 9.0, 'V_opt': 152.0,
        'radii_kpc': np.array([1.0, 2.5, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0]),
        'observed_v': np.array([75.0, 120.0, 144.0, 151.0, 150.0, 149.0, 149.0, 148.0])
    },
    'NGC_1560_Dwarf': {
        'is_disk': True, 'M_total': 9.50e8, 'R_d': 2.6,
        'R_G': 15.0 * kpc, 'h': 0.35 * kpc,
        'R_opt': 8.3, 'V_opt': 78.0,
        'radii_kpc': np.array([1.0, 2.5, 4.0, 6.0, 8.3]),
        'observed_v': np.array([31.0, 52.0, 64.0, 74.0, 78.0])
    },
    'NGC_1052_DF2_UDG': {
        'is_disk': False, 'M_total': 2.00e8, 'r_core': 2.2,
        'R_G': 4.5 * kpc, 'h': 4.5 * kpc,
        'radii_kpc': np.array([1.0, 2.2, 4.5, 7.6])
    },
    'Omega_Centauri_GC': {
        'is_disk': False, 'M_total': 3.50e6, 'r_core': 4.5,
        'R_G': 25.0 * pc, 'h': 25.0 * pc,
        'radii_pc': np.array([2.0, 5.0, 10.0, 15.0, 25.0])
    },
    'M87_Edge_Dropoff': {
        'is_disk': False, 'M_total': 8.80e11, 'R_d': 32.0,
        'R_G': 80.0 * kpc, 'h': 80.0 * kpc,
        'radii_kpc': np.array([20.0, 50.0, 80.0, 110.0, 140.0])
    }
}

def baryonic_aN_disk(r_kpc, cfg):
    """Exact Newtonian acceleration a_N for a 2D Freeman exponential disk + optional bulge."""
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

def compute_field_manual(r_m, a_N, R_G_m, h_m, r_edge_m):
    eta = np.sqrt((R_G_m / h_m) / 15.0)
    dilution = 1.0 / (1.0 + (r_m / r_edge_m) ** 2)
    a = 0.5 * a_N + np.sqrt((0.5 * a_N) ** 2 + a_N * (a0 * eta) * dilution)
    return np.sqrt(a * r_m) / 1000.0

def compute_field_derived(r_m, a_N, R_G_m, q_axis=0.0):
    """
    MOND2 - 100% Smooth (C^infinity) Field Equation (7 lines, zero extra turn-on needed):
    - 2D Freeman disk a_N naturally provides the steep 1-2.5 kpc climb.
    - Smooth flatness q and outer R_G dilution handle the plateau and Gaia DR3 drop-off.
    """
    q_disk = 1.0 / np.sqrt(1.0 + (R_G_m / H_DISK_UNIVERSAL) ** 2)
    q = q_axis + (1.0 - q_axis) * q_disk
    eta = 1.0 / np.sqrt(15.0 * q)

    dilution = 1.0 / (1.0 + (r_m / R_G_m) ** 2)
    a = 0.5 * a_N * (1.0 + np.sqrt(1.0 + (4.0 * a0 * eta * dilution) / a_N))
    return np.sqrt(a * r_m) / 1000.0

def compute_pure_urc(r_kpc, R_opt, V_opt):
    x = r_kpc / R_opt
    shape = np.sqrt((1.97 * (x ** 1.22)) / ((x ** 2 + 0.78 ** 2) ** 1.43) + 0.36 * (x ** 2) / (x ** 2 + 0.6 ** 2))
    norm = np.sqrt((1.97 * (1.0 ** 1.22)) / ((1.0 ** 2 + 0.78 ** 2) ** 1.43) + 0.36 * (1.0 ** 2) / (1.0 ** 2 + 0.6 ** 2))
    return V_opt * (shape / norm)

def enclosed_plummer_mass(r, M_total, r_core):
    return M_total * (r ** 3) / ((r ** 2 + r_core ** 2) ** 1.5)

# 1. Solar System
earth = catalog['Solar_System_Earth']
r_e = earth['r_au'] * AU
aN_e = G * (earth['M_total'] * M_sun) / (r_e ** 2)
v_man = compute_field_manual(r_e, aN_e, earth['R_G'], earth['h'], earth['R_G'])
v_uni = compute_field_derived(r_e, aN_e, earth['R_G'], q_axis=0.0)
print("=== 1. SOLAR SYSTEM (EARTH) ===")
print(f"r = 1.0 AU | Obs = 29.78 | Manual = {v_man:.2f} | Derived = {v_uni:.2f} | Pure URC = N/A\n")

# 2-6. Spiral Galaxies (Exact 2D Freeman Disk a_N + Realistic Baryonic Masses)
for idx, key in enumerate(['Milky_Way', 'Andromeda_M31', 'M33', 'NGC_3198', 'NGC_1560_Dwarf'], start=2):
    cfg = catalog[key]
    print(f"=== {idx}. {key.upper()} (2D FREEMAN DISK) ===")
    for r_kpc, v_obs in zip(cfg['radii_kpc'], cfg['observed_v']):
        aN = baryonic_aN_disk(r_kpc, cfg)
        v_man = compute_field_manual(r_kpc * kpc, aN, cfg['R_G'], cfg['h'], cfg['R_G'])
        v_uni = compute_field_derived(r_kpc * kpc, aN, cfg['R_G'], q_axis=0.0)
        v_urc = compute_pure_urc(r_kpc, cfg['R_opt'], cfg['V_opt'])
        print(f"r = {r_kpc:4.1f} kpc | Obs = {v_obs:5.1f} | Manual = {v_man:6.2f} | Derived = {v_uni:6.2f} | Pure URC = {v_urc:6.2f}")
    print()

# 7. Spherical Systems
print("=== 7. 3D SPHERICAL SYSTEMS (PURE URC = N/A) ===")
for key, unit, scale_unit in [('NGC_1052_DF2_UDG', 'kpc', kpc), ('Omega_Centauri_GC', 'pc', pc)]:
    cfg = catalog[key]
    r_arr = cfg['radii_kpc'] if unit == 'kpc' else cfg['radii_pc']
    for r_val in r_arr:
        r_m = r_val * scale_unit
        M_enc = enclosed_plummer_mass(r_val, cfg['M_total'], cfg['r_core'])
        aN = G * (M_enc * M_sun) / (r_m ** 2)
        s_man = compute_field_manual(r_m, aN, cfg['R_G'], cfg['h'], cfg['R_G']) / np.sqrt(3)
        s_uni = compute_field_derived(r_m, aN, cfg['R_G'], q_axis=1.0) / np.sqrt(3)
        print(f"{key:17s} | r = {r_val:4.1f} {unit:3s} | Manual sigma = {s_man:5.2f} | Derived sigma = {s_uni:5.2f} km/s")

m87 = catalog['M87_Edge_Dropoff']
for r_kpc in m87['radii_kpc']:
    r_m = r_kpc * kpc
    aN = G * (m87['M_total'] * M_sun) / (r_m ** 2)
    s_man = compute_field_manual(r_m, aN, m87['R_G'], m87['h'], m87['R_G']) / np.sqrt(3)
    s_uni = compute_field_derived(r_m, aN, m87['R_G'], q_axis=1.0) / np.sqrt(3)
    print(f"M87_Edge_Dropoff  | r = {r_kpc:5.1f} kpc| Manual sigma = {s_man:6.2f} | Derived sigma = {s_uni:6.2f} km/s")
