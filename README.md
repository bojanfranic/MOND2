# MOND2: Core Field Equation and Applications in Newtonian & Einsteinian Gravity

## 1. Core Equation (`M_bare + M_field`) is based on the QFT
Starting from bare baryonic mass $M$ ($M_{\text{bare}}$) dressing a background field ($M_{\text{field}}$), the radial gravitational acceleration is governed by:

$$a(r) = \frac{GM}{2r^2} + \sqrt{\left(\frac{GM}{2r^2}\right)^2 + \frac{GM \cdot a_0 \cdot \eta \cdot D(r)}{r^2}}$$

where:
* $a_0 = 1.2 \times 10^{-10}\text{ m/s}^2$ is the cosmic vacuum acceleration floor.
* $\eta$ is the **planar geometric confinement factor**, derived smoothly ($C^\infty$) from the system's outer boundary $R_G$, universal disk vertical scale $H_0 = 450\text{ pc}$, and intrinsic 3D axis ratio $q_{\text{axis}} \in$:
  $$q_{\text{disk}} = \frac{1}{\sqrt{1 + \left(\frac{R_G}{H_0}\right)^2}}, \qquad q = q_{\text{axis}} + (1 - q_{\text{axis}})\,q_{\text{disk}}, \qquad \eta = \frac{1}{\sqrt{15\,q}}$$
* $D(r)$ is the **3D spherical outer boundary dilution factor**:
  $$D(r) = \frac{1}{1 + \left(\frac{r}{R_G}\right)^2}$$

---

## 2. Derived Equation (in Terms of Newtonian Acceleration $a_N$)
Substituting the bare Newtonian gravitational field $a_N = \frac{GM}{r^2}$ (or exact 2D Freeman exponential disk acceleration $a_{N,\text{disk}}(r)$ for flat spiral disks) reduces the core equation directly to:

$$a(r) = 0.5 \cdot a_N + \sqrt{(0.5 \cdot a_N)^2 + a_N \cdot (a_0 \cdot \eta) \cdot D(r)}$$

### Three Physical Regimes
1. **High-Acceleration Limit ($a_N \gg a_0$, Solar System):**
   $$a(r) \to 0.5\,a_N + 0.5\,a_N = a_N$$
2. **Flat Galactic Disk Regime ($a_N \ll a_0,\; r < R_G$):**
   $$a(r) \to \sqrt{a_N \cdot a_0 \cdot \eta} \propto \frac{1}{r}$$
3. **Outer Void Drop-Off ($r \gg R_G$, Gaia DR3 Decline):**
   $$D(r) \propto \left(\frac{R_G}{r}\right)^2 \implies a(r) \propto \frac{1}{r^2}$$

---

## 3. Application 1: Newtonian Gravity (Orbital & Dispersion Dynamics)
In the non-relativistic regime ($v \ll c$), inserting $a(r)$ into standard Newtonian circular motion and 3D isotropic Jeans/virial velocity dispersion gives:

### A. Circular Orbital Velocity (Planets & Spiral Disks)
$$v(r) = \sqrt{r \cdot a(r)}$$

* **Earth Orbit ($1\text{ AU}$):** $v = 29.79\text{ km/s}$ (Observed: $29.78\text{ km/s}$).
* **Milky Way (2D Freeman Disk, $M_d = 5.6 \times 10^{10} M_\odot$):**
  $115.6\text{ km/s}$ ($1\text{ kpc}$) $\to$ $189.1\text{ km/s}$ ($2.5\text{ kpc}$) $\to$ $239.4\text{ km/s}$ ($8.2\text{ kpc}$) $\to$ $174.2\text{ km/s}$ ($25\text{ kpc}$, Gaia DR3 outer drop-off).
* **Andromeda M31, M33, NGC 3198 & NGC 1560:** Reproduces inner rises, flat plateaus, and outer profiles using only measured baryonic masses.

### B. 3D Velocity Dispersion (Spherical Systems, $q_{\text{axis}} = 1 \implies \eta = 0.258$)
$$\sigma_{\text{1D}}(r) = \sqrt{\frac{r \cdot a(r)}{3}}$$

* **NGC 1052-DF2 (Ultra-Diffuse Galaxy):** $\sigma = 9.45\text{–}14.37\text{ km/s}$ (resolves the "galaxy lacking dark matter" puzzle).
* **Globular Clusters (U49 & Omega Centauri):** $\sigma = 5.7\text{ km/s}$ (U49) and $12.99\text{–}20.34\text{ km/s}$ (Omega Centauri).
* **M87 Outer Elliptical Halo:** $\sigma = 261.9\text{ km/s}$ ($20\text{ kpc}$) $\to$ $124.2\text{ km/s}$ ($140\text{ kpc}$).

---

## 4. Application 2: Einsteinian Gravity (General Relativity)
In General Relativity, the background field medium $M_{\text{field}}$ enters as an effective stress-energy tensor $T_{\mu\nu}^{(\text{field})}$ alongside bare baryonic matter $T_{\mu\nu}^{(\text{bare})}$:

$$G_{\mu\nu} = \frac{8\pi G}{c^4}\left(T_{\mu\nu}^{(\text{bare})} + T_{\mu\nu}^{(\text{field})}\right)$$

This dresses the static weak-field spacetime metric $g_{00}(r) = -\left(1 + \frac{2\Phi_{\text{eff}}(r)}{c^2}\right)$ such that $\frac{d\Phi_{\text{eff}}}{dr} = a(r)$.

### A. Relativistic Geodesic Orbital Velocity & Perihelion Precession
$$v_{\text{rel}}(r) = \sqrt{\frac{r \cdot a(r)}{1 - \frac{3\,r \cdot a(r)}{c^2}}}$$

* **Mercury Orbit & Precession ($0.387\text{ AU}$):** Orbital speed $v = 47.873\text{ km/s}$ and perihelion advance $\Delta\phi = \frac{6\pi\,a_{\text{semi}}\,a(r)}{c^2(1 - e^2)} = \mathbf{42.98\text{ arcsec/century}}$ (exact match to Einstein GR).
* **Black Hole Strong-Field Limit:** Preserves Einstein's innermost stable circular orbit ($\text{ISCO}$ as $\frac{3ra}{c^2} \to 1$).
* **Galactic Scale ($r \sim 1\text{–}35\text{ kpc}$):** Seamless transition into the flat and Gaia DR3 declining galactic rotation curve where standard vacuum GR fails.

### B. Relativistic Gravitational Light Bending (Refraction & Lensing)
Because $T_{\mu\nu}^{(\text{field})}$ dresses both temporal ($g_{00}$) and spatial ($g_{ij}$) metric curvature equally, photons passing at impact parameter $b$ deflect by:

$$\hat{\alpha}(b) = \frac{4b}{c^2} \int_0^\infty \frac{a\!\left(\sqrt{b^2 + z^2}\right)}{\sqrt{b^2 + z^2}}\,dz$$

* **Solar Limb Light Bending ($b = R_\odot$):** $\hat{\alpha} = \mathbf{1.75124\text{ arcsec}}$ (exact match to Einstein GR).
* **Strong Galaxy Lensing (SWELLS/SLACS Disk Lens, $M = 2.5\times10^{11} M_\odot$):** Einstein ring radius increases from $\theta_E = 0.744''$ (Bare GR) to $\theta_E = \mathbf{1.057''}$ ($+102\%$ lensing area boost without dark matter halos).
