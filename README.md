# Theoretical Background: Ring Diffuse Power Analysis

## Overview

This document describes the method used in `ring_diffuse_power_analysis.m` to compute the fraction of incident beam power **Pd/Po** that diffusely escapes upward into air through an annular ring on the skin surface, along with its decomposition into polar angle (theta) bins.

The simulation is performed with MCmatlab, a Monte Carlo light transport code. The analysis script post-processes the simulation output without requiring any re-simulation.

---

## 1. Geometry and Coordinate Conventions

The MCmatlab simulation cuboid has:

- **x, y**: lateral directions, centred on the beam axis
- **z**: depth, increasing downward, with z = 0 at the top of the cuboid
- Air fills z = 0 to 0.50 cm
- Epidermis starts at z = 0.50 cm, dermis at z = 0.51 cm, subcutis at z = 0.71 cm

The ring of interest is defined in the xy-plane at the skin surface (z = 0.50 cm), with inner radius r = 0.19 cm and outer radius r = 0.21 cm, centred on the beam axis.

Angular convention used in the output:

- **θ**: polar angle from the surface normal. θ = 0° means light travelling perpendicular to the skin (straight up, normal direction); θ = 90° means light travelling parallel to the skin (grazing). This is identical to the standard refracted angle θ_t in air after Snell's law.

No azimuthal (phi) decomposition is performed. The geometry is rotationally symmetric about the beam axis (on-axis beam, flat homogeneous layers), so the azimuthal distribution of escaping power is uniform and carries no information. In particular, any split into "towards beam" vs "away from beam" halves would yield exactly 50/50 by the reflection symmetry of the geometry, regardless of the physics.

---

## 2. Why `NI_zneg` Cannot Be Used

A naive approach would be to integrate `model.MC.NI_zneg` — the normalised irradiance on the top face of the simulation cuboid (z = 0) — over the ring area. This is incorrect because the top face is 0.5 cm above the skin surface, and the air layer between them introduces lateral displacement.

A photon that exits the skin at radius r and angle theta_t from the normal travels a lateral distance

```
Δr = 0.5 cm × tan(theta_t)
```

before reaching z = 0. At grazing angles (theta_t → 90°), this displacement diverges — even at theta_t = 70° the shift is Δr ≈ 1.4 cm, far larger than the ring radius. `NI_zneg` therefore maps photons to the wrong radial positions, making it unsuitable for spatially resolved analysis.

The correct approach is to evaluate the upward flux directly at the skin surface.

---

## 3. Upward Flux at the Skin Surface via the Diffusion Approximation

### 3.1 The Diffusion Approximation

In a medium where scattering strongly dominates absorption (μ_s' >> μ_a), the photon fluence rate φ(r) satisfies the diffusion equation:

```
-D ∇²φ + μ_a φ = S
```

where D = 1 / (3(μ_a + μ_s')) is the diffusion coefficient and S is the source term. For the considered model, the relevant parameters are:

| Parameter | Value |
|-----------|-------|
| μ_a | 1.68 cm⁻¹ |
| μ_s (scattering) | 275.9 cm⁻¹ |
| g (anisotropy) | 0.80 |
| μ_s' = μ_s(1−g) | 55.2 cm⁻¹ |
| D = 1/(3(μ_a + μ_s')) | 0.00586 cm |
| n | 1.34 |

The condition μ_s' >> μ_a (55.2 >> 1.68) is well satisfied, and at r = 0.19–0.21 cm from the beam axis the distance from the source is approximately 10 scattering mean free paths (1/μ_s' ≈ 0.018 cm), so the diffusion approximation is reasonable in this region.

### 3.2 Partial Current Formula

In diffusion theory, the net flux vector is **J** = −D ∇φ. The upward-directed partial current (the flux of energy moving in the −z direction, i.e. towards the surface) at a planar boundary is given by the Fick partial current expression (Ishimaru 1978; Welch & van Gemert 1995):

```
J_up(x, y) = φ/4 − (D/2) · dφ/dz
```

where:

- φ is the fluence rate [W cm⁻² / W_incident], taken from `model.MC.normalizedFluenceRate` at the first skin voxel (iz1, centre at z = 0.505 cm)
- dφ/dz is estimated by finite difference between the first and second skin voxels:

```
dφ/dz ≈ (φ(iz2) − φ(iz1)) / dz
```

with dz = 0.01 cm.

The sign convention is consistent with MCmatlab's z axis (z increases downward): near the top of the skin, φ decreases going deeper (diffuse light is heading back up), so dφ/dz < 0, making the term −(D/2)·dφ/dz positive and adding to the upward flux correctly.

Any voxels where J_up evaluates to a negative value (unphysical, can occur very close to the beam axis where the diffusion approximation breaks down) are clamped to zero.

### 3.3 Why voxels iz1 and iz2, not the air voxel

The air voxel just above the interface (iz = 50, centre at z = 0.495 cm) has near-zero scattering and absorption. The fluence rate there reflects a mix of the downward-going beam and upward-going diffuse light and is not well described by the diffusion equation. Using iz1 and iz2 — both fully inside the skin — keeps the gradient estimate within the valid diffusive regime.

---

## 4. Total Pd/Po

The total power fraction escaping from the ring is obtained in two steps.

**Step 1**: integrate J_up over the ring area to get the upward flux inside the skin heading toward the interface:

```
P_up/Po = dx · dy · Σ_{ring voxels} J_up(x, y)
```

**Step 2**: not all of this flux escapes — some fraction is internally reflected at the skin-air interface. The average Fresnel reflectance for a Lambertian (isotropic) source inside the medium is:

```
R_avg = ∫₀^{θ_c} R(θ_i) cos(θ_i) sin(θ_i) dθ_i  /  ∫₀^{θ_c} cos(θ_i) sin(θ_i) dθ_i
```

where θ_i is the angle of incidence inside the skin, θ_c = arcsin(n_air/n_skin) is the critical angle (≈ 48.3° for n_skin = 1.34), and R(θ_i) is the unpolarised Fresnel reflectance:

```
R(θ_i) = ½ (Rs² + Rp²)

Rs = (n_skin cos θ_i − n_air cos θ_t) / (n_skin cos θ_i + n_air cos θ_t)
Rp = (n_air  cos θ_i − n_skin cos θ_t) / (n_air  cos θ_i + n_skin cos θ_t)
```

with θ_t given by Snell's law: n_skin sin θ_i = n_air sin θ_t. For n_skin = 1.34, this gives R_avg ≈ 0.067 and T_avg = 1 − R_avg ≈ 0.933.

The final result is:

```
Pd/Po = (dx · dy · Σ J_up) × T_avg
```

---

## 5. Theta Angular Decomposition

### 5.1 Snell's Law and the Full Hemisphere

A key physical point is that the angular distribution must be expressed as angles **in air after refraction**, not as angles inside the tissue. Snell's law maps the escape cone inside skin to the full hemisphere in air:

```
n_skin · sin(θ_i) = n_air · sin(θ_t)
```

- A photon inside skin at θ_i = 0° (straight up) exits at θ_t = 0° (straight up).
- A photon at the critical angle θ_i = θ_c ≈ 48.3° exits at θ_t = 90° (grazing).
- Photons at θ_i > θ_c are totally internally reflected and never escape.

The full escape cone inside skin (0° to 48.3°) therefore maps to the **full hemisphere in air** (0° to 90°). This means **all angular bins receive non-zero power** — including near-grazing angles — because Snell's law compresses the angular range inward from the skin side and expands it on the air side.

### 5.2 Angular weight function in air

The radiance inside a scattering medium is nearly isotropic (diffusion approximation). To obtain the angular distribution of escaping power expressed as angles in air, we must account for the Jacobian of the Snell's law transformation when changing the integration variable from θ_i (inside skin) to θ_t (in air).

Differentiating Snell's law:

```
n_skin · cos(θ_i) dθ_i = n_air · cos(θ_t) dθ_t
=> dθ_i/dθ_t = (n_air · cos(θ_t)) / (n_skin · cos(θ_i))
```

Starting from the Lambertian weight inside skin, cos(θ_i) · sin(θ_i) dθ_i, and substituting:

```
w(θ_t) = T_Fresnel(θ_t) · cos(θ_i) · sin(θ_i) · (n_air · cos(θ_t)) / (n_skin · cos(θ_i))
        = T_Fresnel(θ_t) · sin(θ_i) · (n_air / n_skin) · cos(θ_t)
```

Using Snell's law to substitute sin(θ_i) = (n_air/n_skin) · sin(θ_t):

```
w(θ_t) = T_Fresnel(θ_t) · (n_air/n_skin)² · sin(θ_t) · cos(θ_t)
```

The prefactor (n_air/n_skin)² is a constant and cancels upon normalisation, leaving the normalised weight function:

```
w(θ_t) ∝ T_Fresnel(θ_t) · sin(θ_t) · cos(θ_t)
```

where θ_t ∈ [0°, 90°]. T_Fresnel(θ_t) is close to its maximum near θ_t = 0° (normal incidence, low Fresnel reflection) and decreases toward zero as θ_t → 90°, because those photons originated near the critical angle inside skin where T_Fresnel → 0. The sin(θ_t) · cos(θ_t) factor peaks at 45°. The combined distribution therefore peaks somewhere in the mid-range and is non-zero across the full hemisphere.

### 5.3 Two weight functions

Two distinct weight functions are used in the script:

**(A) Power weight** — used to compute the fraction of total escaped power in each angular bin. This integrates radiance over solid angle and therefore includes the sin(θ_t) solid-angle factor:

```
w_power(θ_t) ∝ T_Fresnel(θ_t) · sin(θ_t) · cos(θ_t)
```

**(B) Intensity weight** — used for the Pd/(Po·dΩ) plot. This represents radiant intensity (power per steradian) and does **not** include sin(θ_t), because the solid angle dΩ already accounts for it:

```
w_intensity(θ_t) ∝ T_Fresnel(θ_t) · cos(θ_t)
```

This peaks at θ_t = 0° (normal direction) and falls to zero at θ_t = 90° (grazing), consistent with a Lambertian source.

### 5.4 Solid angle per bin

Each bin spans a full annular strip of the hemisphere (all azimuthal angles, 0 to 2π). The solid angle is computed exactly as:

```
dΩ = 2π · [cos(θ_t_lo) − cos(θ_t_hi)]
```

where θ_t_lo and θ_t_hi are the lower and upper boundaries of the bin. The **exact** cosine-difference form is used rather than the differential approximation 2π·sin(θ_t)·Δθ_t, which would introduce error for finite bin widths.

### 5.5 Binning

The power weight w_power(θ_t) is integrated numerically over each of the 30 bins (3° wide each, covering 0° to 90°) to give the fractional power w_bin(k). The per-bin Pd/Po values are then:

```
Pd_Po_bin(k) = Pd/Po × w_bin(k)
```

The solid-angle-normalised radiant intensity per bin is:

```
Pd_Po_dOmega(k) = Pd/Po × w_intensity_bin(k) / dΩ(k)     [sr⁻¹]
```

---

## 6. Assumptions and Limitations

**Diffusion approximation**: Valid when μ_s' >> μ_a and the point of interest is many scattering mean free paths from any source or boundary. At 2940 nm in skin these conditions are met at r = 0.2 cm, but would fail close to the beam axis or at much longer wavelengths where absorption dominates.

**Isotropic radiance (Lambertian emission)**: The diffusion approximation implies isotropic radiance inside the medium. In reality, close to a source the radiance retains some forward bias (described by higher-order P3 or adding-doubling transport theory). For the geometry here (r = 0.2 cm, approximately 10 MFPs from the beam), the bias is small.

**Single-pass Fresnel factor**: The script applies T_avg as a single Fresnel transmission correction. Multiple internal reflections (photons that are reflected at the surface, scatter back, and make another attempt) are not explicitly accounted for, but their contribution is already partially captured in the NFR distribution computed by the Monte Carlo simulation.

**Epidermis optical properties used throughout**: The diffusion coefficient D is taken from the epidermis, which is the medium at the surface. The epidermis layer is only 0.01 cm (one voxel) thick, but since both iz1 and iz2 are inside the skin the gradient estimate reflects predominantly epidermal transport.

---

## 7. Key Result Interpretation

The output table and plots show that:

- All angular bins from θ = 0° to 90° receive non-zero power. This is a direct consequence of Snell's law: the narrow escape cone inside skin (0° to 48.3° from normal) is refractively expanded to cover the full hemisphere in air.
- The Pd/Po per bin plot peaks in the mid-range (roughly θ = 30°–50°), reflecting the combined effect of the sin(θ_t)·cos(θ_t) solid-angle weighting and the Fresnel transmittance dropping toward zero at grazing angles.
- The Pd/(Po·dΩ) plot (radiant intensity) peaks near θ = 0° (normal direction) and decreases monotonically toward grazing, consistent with the Lambertian cos(θ_t) dependence. The difference between the two plots is purely geometric: bins near grazing (θ ≈ 90°) have large solid angles, making their Pd/Po contribution relatively large even though their radiant intensity is small.
- Near-grazing bins (θ close to 90°) have non-zero but small radiant intensity, because photons exiting at those angles originated near the critical angle inside skin where the Fresnel transmittance approaches zero.

---

## References

- Ishimaru, A. (1978). *Wave Propagation and Scattering in Random Media*. Academic Press.
- Welch, A.J. & van Gemert, M.J.C. (1995). *Optical-Thermal Response of Laser-Irradiated Tissue*. Plenum Press.
- Haskell, R.C. et al. (1994). Boundary conditions for the diffusion equation in radiative transfer. *Journal of the Optical Society of America A*, 11(10), 2727–2741.
- Hansen, A.K. et al. (2018). Towards accurate simulation of two-photon excited luminescence from turbid biological tissue. *Journal of Biomedical Optics*, 23(12), 121622.
