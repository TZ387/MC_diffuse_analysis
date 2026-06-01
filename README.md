# Theoretical Background: Ring Diffuse Power Analysis

## Overview

This document describes the method used in `ring_diffuse_power_analysis.m` to compute the fraction of incident beam power **Pd/Po** that diffusely escapes upward into air through an annular ring on the skin surface, along with its decomposition into azimuthal (phi) and polar (theta) angular bins.

The simulation is performed with MCmatlab, a Monte Carlo light transport code. The analysis script post-processes the simulation output without requiring any re-simulation.

---

## 1. Geometry and Coordinate Conventions

The MCmatlab simulation cuboid has:

- **x, y**: lateral directions, centred on the beam axis
- **z**: depth, increasing downward, with z = 0 at the top of the cuboid
- Air fills z = 0 to 0.50 cm
- Epidermis starts at z = 0.50 cm, dermis at z = 0.51 cm, subcutis at z = 0.71 cm

The ring of interest is defined in the xy-plane at the skin surface (z = 0.50 cm), with inner radius r = 0.19 cm and outer radius r = 0.21 cm, centred on the beam axis.

Angular conventions used in the output (distinct from MCmatlab's internal conventions):

- **phi**: azimuthal angle in the xy-plane, measured from the +x axis, in the range [0°, 360°). "Towards beam" is defined as phi ∈ [0°, 180°), "away from beam" as phi ∈ [180°, 360°).
- **theta_user**: polar angle measured from the skin surface (not from the surface normal). theta_user = 0° means light travelling parallel to the skin; theta_user = 90° means light travelling perpendicular to the skin (straight up). This relates to the standard polar angle from the surface normal theta_n by: **theta_user = 90° − theta_n**.

---

## 2. Why `NI_zneg` Cannot Be Used

A naive approach would be to integrate `model.MC.NI_zneg` — the normalised irradiance on the top face of the simulation cuboid (z = 0) — over the ring area. This is incorrect because the top face is 0.5 cm above the skin surface, and the air layer between them introduces lateral displacement.

A photon that exits the skin at radius r and angle theta_n from the normal travels a lateral distance

```
Δr = 0.5 cm × tan(theta_n)
```

before reaching z = 0. Near the critical angle (theta_n ≈ 48° for n_skin = 1.34), this displacement is Δr ≈ 0.56 cm — larger than the ring radius itself. `NI_zneg` therefore maps photons to the wrong radial positions, making it unsuitable for spatially resolved analysis.

The correct approach is to evaluate the upward flux directly at the skin surface.

---

## 3. Upward Flux at the Skin Surface via the Diffusion Approximation

### 3.1 The Diffusion Approximation

In a medium where scattering strongly dominates absorption (mus' >> mua), the photon fluence rate phi(r) satisfies the diffusion equation:

```
-D ∇²φ + μ_a φ = S
```

where D = 1 / (3(μ_a + μ_s')) is the diffusion coefficient and S is the source term. For considered model, the relevant parameters are:

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

## 5. Phi Decomposition

The phi angle of each ring voxel is computed from its (x, y) centre position:

```
phi = atan2(y, x)  [converted to degrees, range 0°–360°]
```

The ring is split into two halves:

- **Towards beam**: phi ∈ [0°, 180°) — the +x half-plane
- **Away from beam**: phi ∈ [180°, 360°) — the −x half-plane

Each half contributes its own J_up integral, giving Pd_towards/Po and Pd_away/Po independently. This decomposition is spatially exact — it uses the actual (x, y) position of each voxel in the ring, with no approximation.

---

## 6. Theta Angular Decomposition

### 6.1 Angular distribution of escaping light

The phi decomposition tells us how much power escapes from each half of the ring, but not at what angle. For the angular distribution we rely on the diffusion-theory result that the radiance inside a scattering medium is nearly isotropic. A Lambertian (isotropic) internal source produces an angular distribution of escaping radiance in air proportional to:

```
L_esc(θ_n) ∝ T_Fresnel(θ_n) · cos(θ_n)
```

where θ_n is the angle in air measured from the surface normal, and T_Fresnel(θ_n) is the single-pass Fresnel transmittance from skin into air at the corresponding internal angle. The cos(θ_n) factor arises from the projection of isotropic radiance onto the surface normal (Lambert's cosine law for flux).

The differential power in a solid-angle element is:

```
dP ∝ L_esc(θ_n) · sin(θ_n) dθ_n dφ
```

giving the angular weight function:

```
w(θ_n) = T_Fresnel(θ_n) · cos(θ_n) · sin(θ_n)
```

This is computed numerically on a fine grid over θ_n ∈ [0°, θ_c], and normalised so that ∫ w dθ_n = 1.

### 6.2 Critical angle cutoff

No light escapes beyond the critical angle. In the user convention:

```
theta_user = 90° − theta_n
```

The critical angle theta_n = 48.3° corresponds to theta_user = 41.7°. Therefore all bins with theta_user < 41.7° have exactly zero escaped power, as confirmed in the output.

### 6.3 Binning

The weight function w(θ_n) is integrated over each of the 30 bins (3° wide in theta_user, i.e. 3° wide in theta_n) to give the fractional power w_bin(k) in bin k. The total Pd/Po for each phi half is then distributed across bins:

```
Pd_Po_table(k, towards) = Pd_towards/Po × w_bin(k)
Pd_away_Po_table(k, away)   = Pd_away/Po    × w_bin(k)
```

This assumes the same angular shape for both phi halves, which is justified because the Fresnel-modified Lambertian emission pattern is azimuthally symmetric.

---

## 7. Assumptions and Limitations

**Diffusion approximation**: Valid when μ_s' >> μ_a and the point of interest is many scattering mean free paths from any source or boundary. At 2940 nm in skin these conditions are met at r = 0.2 cm, but would fail close to the beam axis or at much longer wavelengths where absorption dominates.

**Isotropic radiance (Lambertian emission)**: The diffusion approximation implies isotropic radiance inside the medium. In reality, close to a specular source the radiance retains some forward bias (described by higher-order P3 or adding-doubling transport theory). For the geometry here (r = 0.2 cm, approximately 10 MFPs from the beam), the bias is small.

**Single-pass Fresnel factor**: The script applies T_avg as a single Fresnel transmission correction. Multiple internal reflections (photons that are reflected at the surface, scatter back, and make another attempt) are not explicitly accounted for, but their contribution is already partially captured in the NFR distribution computed by the Monte Carlo simulation.

**Same angular shape for both phi halves**: The angular weight w_bin(k) is derived from the average Fresnel transmittance and is the same for both phi halves. In reality the two halves may have slightly different angular distributions (e.g. if the near-beam half has a slightly more forward-peaked fluence field), but this is a second-order effect.

**Epidermis optical properties used throughout**: The diffusion coefficient D is taken from the epidermis, which is the medium at the surface. The epidermis layer is only 0.01 cm (one voxel) thick, but since both iz1 and iz2 are inside the skin the gradient estimate reflects predominantly epidermal transport.

---

## 8. Key Result Interpretation

The output table shows that:

- Bins theta_user = 0°–42° contain zero power, as expected from the critical angle cutoff.
- The peak of the distribution is near theta_user ≈ 42°–48° (just above the critical angle threshold), where the Fresnel transmittance rises steeply from zero.
- Power decreases smoothly towards theta_user = 90° (normal direction) because the Lambertian cos(θ_n) weighting reduces the flux contribution at near-normal emission angles despite full Fresnel transmission there.
- The slight asymmetry between "towards" and "away" phi halves (ratio ≈ 1.14) reflects the asymmetry in the fluence rate distribution at the skin surface due to the directional nature of subsurface scattering relative to the ring position.

---

## References

- Ishimaru, A. (1978). *Wave Propagation and Scattering in Random Media*. Academic Press.
- Welch, A.J. & van Gemert, M.J.C. (1995). *Optical-Thermal Response of Laser-Irradiated Tissue*. Plenum Press.
- Haskell, R.C. et al. (1994). Boundary conditions for the diffusion equation in radiative transfer. *Journal of the Optical Society of America A*, 11(10), 2727–2741.
- Hansen, A.K. et al. (2018). Towards accurate simulation of two-photon excited luminescence from turbid biological tissue. *Journal of Biomedical Optics*, 23(12), 121622.
