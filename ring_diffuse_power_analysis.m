%% ring_diffuse_power_analysis.m
%
% Post-processing script for MCmatlab results.
%
% PURPOSE:
%   For a ring at radius r = 0.2 cm from the incident beam axis (width
%   0.02 cm, from r = 0.19 cm to r = 0.21 cm), compute the fraction of
%   incident power that escapes upward into air (Pd/Po), broken down by:
%     - theta bin: 30 bins of 3 deg each, from 0-3 deg up to 87-90 deg
%       where theta = 0 deg means NORMAL to skin surface (perpendicular),
%             theta = 90 deg means PARALLEL to skin surface (grazing).
%
% USAGE:
%   Run your MCmatlab model first, then call this script (or paste it into
%   your model file after the runMonteCarlo line). The model object must
%   be in the workspace as 'model'.
%
% KEY PHYSICS NOTE:
%   model.MC.NI_zneg is the irradiance at the TOP of the simulation cuboid
%   (z = 0), NOT at the skin surface (z = 0.5 cm). Because air fills the
%   top 0.5 cm, photons exiting the skin at oblique angles travel laterally
%   before reaching z = 0. At angles near the critical angle (~47 deg from
%   normal), the lateral shift can exceed 0.5 cm -- far larger than the
%   ring width. Therefore NI_zneg is NOT used here.
%
%   Instead, we use model.MC.normalizedFluenceRate (NFR) at the first two
%   skin voxels (just below the air-skin interface). From the diffusion
%   approximation, the upward partial flux at the interface is:
%
%       J_up(x,y) = NFR(x,y,z1)/4 - (D/2) * dNFR/dz |_{z=z1}
%
%   where:
%     - NFR is already in [W/cm^2 / W_incident] so J_up is in [W/cm^2 / W_incident]
%     - D = 1 / (3*(mua + mus')) is the diffusion coefficient [cm]
%     - dNFR/dz is estimated via finite difference between voxels z1 and z2
%     - z increases downward in MCmatlab convention, so dNFR/dz < 0 at
%       the top of skin (NFR decreases going deeper for diffuse return light),
%       making -D/2 * dNFR/dz > 0, i.e. a positive upward flux (correct)
%
%   Note: J_up computed this way is the flux INSIDE the skin just below
%   the interface. Not all of this exits -- Fresnel reflection turns some
%   photons back. The actually escaping power is J_up * (1 - R_avg), where
%   R_avg is the average internal Fresnel reflectance. However, MCmatlab
%   already tracks individual photon Fresnel events at the interface, so
%   the NFR value at voxel z1 already "knows" about multiple-bounce
%   trapping. For a cleaner estimate we apply only the single-pass
%   transmission factor below.
%
%   VALIDITY: The diffusion approximation requires mus' >> mua and distance
%   from sources/boundaries >> 1/mus'. At 2940 nm in epidermis:
%     mus' ~ 55 cm^-1, mua ~ 1.68 cm^-1, mean free path ~ 0.018 cm
%   At r = 0.19-0.21 cm from beam axis, we are ~10 scattering MFPs away,
%   so the approximation is reasonable.
%
% -------------------------------------------------------------------------

%% ---- USER-ADJUSTABLE PARAMETERS ----------------------------------------
r_inner     = 0.19;    % [cm] inner radius of ring
r_outer     = 0.21;    % [cm] outer radius of ring

% Optical properties of epidermis at 2940 nm (from the model file)
mua_skin    = 1.68;              % [cm^-1] absorption coefficient
mus_skin    = 36.782*5*1.5;      % [cm^-1] scattering coefficient
g_skin      = 0.8;               % [-] scattering anisotropy
n_skin      = 1.34;              % [-] refractive index of epidermis
n_air       = 1.00;              % [-] refractive index of air

N_theta     = 30;      % number of theta bins (3 deg each -> 0 to 90 deg)
% -------------------------------------------------------------------------

%% ---- DERIVED OPTICAL QUANTITIES ----------------------------------------
musp_skin = mus_skin * (1 - g_skin);          % reduced scattering [cm^-1]
D_skin    = 1 / (3 * (mua_skin + musp_skin)); % diffusion coeff [cm]

fprintf('Epidermis at 2940 nm:\n');
fprintf('  mua = %.3f cm^-1,  mus = %.3f cm^-1,  mus'' = %.3f cm^-1\n', ...
    mua_skin, mus_skin, musp_skin);
fprintf('  D   = %.5f cm\n', D_skin);

%% ---- IDENTIFY SKIN SURFACE VOXELS IN Z --------------------------------
% Grid spacing
dz = model.G.Lz / model.G.nz;   % [cm] voxel size in z
dx = model.G.Lx / model.G.nx;
dy = model.G.Ly / model.G.ny;

% z-coordinates of voxel centres
z_vec = ((1:model.G.nz) - 0.5) * dz;   % [cm], size [1, nz]

% From the geometry function: skin starts at Z > 0.50 cm
% Voxel 51 has centre at z = 0.505 cm (first epidermis voxel)
% Voxel 52 has centre at z = 0.515 cm (second epidermis voxel)
z_skin_surface = 0.50;                           % [cm] air-skin interface depth

% Find the first two voxels inside skin
iz1 = find(z_vec > z_skin_surface, 1, 'first'); % first skin voxel index
iz2 = iz1 + 1;                                  % second skin voxel index

fprintf('\nSkin surface at z = %.3f cm\n', z_skin_surface);
fprintf('Using voxels iz1=%d (z=%.4f cm) and iz2=%d (z=%.4f cm) for flux\n', ...
    iz1, z_vec(iz1), iz2, z_vec(iz2));

%% ---- EXTRACT NFR SLICES AT SKIN SURFACE --------------------------------
% NFR is [nx, ny, nz], units [W/cm^2 / W_incident]
NFR1 = double(model.MC.normalizedFluenceRate(:,:,iz1));   % [nx, ny] first skin voxel
NFR2 = double(model.MC.normalizedFluenceRate(:,:,iz2));   % [nx, ny] second skin voxel

% Upward partial flux at the interface [W/cm^2 / W_incident]
% J_up = NFR/4 - (D/2) * dNFR/dz
% dNFR/dz = (NFR2 - NFR1)/dz  (z increases downward, so this is the downward gradient)
dNFR_dz = (NFR2 - NFR1) / dz;           % [W/cm^3 / W_incident]
J_up = NFR1/4 - (D_skin/2) * dNFR_dz;  % [W/cm^2 / W_incident]

% Clamp negative values (unphysical; can occur near beam axis where
% diffusion approximation breaks down)
J_up(J_up < 0) = 0;

fprintf('Max J_up: %.4e, Min J_up (before clamp): %.4e\n', ...
    max(J_up(:)), min(NFR1(:)/4 - (D_skin/2)*dNFR_dz(:)));

%% ---- BUILD XY GRID AND RING MASK ---------------------------------------
x_vec = ((1:model.G.nx) - 0.5) * dx - model.G.Lx/2;
y_vec = ((1:model.G.ny) - 0.5) * dy - model.G.Ly/2;

[X_grid, Y_grid] = meshgrid(x_vec, y_vec);
X_grid = X_grid.';   % [nx, ny]
Y_grid = Y_grid.';

R_grid = sqrt(X_grid.^2 + Y_grid.^2);   % [cm]

ring_mask = (R_grid >= r_inner) & (R_grid < r_outer);

fprintf('\nRing: %d voxels, area = %.4f cm^2 (geometric: %.4f cm^2)\n', ...
    sum(ring_mask(:)), sum(ring_mask(:))*dx*dy, pi*(r_outer^2 - r_inner^2));

%% ---- TOTAL Pd/Po -------------------------------------------------------
% Integrate J_up over the ring area -> power fraction heading upward inside skin
J_up_ring_total = dx * dy * sum(J_up(ring_mask));
fprintf('\nUpward flux inside skin integrated over ring: %.6f [W/W_incident]\n', J_up_ring_total);

% Fresnel average internal reflectance at skin-air interface
% Compute numerically over all angles up to critical angle
theta_c = asin(n_air / n_skin);
N_fine  = 2000;
ti_fine = linspace(0, theta_c - 1e-9, N_fine);   % angle inside skin from normal [rad]
tt_fine = asin((n_skin/n_air) * sin(ti_fine));    % angle in air [rad]

cos_ti = cos(ti_fine);
cos_tt = cos(tt_fine);
Rs = ((n_skin*cos_ti - n_air*cos_tt) ./ (n_skin*cos_ti + n_air*cos_tt)).^2;
Rp = ((n_air*cos_ti  - n_skin*cos_tt) ./ (n_air*cos_ti  + n_skin*cos_tt)).^2;
R_fresnel = 0.5*(Rs + Rp);   % unpolarised reflectance

% Weight by Lambertian emission inside medium: cos(ti)*sin(ti)
w_fine = cos_ti .* sin(ti_fine);
R_avg  = trapz(ti_fine, R_fresnel .* w_fine) / trapz(ti_fine, w_fine);
T_avg  = 1 - R_avg;

fprintf('Average Fresnel reflectance (internal): R_avg = %.4f\n', R_avg);
fprintf('Average Fresnel transmittance:          T_avg = %.4f\n', T_avg);

% Total escaping power fraction
Pd_Po_total = J_up_ring_total * T_avg;
fprintf('\n=== TOTAL Pd/Po for ring (r=%.2f to %.2f cm) ===\n', r_inner, r_outer);
fprintf('Pd/Po = %.6f  (%.4f %%)\n', Pd_Po_total, Pd_Po_total*100);

%% ---- ANGULAR WEIGHT FUNCTION -------------------------------------------
% We want the angular distribution of escaping light expressed as angles
% in AIR (theta_t), not inside the skin.
%
% Snell's law: n_skin * sin(theta_i) = n_air * sin(theta_t)
%
% Key point: the full escape cone inside skin (theta_i = 0 to theta_c)
% maps to the FULL hemisphere in air (theta_t = 0 to 90 deg).
% A photon at the critical angle inside skin exits at theta_t = 90 deg
% (grazing). A photon going straight up exits at theta_t = 0 (normal).
% Snell's law EXPANDS the angular range, so ALL bins in air have non-zero
% weight.
%
% The weight function in air accounts for the Jacobian of the Snell's law
% transformation and the Fresnel transmittance:
%
%   w_power(theta_t)     proportional to  T_Fresnel(theta_t) * sin(theta_t) * cos(theta_t)
%   w_intensity(theta_t) proportional to  T_Fresnel(theta_t) * cos(theta_t)
%
% where theta_t runs from 0 to 90 deg (full hemisphere in air).
%
% w_power     -> fraction of total escaped power in each bin (includes
%                sin(theta_t) solid-angle factor)
% w_intensity -> radiant intensity shape Pd/(Po*dOmega), does NOT include
%                sin(theta_t) since dOmega already accounts for it; peaks
%                at theta_t = 0 (normal), consistent with Lambertian source

theta_c_deg = theta_c * 180/pi;
fprintf('\nCritical angle inside skin: %.2f deg from normal\n', theta_c_deg);
fprintf('=> Via Snell''s law this maps to the full hemisphere in air (0-90 deg)\n');
fprintf('=> ALL theta bins have non-zero weight\n');

% Angular grid in AIR (theta_t = angle from normal in air), full hemisphere
theta_t_fine  = linspace(0, pi/2 - 1e-9, N_fine);   % [rad]

% Corresponding angle inside skin via Snell's law
theta_i_fine2 = asin((n_air/n_skin) * sin(theta_t_fine));  % always <= theta_c

cos_tt2 = cos(theta_t_fine);
cos_ti2 = cos(theta_i_fine2);

% Fresnel transmittance at each theta_t (computed via inside-skin angle)
Rs2 = ((n_skin*cos_ti2 - n_air*cos_tt2) ./ (n_skin*cos_ti2 + n_air*cos_tt2)).^2;
Rp2 = ((n_air*cos_ti2  - n_skin*cos_tt2) ./ (n_air*cos_ti2  + n_skin*cos_tt2)).^2;
T_fresnel_t = 1 - 0.5*(Rs2 + Rp2);

weight_power     = T_fresnel_t .* sin(theta_t_fine) .* cos_tt2;
weight_intensity = T_fresnel_t .* cos_tt2;

W_power_norm     = trapz(theta_t_fine, weight_power);
W_intensity_norm = trapz(theta_t_fine, weight_intensity);

% Bin edges in theta_t (angle from normal in air, 0 = normal, 90 = grazing)
% N_theta equal bins spanning 0 to 90 deg
theta_t_edges = (0:N_theta) * (pi/2 / N_theta);   % [rad]

w_bin_power     = zeros(1, N_theta);
w_bin_intensity = zeros(1, N_theta);
for k = 1:N_theta
    tt_lo = theta_t_edges(k);
    tt_hi = theta_t_edges(k+1);
    idx = (theta_t_fine >= tt_lo) & (theta_t_fine <= tt_hi);
    if any(idx)
        w_bin_power(k)     = trapz(theta_t_fine(idx), weight_power(idx))     / W_power_norm;
        w_bin_intensity(k) = trapz(theta_t_fine(idx), weight_intensity(idx)) / W_intensity_norm;
    end
end

fprintf('Sum of angular bin weights (power):     %.6f (should be 1.0)\n', sum(w_bin_power));
fprintf('Sum of angular bin weights (intensity): %.6f\n', sum(w_bin_intensity));

%% ---- SOLID ANGLE PER BIN -----------------------------------------------
% Each bin covers a full annular strip of the hemisphere (full 2*pi in phi).
% Exact solid angle:
%   dOmega = 2*pi * [cos(theta_t_lo) - cos(theta_t_hi)]
%
% This is exact for any bin width; the differential approximation
% 2*pi*sin(theta_t)*delta_theta is only valid for infinitesimally narrow bins.

dOmega = zeros(1, N_theta);
for k = 1:N_theta
    dOmega(k) = 2*pi * (cos(theta_t_edges(k)) - cos(theta_t_edges(k+1)));
end

%% ---- BUILD OUTPUT ARRAYS -----------------------------------------------
% Pd/Po in each theta bin
Pd_Po_bin = Pd_Po_total * w_bin_power;   % [1], size [1, N_theta]

% Solid-angle-normalised radiant intensity: Pd / (Po * dOmega)  [sr^-1]
Pd_Po_dOmega = Pd_Po_bin ./ dOmega;

%% ---- DISPLAY RESULTS ---------------------------------------------------
theta_t_edges_deg = theta_t_edges * 180/pi;

fprintf('\n');
fprintf('================================================================\n');
fprintf(' Pd/Po breakdown by theta\n');
fprintf(' Ring: r = %.2f to %.2f cm\n', r_inner, r_outer);
fprintf(' Angles in air after Snell''s law refraction\n');
fprintf(' theta = 0 deg: normal to skin;  theta = 90 deg: grazing\n');
fprintf('================================================================\n');
fprintf('  theta [deg]    |    Pd/Po per bin   | Pd/(Po*dOmega) [sr^-1]\n');
fprintf('----------------------------------------------------------------\n');
for k = 1:N_theta
    fprintf('  %5.1f - %5.1f   |   %11.3e     |   %11.3e\n', ...
        theta_t_edges_deg(k), theta_t_edges_deg(k+1), ...
        Pd_Po_bin(k), Pd_Po_dOmega(k));
end
fprintf('----------------------------------------------------------------\n');
fprintf('  TOTAL           |   %11.3e\n', Pd_Po_total);
fprintf('================================================================\n');

%% ---- OPTIONAL: SAVE TO CSV ---------------------------------------------
% header = {'theta_lo_deg','theta_hi_deg','Pd_Po_bin','Pd_Po_per_sr'};
% data   = [theta_t_edges_deg(1:end-1)', theta_t_edges_deg(2:end)', ...
%           Pd_Po_bin', Pd_Po_dOmega'];
% writecell([header; num2cell(data)], 'ring_diffuse_results.csv');

%% ---- PLOT --------------------------------------------------------------
theta_centres_deg = (theta_t_edges_deg(1:end-1) + theta_t_edges_deg(2:end)) / 2;

figure('Name','Ring diffuse power analysis','Color','w','Position',[100 100 900 420]);

subplot(1,2,1);
bar(theta_centres_deg, Pd_Po_bin, 'FaceColor', [0.2 0.5 0.8]);
xlabel('\theta [deg]  (0° = normal,  90° = grazing)');
ylabel('Pd/Po per bin');
title(sprintf('Ring r = %.2f–%.2f cm: power per bin', r_inner, r_outer));
xlim([0 90]);
grid on;

subplot(1,2,2);
% Pd/(Po*dOmega): solid-angle-normalised radiant intensity [sr^-1]
% Peaks near theta = 0 (normal direction), consistent with Lambertian
% cos(theta_t) emission. Near-grazing bins have large dOmega but low
% radiant intensity because Fresnel transmittance -> 0 at the critical angle.
bar(theta_centres_deg, Pd_Po_dOmega, 'FaceColor', [0.8 0.4 0.2]);
xlabel('\theta [deg]  (0° = normal,  90° = grazing)');
ylabel('Pd / (Po \cdot d\Omega)  [sr^{-1}]');
title(sprintf('Radiant intensity (ring r = %.2f–%.2f cm)', r_inner, r_outer));
xlim([0 90]);
grid on;
