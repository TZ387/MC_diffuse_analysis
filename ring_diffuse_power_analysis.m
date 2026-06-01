%% ring_diffuse_power_analysis.m
%
% Post-processing script for MCmatlab results.
%
% PURPOSE:
%   For a ring at radius r = 0.2 cm from the incident beam axis (width
%   0.02 cm, from r = 0.19 cm to r = 0.21 cm), compute the fraction of
%   incident power that escapes upward into air (Pd/Po), broken down by:
%     - phi bracket: "towards beam" (phi in [0,180) deg) vs
%                    "away from beam" (phi in [180,360) deg)
%     - theta bracket: 30 bins of 3 deg each, from 0-3 deg up to 87-90 deg
%       where theta = 0 deg means PARALLEL to skin surface,
%             theta = 90 deg means PERPENDICULAR (normal) to skin surface.
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
PHI_SPLIT_DEG = 180;   % phi boundary between "towards" and "away" halves [deg]
% -------------------------------------------------------------------------

%% ---- DERIVED OPTICAL QUANTITIES ----------------------------------------
musp_skin = mus_skin * (1 - g_skin);          % reduced scattering [cm^-1]
D_skin    = 1 / (3 * (mua_skin + musp_skin)); % diffusion coeff [cm]
c_vacuum  = 3e10;                             % speed of light [cm/s]
v_skin    = c_vacuum / n_skin;                % speed in skin [cm/s]

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

R_grid   = sqrt(X_grid.^2 + Y_grid.^2);                   % [cm]
PHI_grid = mod(atan2(Y_grid, X_grid) * (180/pi), 360);    % [deg] 0-360

ring_mask = (R_grid >= r_inner) & (R_grid < r_outer);

fprintf('\nRing: %d voxels, area = %.4f cm^2 (geometric: %.4f cm^2)\n', ...
    sum(ring_mask(:)), sum(ring_mask(:))*dx*dy, pi*(r_outer^2 - r_inner^2));

%% ---- TOTAL Pd/Po AT SKIN SURFACE (before Fresnel escape) ---------------
% Integrate J_up over the ring area -> power fraction heading upward inside skin
J_up_ring_total   = dx * dy * sum(J_up(ring_mask));
fprintf('\nUpward flux inside skin integrated over ring: %.6f [W/W_incident]\n', J_up_ring_total);

% Fresnel average internal reflectance at skin-air interface
% R_phi_avg: fraction of upward flux that is internally reflected
% Compute numerically over all angles up to critical angle
theta_c = asin(n_air / n_skin);
N_fine  = 2000;
ti_fine = linspace(0, theta_c - 1e-9, N_fine);   % angle in skin from normal
tt_fine = asin((n_skin/n_air) * sin(ti_fine));    % angle in air

cos_ti = cos(ti_fine);
cos_tt = cos(tt_fine);
Rs = ((n_skin*cos_ti - n_air*cos_tt) ./ (n_skin*cos_ti + n_air*cos_tt)).^2;
Rp = ((n_air*cos_ti  - n_skin*cos_tt) ./ (n_air*cos_ti  + n_skin*cos_tt)).^2;
R_fresnel = 0.5*(Rs + Rp);   % unpolarised reflectance
T_fresnel = 1 - R_fresnel;

% Weight by Lambertian emission inside medium: cos(ti)*sin(ti)
w_fine  = cos_ti .* sin(ti_fine);
R_avg   = trapz(ti_fine, R_fresnel .* w_fine) / trapz(ti_fine, w_fine);
T_avg   = 1 - R_avg;

fprintf('Average Fresnel reflectance (internal): R_avg = %.4f\n', R_avg);
fprintf('Average Fresnel transmittance:          T_avg = %.4f\n', T_avg);

% Escaping power fraction
Pd_Po_total = J_up_ring_total * T_avg;
fprintf('\n=== TOTAL Pd/Po for ring (r=%.2f to %.2f cm) ===\n', r_inner, r_outer);
fprintf('Pd/Po = %.6f  (%.4f %%)\n', Pd_Po_total, Pd_Po_total*100);

%% ---- PHI SPLIT ---------------------------------------------------------
phi_towards_mask = ring_mask & (PHI_grid <  PHI_SPLIT_DEG);
phi_away_mask    = ring_mask & (PHI_grid >= PHI_SPLIT_DEG);

Pd_Po_towards_total = dx * dy * sum(J_up(phi_towards_mask)) * T_avg;
Pd_Po_away_total    = dx * dy * sum(J_up(phi_away_mask))    * T_avg;

fprintf('\nPd/Po towards beam (phi 0-180 deg):    %.6f\n', Pd_Po_towards_total);
fprintf('Pd/Po away from beam (phi 180-360 deg): %.6f\n', Pd_Po_away_total);

%% ---- ANGULAR WEIGHT FUNCTION -------------------------------------------
% We want the angular distribution of escaping light expressed as angles
% in AIR (theta_t), not inside the skin.
%
% Snell's law: n_skin * sin(theta_i) = n_air * sin(theta_t)
%
% Key point: the full escape cone inside skin (theta_i = 0 to theta_c)
% maps to the FULL hemisphere in air (theta_t = 0 to 90 deg).
% A photon at the critical angle inside skin (theta_i = theta_c) exits
% at theta_t = 90 deg (grazing). A photon going straight up (theta_i = 0)
% exits straight up (theta_t = 0). Snell's law EXPANDS the angular range,
% so ALL bins in air have non-zero weight.
%
% The weight function in air must account for the Jacobian of the
% Snell's law angle transformation. Starting from the Lambertian source
% inside skin (isotropic radiance -> flux weight cos(theta_i)*sin(theta_i)):
%
%   w(theta_i) d(theta_i) -> w(theta_t) d(theta_t)
%
% Differentiating Snell's law:
%   n_skin * cos(theta_i) d(theta_i) = n_air * cos(theta_t) d(theta_t)
%   => d(theta_i)/d(theta_t) = (n_air * cos(theta_t)) / (n_skin * cos(theta_i))
%
% So the weight per unit theta_t in air is:
%
%   w(theta_t) = T_Fresnel(theta_t) * cos(theta_i) * sin(theta_i)
%                * (n_air * cos(theta_t)) / (n_skin * cos(theta_i))
%              = T_Fresnel(theta_t) * sin(theta_i) * (n_air/n_skin) * cos(theta_t)
%
% Using Snell: sin(theta_i) = (n_air/n_skin) * sin(theta_t), so:
%
%   w(theta_t) = T_Fresnel(theta_t) * (n_air/n_skin)^2 * sin(theta_t) * cos(theta_t)
%
% The (n_air/n_skin)^2 factor is a constant and cancels in normalisation.
% So the normalised weight is:
%
%   w(theta_t) proportional to T_Fresnel(theta_t) * sin(theta_t) * cos(theta_t)
%
% where theta_t runs from 0 to 90 deg (full hemisphere in air).
% This is the same functional form as before but now correctly covering
% the full hemisphere, with T_Fresnel evaluated at each theta_t.

theta_c_deg = theta_c * 180/pi;
fprintf('\nCritical angle inside skin: %.2f deg from normal\n', theta_c_deg);
fprintf('=> Via Snell''s law this maps to the full hemisphere in air (0-90 deg)\n');
fprintf('=> ALL theta_user bins have non-zero weight\n');

% Angular grid in AIR (theta_t = angle from normal in air), full hemisphere
theta_t_fine = linspace(0, pi/2 - 1e-9, N_fine);   % [rad]

% Corresponding angle inside skin via Snell's law
theta_i_fine2 = asin((n_air/n_skin) * sin(theta_t_fine));  % always <= theta_c

cos_tt2 = cos(theta_t_fine);
cos_ti2 = cos(theta_i_fine2);

% Fresnel transmittance, computed at the inside-skin angle theta_i
% (equivalent to evaluating at theta_t via Snell's law)
Rs2 = ((n_skin*cos_ti2 - n_air*cos_tt2) ./ (n_skin*cos_ti2 + n_air*cos_tt2)).^2;
Rp2 = ((n_air*cos_ti2  - n_skin*cos_tt2) ./ (n_air*cos_ti2  + n_skin*cos_tt2)).^2;
T_fresnel_t = 1 - 0.5*(Rs2 + Rp2);

% Two weight functions are needed:
%
% (A) POWER weight: used to compute what fraction of total escaped power
%     falls in each bin. This integrates radiance over solid angle, so it
%     includes the sin(theta_t) solid-angle factor:
%
%       w_power(theta_t) = T_Fresnel * cos(theta_t) * sin(theta_t)
%
%     Integrating over theta_t gives power per unit phi. The (n_air/n_skin)^2
%     Snell Jacobian factor cancels in normalisation and is omitted.
%
% (B) INTENSITY weight: used for Pd/(Po*dOmega), i.e. radiant intensity
%     per unit incident power. This is the radiance-like quantity and does
%     NOT include sin(theta_t), because dOmega already accounts for it:
%
%       w_intensity(theta_t) = T_Fresnel * cos(theta_t)
%
%     This peaks at theta_t = 0 (normal direction, theta_user = 90 deg)
%     and falls to zero at theta_t = 90 deg (grazing, theta_user = 0 deg),
%     consistent with a Lambertian source as seen in the reference figure.

weight_power     = T_fresnel_t .* sin(theta_t_fine) .* cos_tt2;  % for Pd/Po per bin
weight_intensity = T_fresnel_t .* cos_tt2;                        % for Pd/(Po*dOmega)

W_power_norm     = trapz(theta_t_fine, weight_power);
W_intensity_norm = trapz(theta_t_fine, weight_intensity);

% Bin edges in user convention (theta_user = 90 - theta_t_deg)
%   theta_user = 0  -> grazing (theta_t = 90 deg)
%   theta_user = 90 -> normal  (theta_t = 0 deg)
% Bin k: theta_user in [theta_user_edges(k), theta_user_edges(k+1)]
%      = theta_t    in [90-theta_user_edges(k+1), 90-theta_user_edges(k)] deg

theta_user_edges = (0:N_theta) * (90/N_theta);   % [deg]

w_bin_power     = zeros(1, N_theta);   % fractional power per bin
w_bin_intensity = zeros(1, N_theta);   % unnormalised intensity shape per bin
for k = 1:N_theta
    tt_lo = (90 - theta_user_edges(k+1)) * pi/180;
    tt_hi = (90 - theta_user_edges(k))   * pi/180;
    idx = (theta_t_fine >= tt_lo) & (theta_t_fine <= tt_hi);
    if any(idx)
        w_bin_power(k)     = trapz(theta_t_fine(idx), weight_power(idx))     / W_power_norm;
        w_bin_intensity(k) = trapz(theta_t_fine(idx), weight_intensity(idx)) / W_intensity_norm;
    end
end

% For backward compatibility keep w_bin as the power weight
w_bin = w_bin_power;

fprintf('Sum of angular bin weights (power):     %.6f (should be 1.0)\n', sum(w_bin_power));
fprintf('Sum of angular bin weights (intensity): %.6f\n', sum(w_bin_intensity));

%% ---- SOLID ANGLE PER BIN -----------------------------------------------
% Each theta_user bin covers a full annular strip of the hemisphere.
% In terms of theta_t (angle from normal in air, theta_t = 90 - theta_user):
%
%   dOmega = 2*pi * integral_{theta_t_lo}^{theta_t_hi} sin(theta_t) d(theta_t)
%           = 2*pi * [cos(theta_t_lo) - cos(theta_t_hi)]
%
% This is the EXACT expression for the solid angle of an annular strip.
% It is NOT simply 2*pi*sin(theta_t)*delta_theta, which is only the
% differential approximation valid for infinitesimally narrow bins.
%
% The bins are defined in theta_user (= 90 - theta_t), so equal steps
% in theta_user correspond to equal steps in theta_t. For bin k:
%   theta_t in [90-theta_user_edges(k+1), 90-theta_user_edges(k)]
%           = [tt_lo, tt_hi]  with tt_lo < tt_hi
%
% For a phi-split half (towards or away, each spanning pi in phi),
% the solid angle is exactly half: dOmega_half = pi*[cos(tt_lo)-cos(tt_hi)].

dOmega_full = zeros(1, N_theta);   % [sr] full 2*pi annular strip (both phi halves)
dOmega_half = zeros(1, N_theta);   % [sr] half strip (one phi bracket, pi wide)
for k = 1:N_theta
    tt_lo = (90 - theta_user_edges(k+1)) * pi/180;   % lower theta_t [rad]
    tt_hi = (90 - theta_user_edges(k))   * pi/180;   % upper theta_t [rad]
    dOmega_full(k) = 2*pi * (cos(tt_lo) - cos(tt_hi));
    dOmega_half(k) =   pi * (cos(tt_lo) - cos(tt_hi));
end

%% ---- BUILD OUTPUT TABLE ------------------------------------------------
% Distribute the total escaping power into angular and phi bins
Pd_Po_table = zeros(N_theta, 2);
for k = 1:N_theta
    Pd_Po_table(k, 1) = Pd_Po_towards_total * w_bin(k);
    Pd_Po_table(k, 2) = Pd_Po_away_total    * w_bin(k);
end

% Solid-angle-normalised radiant intensity: Pd / (Po * dOmega)  [sr^-1]
% Each phi half spans pi sr per annular strip -> use dOmega_half
% The total (both halves) spans 2*pi sr per annular strip -> use dOmega_full
Pd_Po_dOmega_towards = Pd_Po_table(:,1)' ./ dOmega_half;    % [sr^-1]
Pd_Po_dOmega_away    = Pd_Po_table(:,2)' ./ dOmega_half;    % [sr^-1]
Pd_Po_dOmega_total   = sum(Pd_Po_table,2)' ./ dOmega_full;  % [sr^-1]

%% ---- DISPLAY RESULTS ---------------------------------------------------
fprintf('\n');
fprintf('==========================================================================\n');
fprintf(' Pd/Po breakdown by theta and phi\n');
fprintf(' Ring: r = %.2f to %.2f cm,  phi split at %.0f deg\n', r_inner, r_outer, PHI_SPLIT_DEG);
fprintf(' Flux from NFR at skin surface; angles in air after Snell's law refraction\n');
fprintf('==========================================================================\n');
fprintf(' theta_user [deg]  |  Towards beam  |  Away from beam  |  Total bin\n');
fprintf('--------------------------------------------------------------------------\n');
for k = 1:N_theta
    lo  = theta_user_edges(k);
    hi  = theta_user_edges(k+1);
    tot = Pd_Po_table(k,1) + Pd_Po_table(k,2);
    fprintf('  %5.1f - %5.1f     |  %11.3e  |    %11.3e   | %11.3e\n', ...
        lo, hi, Pd_Po_table(k,1), Pd_Po_table(k,2), tot);
end
fprintf('--------------------------------------------------------------------------\n');
fprintf('  TOTAL              |  %11.3e  |    %11.3e   | %11.3e\n', ...
    Pd_Po_towards_total, Pd_Po_away_total, Pd_Po_total);
fprintf('==========================================================================\n');

%% ---- OPTIONAL: SAVE TO CSV ---------------------------------------------
% header = {'theta_lo_deg','theta_hi_deg','Pd_Po_towards','Pd_Po_away','Pd_Po_bin_total'};
% data   = [theta_user_edges(1:end-1)', theta_user_edges(2:end)', ...
%           Pd_Po_table(:,1), Pd_Po_table(:,2), sum(Pd_Po_table,2)];
% writecell([header; num2cell(data)], 'ring_diffuse_results.csv');

%% ---- PLOT --------------------------------------------------------------
figure('Name','Ring diffuse power analysis','Color','w','Position',[100 100 1500 420]);

theta_centres = theta_user_edges(1:end-1) + 0.5*(90/N_theta);

% Note on x-axis convention: following the reference figure convention,
% theta = 0 deg is normal to skin (perpendicular) and theta = 90 deg is
% parallel to skin (grazing). The x-axis is therefore reversed relative
% to the internal theta_user variable (where 90 = normal).
xlabel_str = '\theta  [deg]   (0° = normal,  90° = parallel to skin)';

subplot(1,3,1);
bar(theta_centres, Pd_Po_table, 'stacked');
legend('Towards beam (\phi = 0-180°)', 'Away from beam (\phi = 180-360°)', ...
    'Location','northeast');
xlabel(xlabel_str);
ylabel('Pd/Po per bin');
title(sprintf('Ring r = %.2f–%.2f cm: Pd/Po per bin', r_inner, r_outer));
set(gca, 'XDir', 'reverse');
xlim([0 90]);
grid on;

subplot(1,3,2);
bar(theta_centres, w_bin);
xlabel(xlabel_str);
ylabel('Fraction of total escaped power');
title(sprintf('Angular weight function\n(Fresnel + Snell, n_{skin}=%.2f)', n_skin));
set(gca, 'XDir', 'reverse');
xlim([0 90]);
grid on;

subplot(1,3,3);
% Pd/(Po*dOmega): solid-angle-normalised radiant intensity [sr^-1]
% Shows power escaping per steradian per unit incident power.
% Peaks near theta = 0 (normal direction), consistent with Lambertian
% emission which goes as cos(theta_t) where theta_t = 90 - theta_user.
% The total curve uses dOmega_full (2*pi strip); the phi-split curves use
% dOmega_half (pi strip), so all three share the same y-axis units.
bar(theta_centres, [Pd_Po_dOmega_towards; Pd_Po_dOmega_away]', 'stacked');
legend('Towards beam (\phi = 0-180°)', 'Away from beam (\phi = 180-360°)', ...
    'Location','northeast');
xlabel(xlabel_str);
ylabel('Pd / (Po \cdot d\Omega)  [sr^{-1}]');
title(sprintf('Solid-angle-normalised radiant intensity\nd\\Omega = \\pi[cos(\\theta_{t,lo})-cos(\\theta_{t,hi})]  per half'));
set(gca, 'XDir', 'reverse');
xlim([0 90]);
grid on;
