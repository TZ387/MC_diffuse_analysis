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
NFR1 = double(model.MC.NFR(:,:,iz1));   % [nx, ny] first skin voxel
NFR2 = double(model.MC.NFR(:,:,iz2));   % [nx, ny] second skin voxel

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
% Angular distribution of escaping light (in air, from surface normal)
% is proportional to T_Fresnel(theta_n) * cos(theta_n) * sin(theta_n)
% where theta_n is the angle in air from the surface normal.
%
% User convention: theta_user = 90 - theta_n_deg
%   theta_user = 0  -> grazing (parallel to skin)
%   theta_user = 90 -> normal  (perpendicular to skin)
%
% Note: nothing escapes for theta_n > theta_c (total internal reflection),
% which maps to theta_user < (90 - theta_c_deg).

theta_c_deg = theta_c * 180/pi;
fprintf('\nCritical angle from normal: %.2f deg\n', theta_c_deg);
fprintf('=> Light only escapes for theta_user > %.2f deg (from surface)\n', 90 - theta_c_deg);

% Angular grid in air (from normal)
theta_n_fine = linspace(0, theta_c - 1e-9, N_fine);
theta_i_fine2 = asin((n_air/n_skin) * sin(theta_n_fine));  % angle in skin
cos_tn = cos(theta_n_fine);
cos_ti2 = cos(theta_i_fine2);
Rs2 = ((n_skin*cos_ti2 - n_air*cos_tn) ./ (n_skin*cos_ti2 + n_air*cos_tn)).^2;
Rp2 = ((n_air*cos_ti2  - n_skin*cos_tn) ./ (n_air*cos_ti2  + n_skin*cos_tn)).^2;
T_fresnel_air = 1 - 0.5*(Rs2 + Rp2);

weight_fine = T_fresnel_air .* cos_tn .* sin(theta_n_fine);
W_total_norm = trapz(theta_n_fine, weight_fine);

% Bin edges in user convention (theta_user: 0 to 90, step 3 deg)
theta_user_edges = (0:N_theta) * (90/N_theta);      % [deg]
theta_n_edges    = 90 - theta_user_edges;            % corresponding normal angles [deg]
% Bin k covers theta_user in [theta_user_edges(k), theta_user_edges(k+1)]
%             = theta_n    in [theta_n_edges(k+1),  theta_n_edges(k)    ]

w_bin = zeros(1, N_theta);
for k = 1:N_theta
    tn_lo = theta_n_edges(k+1) * pi/180;
    tn_hi = theta_n_edges(k)   * pi/180;
    if tn_lo >= theta_c
        w_bin(k) = 0;
    else
        tn_hi_eff = min(tn_hi, theta_c - 1e-9);
        idx = (theta_n_fine >= tn_lo) & (theta_n_fine <= tn_hi_eff);
        if any(idx)
            w_bin(k) = trapz(theta_n_fine(idx), weight_fine(idx)) / W_total_norm;
        end
    end
end

fprintf('Sum of angular bin weights: %.6f (should be 1.0)\n', sum(w_bin));

%% ---- BUILD OUTPUT TABLE ------------------------------------------------
% Distribute the total escaping power into angular and phi bins
Pd_Po_table = zeros(N_theta, 2);
for k = 1:N_theta
    Pd_Po_table(k, 1) = Pd_Po_towards_total * w_bin(k);
    Pd_Po_table(k, 2) = Pd_Po_away_total    * w_bin(k);
end

%% ---- DISPLAY RESULTS ---------------------------------------------------
fprintf('\n');
fprintf('==========================================================================\n');
fprintf(' Pd/Po breakdown by theta and phi\n');
fprintf(' Ring: r = %.2f to %.2f cm,  phi split at %.0f deg\n', r_inner, r_outer, PHI_SPLIT_DEG);
fprintf(' Flux estimated from NFR at skin surface via diffusion approximation\n');
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
figure('Name','Ring diffuse power analysis','Color','w','Position',[100 100 1000 420]);

theta_centres = theta_user_edges(1:end-1) + 0.5*(90/N_theta);

subplot(1,2,1);
bar(theta_centres, Pd_Po_table, 'stacked');
legend('Towards beam (\phi = 0-180°)', 'Away from beam (\phi = 180-360°)', ...
    'Location','northwest');
xlabel('\theta_{user}  [deg]   (0° = parallel,  90° = normal to skin)');
ylabel('Pd/Po per bin');
title(sprintf('Ring r = %.2f–%.2f cm: angular Pd/Po breakdown', r_inner, r_outer));
grid on;

subplot(1,2,2);
bar(theta_centres, w_bin);
xlabel('\theta_{user}  [deg]');
ylabel('Fraction of total escaped power');
title(sprintf('Angular weight function\n(Fresnel-modified Lambertian, n_{skin}=%.2f)', n_skin));
grid on;
