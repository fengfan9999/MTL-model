function [besselbeam_intensity, fs, ux_input] = besselbeam_comsol(scaling_factor, element_no_half, f0, cycle, stress, time_array_step, excitation_time)
%=============================================================================================================
% Inputs:
% a. scaling_factor: coefficient of bessel_beam, 700 < sf < 1500 is good.
% b. element_no_half: half number of transducer, if element_no_half = 10 that
% means the transducer has 21 elements.
% c. f0: central_frequency [Hz], Max support = 4MHz or 4e6Hz.
% d. cycle: cycle number of push beam tone burst, previous one was 300.
% e. stress: initial stress from transducer [Pa].
% f. time_array_step: separations of recorded intensity data of induced
% pushing beam.
% g. excitation_time: time of pushing beam during excited.
% ** remove comment mark in the part of calculation of ARF if you want to
% ** save a shear wave animation avi.

% outpus
% a. besselbeam_intensity: push beam intensity distribution which will be used
% in comsol.
% b. ux_input: use it when saving an animation, or no need to use.
% c. fs: sampling frequency.

%=============================================================================================================
% Bessel Beam Based ARF
ROI_length = 0.04; % [m] 5cm
ROI_depth = 0.04; % [m] 5cm
PML_size = 10; % [grid]
ROI_length_size = 180; % [grid]
ROI_depth_size = 180; % [grid]
Nx = ROI_depth_size + 2 * PML_size; % [grid]
Ny = ROI_length_size + 2 * PML_size; % [grid]
dx = ROI_depth / ROI_depth_size; % [m]
dy = ROI_length / ROI_length_size; % [m]
kgrid = kWaveGrid(Nx, dx, Ny, dy);

% define medium

rho = 1000; % [kg/m^3]
c = 1540; % [m/s]
alpha = 0; % [dB/cm/MHz]

medium.density = rho * ones(Nx, Ny);
medium.sound_speed = c * ones(Nx, Ny);
medium.alpha_power = 1.01;
medium.alpha_coeff = alpha * ones(Nx, Ny);

% time array
kgrid.t_array = 0 : time_array_step : excitation_time;

% define transducer
% a. tdr position

element_size = (2 * element_no_half + 1) * dy *1000; % [mm]
fprintf('transducer diameter is %4.2f mm\n', element_size);

tdr_row = PML_size + 1; % [grid]
source.p_mask = zeros(Nx, Ny);
source.p_mask(tdr_row, Ny/2 - element_no_half : Ny/2 + element_no_half) = 1;

fs = 1/kgrid.dt; % [Hz] sampling frequency

% define bessel beam
w = 2 * pi * f0; % radial frequency
k = 2 * pi * f0 / c; % wavenumber
radius = linspace(-element_no_half * dy, element_no_half * dy, 2 * element_no_half + 1); % [m]
t = linspace(0, kgrid.dt * (kgrid.Nt - 1), kgrid.Nt);
pressure = besselj(0, scaling_factor * radius);

for nn = 1 : length(pressure)
    
        source.p(nn, :) = stress .* (pressure(nn)) .* toneBurst(fs, f0, cycle, 'Envelope', 'Gaussian'); % non-apodized
    
end

sensor.mask = zeros(Nx, Ny);
sensor.mask(tdr_row : Nx - PML_size, Ny/2 - element_no_half : Ny/2 + element_no_half) = 1;
sensor.record = {'u', 'p_max_all', 'I', 'p'};

% run simulation

input_args = {'DisplayMask', source.p_mask, 'PlotScale', 'auto', 'PMLSize', PML_size, 'MeshPlot', true, 'PMLAlpha', 10};
sensor_data_ARF = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:});

% visulization
figure;
imagesc(kgrid.y_vec * 1000, kgrid.x_vec * 1000, sensor_data_ARF.p_max_all);
colormap jet;
%h = colorbar;
axis image;
%ylabel(h, '[pa]');
xlabel('y position [mm]');
ylabel('x position [mm]');

% intensity
besselbeam_intensity = reshape(((sensor_data_ARF.p(:, kgrid.Nt - 100)).^2)./(rho * c), [ROI_depth_size, 2 * element_no_half + 1]);

% % calculation of ARF
% Nt = kgrid.Nt;
% ARF = 2 * alpha * sensor_data_ARF.Ix / c;
% 
% for nn = 1 : ROI_depth_size * (element_no_half * 2 + 1)
%     
%     [F(nn, :), lowenvelope] = envelope(ARF(nn, :), 100, 'peak');
%     
% end
% 
% % input for animation
% ux_input = (dx/(2 * rho * c)) .* F; % negative means the body force direction is pushed in objective

end
