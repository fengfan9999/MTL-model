%% Bessel Beam Based ARF
% clear; clc; close all;
condition = 1;

ROI_length = 0.04; % [m] 0.04
ROI_depth = 0.04; % [m] 0.04
PML_size = 10; % [grid]
ROI_length_size = 180; % [grid] 180
ROI_depth_size = 260; % [grid] 260 5MHz
Nx = ROI_depth_size + 2 * PML_size; % [grid]
Ny = ROI_length_size + 2 * PML_size; % [grid]
dx = ROI_depth / ROI_depth_size; % [m]
dy = ROI_length / ROI_length_size; % [m]
kgrid = kWaveGrid(Nx, dx, Ny, dy);
ensemble_number = 45;

% define medium

rho = 1000; % [kg/m^3]
c = 1540; % [m/s]
alpha = 0.7; % [dB/cm/MHz]

medium.density = rho * ones(Nx, Ny);
medium.sound_speed = c * ones(Nx, Ny);
medium.alpha_power = 1.01; % 1.01
medium.alpha_coeff = alpha * ones(Nx, Ny);

% time array
cfl = 0.3;
t_end = 100e-6; % [sec] 100e-6
% kgrid.makeTime(medium.sound_speed);
t_array_step = 4.5e-8; % [sec]6.5e-8 [180 180]  4.5e-8 [260 180]
kgrid.t_array = 0 : t_array_step : t_end;

% define transducer
% a. tdr position 

element_no_half = 21;
element_size = (2 * element_no_half + 1) * dy *1000; % [mm]
fprintf('transducer diameter is %4.2f mm\n', element_size);

tdr_row = PML_size + 1; % [grid]
source.p_mask = zeros(Nx, Ny);
source.p_mask(tdr_row, Ny/2 - element_no_half : Ny/2 + element_no_half) = 1;

f0 = 5e6; % [Hz] 5 MHz
fs = 1/kgrid.dt; % [Hz]
%cycle = round((ROI_depth/c) * f0);

% define bessel beam
% ==================================================================================
w = 2 * pi * f0; % radial frequency
k = 2 * pi * f0 / c; % wavenumber
scaling_factor = 0; 

%radius = element_no_half * dy; % [m]
radius = linspace(-element_no_half * dy, element_no_half * dy, 2 * element_no_half + 1); % [m]

%beta = sqrt(k.^2 - scaling_factor.^2); % propagation constant
t = linspace(0, kgrid.dt * (kgrid.Nt - 1), kgrid.Nt);

%z = 0 : kgrid.dt * c : kgrid.dt * (kgrid.Nt - 1) * c; % [mm] respective spatial point at each time point
%z = 0;

%source.p_mask = zeros(Nx, Ny);
%source.p_mask(PML_size + 10, Ny/ 2 - element_no_half : Ny/2 + element_no_half) = 1;
pressure = besselj(0, scaling_factor * radius);

%cycle = round((t_end * f0));
cycle = (t_end * f0);
% cycle = 2000;
stress = 1; % [Pa]

if condition == 1
    % bessel aperture
    for nn = 1 : length(pressure)
    
        source.p(nn, :) = stress .* (pressure(nn)) .* toneBurst(fs, f0, cycle, 'Envelope', 'Gaussian'); 


    end

elseif condition == 2
    
     % focus aperture
       source.p = stress .* toneBurst(fs, f0, cycle, 'Envelope', 'Gaussian'); 
       source.p = focus(kgrid, source.p, source.p_mask,[0.0133,0], c); % focus
else
    t1 = 0 : 1/fs : (cycle/(6 * f0));
    t2 = (cycle/(6 * f0) + kgrid.dt) : 1/fs : (2 * cycle/(6 * f0));
    t3 = (2 * cycle/(6 * f0) + kgrid.dt) : 1/fs : (3 * cycle/(6 * f0));
    t4 = (3 * cycle/(6 * f0) + kgrid.dt) : 1/fs : (4 * cycle/(6 * f0));
    t5 = (4 * cycle/(6 * f0) + kgrid.dt) : 1/fs : (5 * cycle/(6 * f0));
    t6 = (5 * cycle/(6 * f0) + kgrid.dt) : 1/fs : (6 * cycle/(6 * f0)) + kgrid.dt;

    t1_ex = [t1 zeros(size(t2)) zeros(size(t3)) zeros(size(t4)) zeros(size(t5)) zeros(size(t6))];
    t2_ex = [zeros(size(t1)) t2 zeros(size(t3)) zeros(size(t4)) zeros(size(t5)) zeros(size(t6))];
    t3_ex = [zeros(size(t1)) zeros(size(t2)) t3 zeros(size(t4)) zeros(size(t5)) zeros(size(t6))];
    t4_ex = [zeros(size(t1)) zeros(size(t2)) zeros(size(t3)) t4 zeros(size(t5)) zeros(size(t6))];
    t5_ex = [zeros(size(t1)) zeros(size(t2)) zeros(size(t3)) zeros(size(t4)) t5 zeros(size(t6))];
    t6_ex = [zeros(size(t1)) zeros(size(t2)) zeros(size(t3)) zeros(size(t4)) zeros(size(t5)) t6];
    t_ex = 0 : 1/fs : cycle/f0;

    source_physics_1 = sin(2 * pi * f0 * t1_ex);
    source_physics_2 = sin(2 * pi * f0 * t2_ex);
    source_physics_3 = sin(2 * pi * f0 * t3_ex);
    source_physics_4 = sin(2 * pi * f0 * t4_ex);
    source_physics_5 = sin(2 * pi * f0 * t5_ex);
    source_physics_6 = sin(2 * pi * f0 * t6_ex);

    
    source_1 = focus(kgrid, source_physics_1, source.p_mask,[-0.0133,0], c);
    source_2 = focus(kgrid, source_physics_2, source.p_mask,[-0.0066,0], c);
    source_3 = focus(kgrid, source_physics_3, source.p_mask,[0.0001,0], c);
    source_4 = focus(kgrid, source_physics_4, source.p_mask,[0.0068,0], c);
    source_5 = focus(kgrid, source_physics_5, source.p_mask,[0.0135,0], c);
    source_6 = focus(kgrid, source_physics_6, source.p_mask,[0.02,0], c);

%     source_4 = focus(kgrid, source_physics_3, source.p_mask,[0.0133,0], c);

    source.p = source_1(:, 1 : kgrid.Nt) + source_2(:, 1 : kgrid.Nt) + source_3(:, 1 : kgrid.Nt) + source_4(:, 1 : kgrid.Nt) + source_5(:, 1 : kgrid.Nt) + source_6(:, 1 : kgrid.Nt);
end
% ===================================================================================

% define sensor

sensor.mask = zeros(Nx, Ny);
sensor.mask(tdr_row : Nx - PML_size, Ny/2 - element_no_half : Ny/2 + element_no_half) = 1;
sensor.record = {'u', 'p_max_all', 'I', 'p'};

% run simulation

input_args = {'DisplayMask', source.p_mask, 'PlotScale', 'auto', 'PMLSize', PML_size, 'MeshPlot', true, 'PMLAlpha', 10};
sensor_data_ARF = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:});

% visulization
figure;
imagesc(kgrid.y_vec * 1000, kgrid.x_vec * 1000, 20 * log10(sensor_data_ARF.p_max_all./(max(sensor_data_ARF.p_max_all, [], 'all'))));
% colormap jet;
caxis([-10 0]);
h = colorbar;
axis image;
ylabel(h, 'dB');
xlabel('mm');
ylabel('mm');

% apodization
figure;
plot(kgrid.y_vec(11 : end - 10) * 1000, sensor_data_ARF.p_max_all(11, 11 : end - 10)./(max(sensor_data_ARF.p_max_all(11, :), [], 'all')));

% central SPL dB
figure;
p_max_all = sensor_data_ARF.p_max_all;
q = 20 * log10(p_max_all(11 : 270, 100)./(max(p_max_all(11 : 270, 100), [], 'all')));
[up, low] = envelope(q, 100, 'peak');
x_vec = kgrid.x_vec;
plot(x_vec(11 : 270) * 1000 + 20, q);
hold on
plot(x_vec(11 : 270) * 1000 + 20, up);
hold off
ylim([-15, 2]);
xlabel('distance (mm)');
ylabel('pressure (dB)');
legend('signal', 'envelope');
title('compressional spl');
% figure;
% q = 20 * log10(sensor_data_ARF.p_max_all(11 : 270, 100)./(max(sensor_data_ARF.p_max_all(11 : 270, 100), [], 'all')));
% [up, low] = envelope(q, 100, 'peak');
% plot(kgrid.x_vec(11 : 270) * 1000 + 20, q);
% hold on
% plot(kgrid.x_vec(11 : 270) * 1000 + 20, up);
% hold off
% ylim([-15, 2]);
% xlabel('distance (mm)');
% ylabel('pressure (dB)');
% legend('signal', 'envelope');

% inclusion
% inclusion_radius = 12;
% row = 100;
% inclusion = makeDisc(Nx, Ny, row, PML_size + 100, inclusion_radius, true); % default PML_size + 70

%% calculation ARF

Nt = kgrid.Nt;
% F = rho * sensor_data_ARF.ux; % body force for each point
% for nn = 1 : ROI_depth_size * (element_no_half * 2 + 1)
%     
%     F(nn, :) = abs(hilbert(F(nn, :))); % envelop of ARF
%     
% end

ARF = 2 * alpha * sensor_data_ARF.Ix / c;
for nn = 1 : ROI_depth_size * (element_no_half * 2 + 1)
    
    [F(nn, :), lowenvelope] = envelope(ARF(nn, :), 100, 'peak');
    
end
% input for step 2

ux_input = (dx/(2 * rho * c)) .* F; % negative means the body force direction is pushed in objective

%% Shear wave
clear medium; clear source; clear sensor;

% define far field/rayleigh distance

wavelength = c/f0;
a_tdr = element_no_half * dy;
z_max = (0.339 * (2 * a_tdr).^2)/wavelength;

% define medium

c_shear = 4; % [m/s]
c_shear_inclusion = 5; % [m/s]

inclusion_radius = 0;
row = 80;
inclusion = makeDisc(Nx, Ny, row, PML_size + 100, inclusion_radius, true); % default PML_size + 70

medium.density = rho * ones(Nx, Ny);
medium.sound_speed_compression = c * ones(Nx, Ny);
medium.sound_speed_shear = c_shear * ones(Nx, Ny);
medium.sound_speed_shear(inclusion == 1) = c_shear_inclusion;

% create time array

%cfl = 0.3;
%t_end = 2e-3; % [sec] 1.2e-2 for 3/4.5      1.7e-3
kgrid.t_array = 0 : t_array_step : 2e-3;

% source
source.u_mask = zeros(Nx, Ny);
source.u_mask(tdr_row : Nx - PML_size, Ny/2 - element_no_half : Ny/2 + element_no_half) = 1;
source.ux = ux_input;

% sensor
sensor.mask = zeros(Nx, Ny);
sensor.mask(tdr_row : Nx - PML_size, PML_size + 1 : Ny - PML_size) = 1;
sensor.record = {'u', 'I'};

% simulation
input_args = {'DisplayMask', 'off', 'PlotScale', 'auto', 'PMLSize', PML_size, 'DataCast', 'single'};
sensor_data_shear = pstdElastic2D(kgrid, medium, source, sensor, input_args{:});

%% animation
% shear wave
time_step = kgrid.dt;
shear_velocity = sensor_data_shear.ux;
particle_pressure = rho * c * shear_velocity;
particle_intensity = sensor_data_shear.Iy;
frame_step = 400;

pressure_plot = reshape(particle_pressure(:, 20000), [260, 180]);
figure; 
imagesc((1 : 260) * dy * 1000, (1 : 180) * dx * 1000, pressure_plot./max(pressure_plot, [], 'all'));
% colormap jet;
% caxis([-3e-11, 3e-11]);
caxis([0, 1]);
axis image;
h = colorbar;
% colormap jet;
xlabel('[mm]');
ylabel('[mm]');
ylabel(h, 'Normalized Pressure');

figure;
plot((1 : 180) * dx * 1000, 20 * log10(pressure_plot(:, 112)./max(pressure_plot(:,112))));
ylim([-30 5]);
xlabel('[mm]');
ylabel('[dB]');

figure;
p_max_all = sensor_data_ARF.p_max_all;
q = 20 * log10(p_max_all(11 : 170, 100)./(max(p_max_all(11 : 170, 100), [], 'all')));
[up, low] = envelope(q, 100, 'peak');
x_vec = kgrid.x_vec;
plot(x_vec(11 : 170) * 1000 + 20, q);
hold on
plot(x_vec(11 : 170) * 1000 + 20, up);
hold off
ylim([-15, 2]);
xlabel('distance (mm)');
ylabel('pressure (dB)');
legend('signal', 'envelope');
title('compressional spl');
%% video
% v = VideoWriter('3MHz_besselbeam_no_attenuated.avi');
% v.FrameRate = 10;
% open(v);
% 
% for nn = 1 : frame_step : kgrid.Nt
%     
%     imagesc((1 : 180) * dy * 1000, (1 : 180) * dx * 1000, reshape(particle_pressure(:, nn), [180, 180]));
%     colormap hot;
%     axis image;
% %     caxis([5e-12, 3.36e-11]);
% %     h = colorbar;
%     xlabel('[mm]');
%     ylabel('[mm]');
% %     ylabel(h, '[Pa]');
%     pause(1)
%     
%     frame = getframe(gcf)
%     writeVideo(v, frame);
% 
% end
% 
% close(v);