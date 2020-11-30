% filtering data
% [sig_1_filtered_LR, sig_2_filtered_LR] = filter_data(sig1_int,sig2_int);

% signals are given
% field
ROI_length = 0.04; % [m] 5cm
ROI_depth = 0.04; % [m] 5cm
PML_size = 10; % [grid]
ROI_length_size = 180; % [grid] 148
ROI_depth_size = 180; % [grid]
Nx = ROI_depth_size + 2 * PML_size; % [grid]
Ny = ROI_length_size + 2 * PML_size; % [grid]
dx = ROI_depth / ROI_depth_size; % [m]
dy = ROI_length / ROI_length_size; % [m]

% calculation
[x, y, z] = size(sig1_int);
fr = y; % time steps
t_array_step = (6.5e-8)*3; % [sec]
Ts = t_array_step * fr; % [sec] displacement time
ensemble_number = 45;
xx = (0 : fr-1)./fr .* Ts; % no. of displacement points 
stepx = xx(2) - xx(1);
lags = -(length(xx) - 1) : (length(xx) - 1);

% ROI_depth_size = 90;

for ttt = 1 : ROI_depth_size;
    for mmm = 1 : ensemble_number;
    
%     [C(ttt, :, mmm), L(ttt, :, mmm)] = xcorr(sig_2_filtered_LR(ttt,:, mmm), sig_1_filtered_LR(ttt,:, mmm), 'coeff');
    [C(ttt, :, mmm), L(ttt, :, mmm)] = xcorr(sig2_int(ttt,:, mmm), sig1_int(ttt,:, mmm), 'coeff');

    end
end


[pval, ploc] = max(C(:, :, 1 : ensemble_number), [], 2);
delay = lags(ploc(:, :, 1 : ensemble_number));

delay = 1 .* permute(delay(:, 1, :), [1 3 2]) .* stepx; % time delay for two respective displacements [sec]
%% image
5 * dy./delay;
figure(5); imagesc((1:ensemble_number)*dy*1000, (1:180)*dx*1000, ans); colormap jet; h = colorbar; axis image; ylabel(h, '[m/s]');xlabel('[mm]');ylabel('[mm]');
