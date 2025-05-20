addPaths()

%%
do_write_video = 1;

crop_height_vals = 1:300;
crop_width_vals = 1:400;
if(do_write_video)
    video_writer_helper = @ (v, fig_handle)  writeVideo(v, crop_getframe(getframe(fig_handle), crop_height_vals, crop_width_vals));
else
    video_writer_helper = @(a, b) 1/1;
end

out_path = 'out2/';
out_path_helper = @(x) sprintf('out2/%s', x);

if(~exist(out_path, 'dir'))
    mkdir(out_path)
end
%%

rotation_axis_phases = [];
rotation_axis_colorstr = {'k'};

get_fig = @() figure('Position', [100 100 400 400], 'Color', 'white');
%% precession
v = VideoWriter(out_path_helper('precession'));
v.FrameRate = 15;
open(v)

fig_tmp = get_fig();
ax = axes; 
mag1_smt = (yrot(30) * [0 0 1]')';

N_t = 120;

for tt = 1:N_t
    cla(ax);

    rotation_vector_mt = repmat([0; 0; 1], [1 tt]);
    mag1_to_plot = cat(1, mag1_smt, reshape(rotation_vector_mt, [1 size(rotation_vector_mt)]));
    ax = animate_mag_matrix( mag1_to_plot, rotation_axis_phases, rotation_axis_colorstr, 1, [], {'f3', 'i1'});
    title('Precession About z', 'FontSize', 20)
    video_writer_helper(v, fig_tmp)


    cur_magnetization = mag1_smt(:, :, end)';
    magnetization_angle = angle(cur_magnetization(1) + 1i * cur_magnetization(2)) * 180/pi;

    R = zrot(-360 * 4 / N_t);
    mag1_smt = cat(3, mag1_smt, (R * mag1_smt(:, :, end)')');
end

close(v)
%% laboratory frame
v = VideoWriter(out_path_helper('laboratory_frame'));
v.FrameRate = 30;
open(v)

fig_tmp = get_fig();
ax = axes; 
mag1_smt = [0 0 1];
magnetization_angle = 0;
B_1_eff = .707 * [cosd(magnetization_angle + 90 + 25.7143 * 1); sind(magnetization_angle + 90 + 25.7143 * 1); 0];
B_0_eff = [0; 0; .707];

B_eff_smt = cat(1, reshape(B_1_eff, [1 3]), reshape(B_0_eff, [1 3]));

N_t = 210;

for tt = 1:N_t
    cla(ax);

    to_plot_smt = cat(1, mag1_smt, B_eff_smt);
    ax = animate_mag_matrix( to_plot_smt, rotation_axis_phases, rotation_axis_colorstr, 1, [], {'f3', 'i1', 'i1'});
    title('Excitation in Laboratory Frame', 'FontSize', 20)
    video_writer_helper(v, fig_tmp)

    cur_magnetization = mag1_smt(:, :, end)';
    magnetization_angle = angle(cur_magnetization(1) + 1i * cur_magnetization(2)) * 180/pi;

    N_cycles_larmor = 15;
    R = zrot(-360 * N_cycles_larmor / N_t) * throt(-90 / N_t, magnetization_angle + 90);
    mag1_smt = cat(3, mag1_smt, (R * mag1_smt(:, :, end)')');

    

    B_1_eff = .707 * [cosd(magnetization_angle + 90 + 25.7143); sind(magnetization_angle + 90 + 25.7143); 0];
    B_0_eff = [0; 0; .707];
    B_eff_sm = cat(1, reshape(B_1_eff, [1 3]), reshape(B_0_eff, [1 3]));
    B_eff_smt = cat(3, B_eff_smt, reshape(B_eff_sm, [2 3 1]));
end

% write last frame
cla(ax)
to_plot_smt = cat(1, mag1_smt, B_eff_smt);
ax = animate_mag_matrix( to_plot_smt, rotation_axis_phases, rotation_axis_colorstr, 1, [], {'f3', 'i1', 'i1'});
title('Excitation in Laboratory Frame', 'FontSize', 20)
video_writer_helper(v, fig_tmp)

N_pause = 60;
for tt = 1:N_pause
    video_writer_helper(v, fig_tmp)
end

close(v)

%% exc
fig_tmp = get_fig();
ax = axes;
v = VideoWriter(out_path_helper('excitation'));
v.FrameRate = 15;
open(v)
mag1_smt = [0 0 1];
B_eff_smt = reshape([0 0.707 0], [1 3 1]);

N_t = v.FrameRate * 7;

for tt = 1:N_t
    cla(ax);
    to_plot_smt = cat(1, mag1_smt, B_eff_smt);
    ax = animate_mag_matrix( to_plot_smt, rotation_axis_phases, rotation_axis_colorstr, 1, [], {'f3', 'i1'});
    title('Excitation in Rotating Frame', 'FontSize', 20)
    video_writer_helper(v, fig_tmp)
    mag1_smt = cat(3, mag1_smt, (yrot(-90 / N_t) * mag1_smt(:, :, end)')');
    B_eff_smt = cat(3, B_eff_smt, reshape([0 1 0], [1 3 1]));
end

% write last frame
cla(ax)
to_plot_smt = cat(1, mag1_smt, B_eff_smt);
ax = animate_mag_matrix( to_plot_smt, rotation_axis_phases, rotation_axis_colorstr, 1, [], {'f3', 'i1'});
title('Excitation in Rotating Frame', 'FontSize', 20)
video_writer_helper(v, fig_tmp)
N_pause = 30;
for tt = 1:N_pause
    video_writer_helper(v, fig_tmp)
end


close(v)

%% excitation with different beff
b_eff_offres = [-.9 0 .9];
b1_amp = .4; % so that scaling between offres and b1 are correct
offres_degrees = b_eff_offres ./ b1_amp * 90;
titles = {'Negative Off-Resonance', 'On-Resonant', 'Positive Off-Resonance'};

for ii = 1:numel(b_eff_offres)
    fig_tmp = get_fig();
    ax = axes;
    
    v = VideoWriter(out_path_helper(sprintf('excitation_with_gradient_%d', ii)));
    v.FrameRate = 15;
    open(v)
    mag1_smt = [0 0 1];
    if(ii == 2)
        B_eff_init = reshape([-b1_amp 0 0], [1 3]);
    else
        B_eff_init = reshape([ -b1_amp 0 b_eff_offres(ii); -b1_amp 0 0; 0 0 b_eff_offres(ii)], [3 3]);        
    end
    B_eff_smt = B_eff_init;
    
    N_t = v.FrameRate * 5;
    
    for tt = 1:N_t
        cla(ax);
        to_plot_smt = cat(1, mag1_smt, B_eff_smt);
        ax = animate_mag_matrix( to_plot_smt, rotation_axis_phases, rotation_axis_colorstr, 1, [], {'f3', 'i3', 'z2', 'Q2'});
        title(titles{ii}, 'FontSize', 20)
        video_writer_helper(v, fig_tmp)
        mag1_smt = cat(3, mag1_smt, (zrot(-offres_degrees(ii) / N_t) * xrot(90 / N_t) * mag1_smt(:, :, end)')');
        B_eff_smt = cat(3, B_eff_smt, B_eff_init);
    end

    N_pause = v.FrameRate * 1.5;
    for tt = 1:N_pause
        cla(ax);
        to_plot_smt = cat(1, mag1_smt, B_eff_smt);
        ax = animate_mag_matrix( to_plot_smt, rotation_axis_phases, rotation_axis_colorstr, 1, [], {'f3', 'i3', 'z2', 'Q2'});
        title(titles{ii}, 'FontSize', 20)
        video_writer_helper(v, fig_tmp)
    end
    
    close(v)
end

%% inversion
fig_tmp = get_fig();
ax = axes; 
mag1_smt = [0 0 1; 0 .707 0];

v = VideoWriter(out_path_helper('inversion'));
v.FrameRate = 15;
open(v)

N_t = 45;

for tt = 1:N_t
    cla(ax);
    ax = animate_mag_matrix( mag1_smt, rotation_axis_phases, rotation_axis_colorstr, 1, [], {'f3', 'i1'});
    title('Inversion', 'FontSize', 20)
    video_writer_helper(v, fig_tmp)
    mag1_smt = cat(3, mag1_smt, (yrot(-180 / N_t) * mag1_smt(:, :, end)')');
end

N_pause = 30;
cla(ax);
ax = animate_mag_matrix( mag1_smt, rotation_axis_phases, rotation_axis_colorstr, 1, [], {'f3', 'i1'});
title('Inversion', 'FontSize', 20)
for tt = 1:N_pause
    video_writer_helper(v, fig_tmp)
end

close(v)
%% excitation 2 same timing as inversion
fig_tmp = get_fig();
ax = axes; 
mag1_smt = [0 0 1; 0 .707 0];

v = VideoWriter(out_path_helper('excitation2'));
v.FrameRate = 15;
open(v)

N_t = 45;

for tt = 1:N_t
    cla(ax);
    ax = animate_mag_matrix( mag1_smt, rotation_axis_phases, rotation_axis_colorstr, 1, [], {'f3', 'i1'});
    title('Excitation', 'FontSize', 20)
    video_writer_helper(v, fig_tmp)
    mag1_smt = cat(3, mag1_smt, (yrot(-90 / N_t) * mag1_smt(:, :, end)')');
end

N_pause = 30;
cla(ax);
ax = animate_mag_matrix( mag1_smt, rotation_axis_phases, rotation_axis_colorstr, 1, [], {'f3', 'i1'});
title('Excitation', 'FontSize', 20)
for tt = 1:N_pause
    video_writer_helper(v, fig_tmp)
end
close(v)

%% refocusing
for refocusing_rotation_option = [1 2]

    fig_tmp = get_fig();
    ax = axes; 
    
    if(refocusing_rotation_option == 1)
        f_rotation = @(x) yrot(x);
        refocusing_title = 'Refocusing $180_y$';
        refocusing_axis_sm = [0 0.707 0];
    else
        f_rotation = @(x) xrot(x);
        refocusing_title = 'Refocusing $180_x$';
        refocusing_axis_sm = [0.707 0 0];
    end
    
    mag1_smt = [1 0 0; 1 0 0; 1 0 0; 1 0 0 ; 1 0 0];
    N_spins = size(mag1_smt, 1);
    final_offres_angles = [-80 -40 0 40 80];
    
    % initial
    cla(ax);
    ax = animate_mag_matrix(mag1_smt, rotation_axis_phases, rotation_axis_colorstr, 1, [], {});
    
    v = VideoWriter(out_path_helper(sprintf('refocusing_%d_step0', refocusing_rotation_option)));
    v.FrameRate = 15;
    open(v)
    
    % offresonance
    N_t = 20;
    for tt = 1:N_t
        cla(ax);
    
        ax = animate_mag_matrix( mag1_smt, rotation_axis_phases, rotation_axis_colorstr, 1, [], {});
        title(refocusing_title, 'FontSize', 20)
        video_writer_helper(v, fig_tmp)
        to_cat_sm = zeros(size(mag1_smt, [1 2]));
        for ss = 1:N_spins
            to_cat_sm(ss, :) = (zrot(final_offres_angles(ss) / N_t) * reshape(mag1_smt(ss, :, end), [3 1]))';
        end
        mag1_smt = cat(3, mag1_smt, reshape(to_cat_sm, [N_spins 3 1]));
    end
    
    close(v)
    
    
    v = VideoWriter(out_path_helper(sprintf('refocusing_%d_step1', refocusing_rotation_option)));
    v.FrameRate = 15;
    open(v)
    N_t = 60;
    for tt = 1:N_t

        to_append = repmat(refocusing_axis_sm, [1 1 size(mag1_smt, 3)]);

        cla(ax);
        ax = animate_mag_matrix(cat(1, mag1_smt, to_append), rotation_axis_phases, ...
            rotation_axis_colorstr, 1, [], {'f3', 'f3', 'f3', 'f3', 'f3', 'i1'});
        title(refocusing_title, 'FontSize', 20)
    
        video_writer_helper(v, fig_tmp)
        mag1_smt = cat(3, mag1_smt, (f_rotation(-180 / N_t) * mag1_smt(:, :, end)')');
    end
    
    
    % offresonance
    N_t = 20;
    for tt = 1:N_t
        cla(ax);
    
        ax = animate_mag_matrix( mag1_smt, rotation_axis_phases, rotation_axis_colorstr, 1, [], {});
    
        title(refocusing_title, 'FontSize', 20)
        video_writer_helper(v, fig_tmp)
        to_cat_sm = zeros(size(mag1_smt, [1 2]));
        for ss = 1:N_spins
            to_cat_sm(ss, :) = (zrot(final_offres_angles(ss) / N_t) * reshape(mag1_smt(ss, :, end), [3 1]))';
        end
        mag1_smt = cat(3, mag1_smt, reshape(to_cat_sm, [N_spins 3 1]));
    end
    cla(ax);

    ax = animate_mag_matrix( mag1_smt, rotation_axis_phases, rotation_axis_colorstr, 1, [], {});

    title(refocusing_title, 'FontSize', 20)
    video_writer_helper(v, fig_tmp)
    close(v)
end
%% adiabatic
ax = axes; 

b1_maxes = [.4, .11];
titles = {{'Adiabatic Inversion'}, {'Adiabatic Condition Unmet'}};
for ii = 1:numel(b1_maxes)
    mag1_smt = [0 0 1];
    b1_max = b1_maxes(ii);
    N_t = 200;
    
    T_adiabatic = 10;
    f_max = 2;
    dT = T_adiabatic / N_t;
    [mag, omega, phi] = get_sech_pulse(b1_max, N_t, T_adiabatic, -f_max, 4);
    
    B1_eff_mt = [mag * 4.258;
        zeros(1, N_t);    
        omega / (2 * pi)];
    
    v = VideoWriter(out_path_helper(sprintf('adiabatic_b1max_%.2f', b1_max)));
    v.FrameRate = 25;
    open(v)
    fig_tmp = get_fig();
    ax = axes; 
    
    alpha_t = sos(B1_eff_mt, 1) * dT * 360;
    zeta_t = -(atan2(B1_eff_mt(3, :), B1_eff_mt(1, :)) * 180/pi - 90);
    
    for tt = 1:N_t - 1
        cla(ax);
        % dwight figure 2.5
        theta = 0;
        alpha = alpha_t(tt);
        zeta = zeta_t(tt);
        R = zrot(-theta) * yrot(-zeta) * zrot(alpha) * yrot(zeta) * zrot(theta);
    
        mag1_smt = cat(3, mag1_smt, (R * mag1_smt(:, :, end)')');
    
        mag_to_plot_smt = cat(1, mag1_smt, reshape(B1_eff_mt(:, 1:tt + 1), [1 3 tt + 1]) / max(abs(B1_eff_mt(:))));
        ax = animate_mag_matrix( mag_to_plot_smt, rotation_axis_phases, rotation_axis_colorstr, 1, [], {'f3', 'i1'});
        ht = title(titles{ii}, 'FontSize', 20);

        video_writer_helper(v, fig_tmp)
    end
    
    close(v)
end

%% binomial pulse
titles = {'Binomial On-Resonance', 'Binomial Off-Resonance'};
offres_rotations = [0, 180];
for ii = 1:2

    offres_rot = offres_rotations(ii);

    fig_tmp = get_fig();
    ax = axes; 
    mag1_smt = [0 0 1];
    
    v = VideoWriter(out_path_helper(sprintf('binomial_offres_rot_%d_step_0', offres_rot)));
    v.FrameRate = 15;
    open(v)
    
    N_t = 26;
    
    for tt = 1:N_t
        cla(ax);
        ax = animate_mag_matrix( mag1_smt, rotation_axis_phases, rotation_axis_colorstr, 1, [], {});
        title(titles{ii}, 'FontSize', 20)
        video_writer_helper(v, fig_tmp)
        mag1_smt = cat(3, mag1_smt, (yrot(-45 / N_t) * mag1_smt(:, :, end)')');
    end

    for tt = 1:N_t
        cla(ax);
        ax = animate_mag_matrix( mag1_smt, rotation_axis_phases, rotation_axis_colorstr, 1, [], {});
        title(titles{ii}, 'FontSize', 20)
        video_writer_helper(v, fig_tmp)
        mag1_smt = cat(3, mag1_smt, (zrot(-offres_rot / N_t) * mag1_smt(:, :, end)')');
    end

    % include last frame
    cla(ax);
    ax = animate_mag_matrix( mag1_smt, rotation_axis_phases, rotation_axis_colorstr, 1, [], {});
    title(titles{ii}, 'FontSize', 20)
    video_writer_helper(v, fig_tmp)

    close(v)

    v = VideoWriter(out_path_helper(sprintf('binomial_offres_rot_%d_step_1', offres_rot)));
    v.FrameRate = 15;
    open(v)

    for tt = 1:N_t
        cla(ax);
        ax = animate_mag_matrix( mag1_smt, rotation_axis_phases, rotation_axis_colorstr, 1, [], {});
        title(titles{ii}, 'FontSize', 20)
        video_writer_helper(v, fig_tmp)
        mag1_smt = cat(3, mag1_smt, (yrot(-45 / N_t) * mag1_smt(:, :, end)')');
    end

    % include last frame
    cla(ax);
    ax = animate_mag_matrix( mag1_smt, rotation_axis_phases, rotation_axis_colorstr, 1, [], {});
    title(titles{ii}, 'FontSize', 20)
    video_writer_helper(v, fig_tmp)
    
    close(v)
end


%% plot some static vectors for slr example
fig_tmp = get_fig();
ax = axes; 
mag1_smt = [0 0 1; 0 1 0];

ax = animate_mag_matrix( mag1_smt, rotation_axis_phases, rotation_axis_colorstr, 1, [], {'i3', 'Q3'}, 0);

export_fig(fig_tmp, 'out2/basis_vectors.png')