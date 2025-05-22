
set(groot, 'defaulttextinterpreter', 'latex');
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'DefaultLineLineWidth', 1.5);
%%
%	Modified Philip Lee 2025
%
anim_step = 4; % animstep from 0, 1, 2, 3, 4

x = -10:.025:10;		% cm, MANY spins!
z = -5.5:.025:5.5;	

Nsl=11;				% #slices
if(anim_step >= 4)
    Slthick=4*max(z)/Nsl;		% cm, slice thickness
else
    Slthick=2*max(z)/Nsl;		% cm, slice thickness
end
xres = 0.4;			% cm
Nx = 20/xres;
bwpix = 0.5;			% kHz
BWRF=1.5;				% kHz, RF bandwidth
gamma = 4.258;			% kHz/G, Gyromagnetic ratio
Gz = BWRF/(gamma*Slthick);	% G/cm, Slice-select gradient


[xx,zz] = meshgrid(x,z);
rr = sqrt(xx.^2+zz.^2);

f=0*xx;	% -- Allocate frequencies.



%%
fgradz = f + zz * gamma * Gz;	% Frequency with gradient.
fnograd = f;
bw = .5;
bh = .5;


%% make zf plot

titles = {'No Slice Select Gradient', 'Slice Select Gradient On', ...
    'Slice Select Gradient On', 'Slice Select Gradient On', 'Smaller Gradient Amplitude'};

orange_color = max([0.8500, 0.3250, 0.0980] - .1, 0);
yellow_color = [0.9290, 0.6940, 0.1250];
blue_color =  	[0.2010, 0.6450, 0.8330] - .2;
indices_to_plot = 1:5:numel(zz);

offres_cutoff = .1; %[kHz]
offres_indices = intersect(find(abs(fnograd) > offres_cutoff), indices_to_plot);
onres_indices = intersect(find(abs(fnograd) < offres_cutoff), indices_to_plot);

slice_center = 0.1;
slice_thickness = 0.4;
excited_onres_indices = intersect(onres_indices, find(abs(zz - slice_center) < slice_thickness));

fig_slice_select_plot = figure('Position', [100 100 800 600], 'Color', 'white')

if(anim_step == 0)
    scatter(zz(onres_indices), -fnograd(onres_indices), 'MarkerEdgeColor', blue_color, 'MarkerFaceColor', blue_color); hold on;
    scatter(zz(excited_onres_indices), -fnograd(excited_onres_indices), 'MarkerEdgeColor', blue_color, 'MarkerFaceColor', blue_color);
    hlegend = legend('Spins', 'Location', 'south east');
elseif(anim_step >= 1)
    scatter(zz(onres_indices), fgradz(onres_indices), 'MarkerEdgeColor', blue_color, 'MarkerFaceColor', blue_color); hold on;
    hlegend = legend('Spins', 'Excited Spins', 'Location', 'south east');
end
hlegend.AutoUpdate = 'off';
hlegend.FontSize = 20;
%scatter(zz(offres_indices), -fnograd(offres_indices), 'MarkerEdgeColor', yellow_color); 
ax = gca;
ax.YAxis.FontSize = 20;
ax.XAxis.FontSize = 20;
%hlegend = legend('On-Resonance Spins');
%hlegend.FontSize = 14;
ylabel('Frequency [kHz]', 'FontSize', 24)
xlabel('z [cm]', 'FontSize', 24)
xlim_width = 12;
ylim_width = 16;
ylim([-ylim_width/2, ylim_width/2])
xlim([-xlim_width/2, xlim_width/2])
grid on; grid minor; hold on;
title(titles{anim_step + 1}, 'FontSize', 28)


% draw slice select gradient
if(anim_step >= 2)
    
    if(anim_step == 2)
        patch_alphas = [.1];
        bin_mask_sl_to_plot = [0];
    else
        patch_alphas = [.3, .2, .1];
        bin_mask_sl_to_plot = [2, 1, 0]; % make sure 0 is last so that we can lineup with labels
    end
    for ii = 1:numel(bin_mask_sl_to_plot)
        bin_mask_sl = bin_mask_sl_to_plot(ii);
        x1 = linspace(min(z) * 1.5, max(z) * 1.5, 100);
        %slope = Gz * gamma;
        slope = 0;
        y1 = slope * x1 + bin_mask_sl * BWRF + BWRF/2;
        y2 = slope * x1 + bin_mask_sl * BWRF - BWRF/2;

        plot(x1, y1, '--', 'Color', 'k'); hold on;
        plot(x1, y2, '--', 'Color', 'k');

        patch_xcoords = [x1(1) x1(1) x1(end) x1(end)];
        patch_ycoords = [y1(1) y2(1) y2(end) y1(end)];
        patch(patch_xcoords, patch_ycoords, 'k', 'FaceAlpha', patch_alphas(ii));

        x_index = 25;
        %line_angle = 90 -atand(slope * xlim_width/ylim_width); % need to account for aspect ratio of the plot
        line_angle = 0;
        text_x_loc = x1(x_index);
        text_y_loc = (y1(x_index) + y2(x_index))/2 - .1;
        text(text_x_loc, text_y_loc, sprintf('Slice %d', bin_mask_sl + 1), ...
            'FontSize', 16, 'Rotation', line_angle, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    end
end

if(anim_step >= 1)
    % draw a slope annotation
    index_to_place_slope_marker_1 = 150;
    index_to_place_slope_marker_2 = 71;

    x_offset = [0.1 0];
    if(anim_step <= 3)
        offset1 = 1.5;
        offset2 = -.8;
    else
        offset1 = -.2;
        offset2 = -1;
    end
    p1_brace_slope = [zz(index_to_place_slope_marker_1) fgradz(index_to_place_slope_marker_1)] + x_offset;
    p2_brace_slope = [zz(index_to_place_slope_marker_2) (fgradz(index_to_place_slope_marker_2))] + x_offset;


    hbrace = drawbrace( p1_brace_slope, ...
        p2_brace_slope);
    hbrace.Color = 'k';

    text((p1_brace_slope(1) + p2_brace_slope(2))/2 + offset1, (p2_brace_slope(2) + p1_brace_slope(2))/2 + offset2, ...
        '$f = \gamma{}G_{s}z$ ', 'Interpreter', 'latex', 'HorizontalAlignment', 'left', 'FontSize', 20)

end

if(anim_step >= 2)
    % draw a slab width annotation
    if(anim_step <= 3)
        index_to_place_bw_marker = 47;
        brace_width = 1.1;
        offset1 = 0;
        offset2 = 0;
    else
        index_to_place_bw_marker = 45;
        brace_width = 1.9;
        offset1 = .3;
        offset2 = .8;
    end

    p1_brace_dz =[x1(index_to_place_bw_marker) y2(index_to_place_bw_marker)];
    p2_brace_dz =[x1(index_to_place_bw_marker) + brace_width, y2(index_to_place_bw_marker)];

    plot([x1(index_to_place_bw_marker) x1(index_to_place_bw_marker)], [y1(index_to_place_bw_marker) y2(index_to_place_bw_marker)], 'k:');
    plot([x1(index_to_place_bw_marker) x1(index_to_place_bw_marker)] + brace_width, [y1(index_to_place_bw_marker) y2(index_to_place_bw_marker)], 'k:');

    hbrace = drawbrace(p2_brace_dz, p1_brace_dz);
    hbrace.Color = 'k';
    text((p1_brace_dz(1) + p2_brace_dz(1))/2 - .1, (p2_brace_dz(2) + p1_brace_dz(2))/2 - .5 - offset1, ...
        {'Slice Thickness $\Delta z = \frac{\mathrm{RF}_\mathrm{BW}}{\gamma G_s}$'},...
        'Interpreter', 'latex', 'HorizontalAlignment', 'left', 'FontSize', 20, 'VerticalAlignment', 'top', 'Color', orange_color)
    
    % draw RF BW annotation
    x_offset = [1.2 + offset2 0];
    p1_brace_rfbw = [x1(index_to_place_bw_marker) y2(index_to_place_bw_marker)+BWRF] + x_offset;
    p2_brace_rfbw = [x1(index_to_place_bw_marker) (y2(index_to_place_bw_marker))] + x_offset;

    plot([x1(index_to_place_bw_marker) p1_brace_rfbw(1)], [p1_brace_rfbw(2) p1_brace_rfbw(2)], 'k:');
    plot([x1(index_to_place_bw_marker) p2_brace_rfbw(1)], [p2_brace_rfbw(2) p2_brace_rfbw(2)], 'k:');

    % draw RF brace
    hbrace = drawbrace( p1_brace_rfbw, p2_brace_rfbw);
    hbrace.Color = 'k';
    text(p1_brace_rfbw(1) + .3, (p2_brace_rfbw(2) + p1_brace_rfbw(2))/2 - .1, 'RF Bandwidth', ...
        'Interpreter', 'latex', 'HorizontalAlignment', 'left', 'FontSize',20, 'VerticalAlignment', 'middle', 'Color', orange_color)
end

ylim([-ylim_width/2, ylim_width/2])
xlim([-xlim_width/2, xlim_width/2])
%
do_save = 1;
if(do_save)
    export_fig(fig_slice_select_plot, sprintf('out2/slice_select_explain_step%d.png', anim_step), '-m2.');
end