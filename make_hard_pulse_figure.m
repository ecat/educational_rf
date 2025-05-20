

my_pulse = msinc(256, 2);
my_pulse = my_pulse / max(abs(my_pulse(:)));

N_hard_pulses = 16;

b_amps = my_pulse(1: floor(256 / N_hard_pulses) : end);

padding = .4; 
b1_width = 1;
g_width = 1;
grad_amps = .2;

gradient_color = '#FBE5D6';
rf_color = '#8C1515';
figure('Position', [100 100 1200 400], 'Color', 'white'); 

for tt = 1:N_hard_pulses
    
    t1 = (tt - 1) * (b1_width + g_width + padding * 2);
    t2 = t1 + b1_width; 
    t3 = t2 + padding;
    t4 = t3 + g_width;

    y1 = b_amps(tt);
    y2 = grad_amps;

    if(y1 > 0)
        yloc = 0;
    else
        yloc = y1;
    end

    r1 = rectangle('Position', [t1 yloc (t2 - t1) abs(y1)], 'FaceColor', rf_color, 'EdgeColor', 'w', 'LineWidth', 1); hold on;

    if(y2 > 0)
        yloc = 0;
    else
        yloc = y2;
    end
    r2 = rectangle('Position', [t3 yloc (t4 - t3) abs(y2)], 'FaceColor', gradient_color, 'EdgeColor', 'w', 'LineWidth', 1);

end


%legend([r1, r2], 'B1', 'Gradient', 'AUtoUpdate', 'off', 'FontSize', 18);
ax = gca;
plot([0, t3 + 1], [0 0], 'k')
axis off;
box off;