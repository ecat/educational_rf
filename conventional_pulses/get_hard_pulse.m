function [rf, duration] = get_hard_pulse(max_b1, Nt, flip_angle)
% max_b1 in [G] 
% flip_angle in degrees

gamma = 4.258; % kHz / G
target_cycles = flip_angle / 360;

duration = target_cycles/ (gamma * max_b1);

rf = ones(Nt, 1) * max_b1;

end

