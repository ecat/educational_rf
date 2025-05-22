function [mag, omega, phi] = get_wurst_pulse(b1_max, nt, Tp, f_max, N)
% Equation 4, 5 in https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2688743/
% b1_max units don't matter
% nt integer
% Tp in ms
% f_max in kHz
% N unitless

dt = Tp/nt;
t = linspace(0, Tp, nt);
tau = 2 * t / Tp - 1;
mag = b1_max * (1 - abs(cos(pi * (t/Tp)).^N));

omega = 2 * pi * f_max * tau; % krad/s
phi = cumsum(omega) * dt;

%phi = 2 * pi * f_max * log(sech(beta * tau)); % krad
%omega = [0, diff(phi)] / dt;

