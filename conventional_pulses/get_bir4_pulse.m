function [mag, phi] = get_bir4_pulse(b1_max, nt, dt, w_max, fa, beta, tan_kappa)
% b1_max in Gauss
Tp = nt * dt;
time = (0:1:nt-1) .* dt;

kappa = atan(tan_kappa);
fb = zeros(nt, 1);
fw = zeros(nt, 1);
dphi = pi + fa / 2;

for t_section = 1:4
    ind_section = ((t_section-1)*nt/4+1):(t_section*nt/4);
    tao = time(ind_section) / (Tp / 4);
    if t_section == 1
        fb(ind_section) = tanh(beta*(1-tao));
        fw(ind_section) = tan(kappa * tao) / tan_kappa;
    elseif t_section == 2
        fb(ind_section) = tanh(beta*(tao-1)) * exp(1i * dphi);
        fw(ind_section) = tan(kappa * (tao-2)) / tan_kappa;
    elseif t_section == 3
        fb(ind_section) = tanh(beta*(3-tao)) * exp(1i * dphi);
        fw(ind_section) = tan(kappa * (tao-2)) / tan_kappa;
    elseif t_section == 4
        fb(ind_section) = tanh(beta*(tao-3));
        fw(ind_section) = tan(kappa * (tao-4)) / tan_kappa;
    end
end

phi = cumsum(w_max * fw * dt);
phi((nt/4+1):(nt*3/4)) = phi((nt/4+1):(nt*3/4)) + dphi;
mag = abs(fb) * b1_max;


%% plot pulse waveform
time = (0:1:nt-1) .* dt;
figure
subplot(2,1,1)
plot(time*1e3, mag)
xlabel('time (ms)')
xlim([0 Tp*1e3])
title('amplitude (G)')
subplot(2,1,2)
plot(time*1e3, phi)
xlabel('time (ms)')
xlim([0 Tp*1e3])
title('phase (rad)')

end