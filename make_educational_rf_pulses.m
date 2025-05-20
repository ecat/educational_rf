out_folder = './pulses/';

addPaths()
%% make a hard pulse
N = 1000;
dT = 1 / 1e3;

max_B1s_G = [0.03, 0.06, 0.12]; % gives 45, 90, 180 hard pulse

for max_B1_for_fa_90_TBW_2_G = max_B1s_G
    
    B1_t = ones(N, 1) * max_B1_for_fa_90_TBW_2_G;
    B1_reim = cat(2, real(B1_t), imag(B1_t));
    
    tag = sprintf('hard_max_b1_%.2f', max_B1_for_fa_90_TBW_2_G);

    fid = fopen(sprintf('%s/%s_pulse.txt', out_folder, tag), 'w');
    fprintf(fid, '%f %f\n', B1_reim.');
    fclose(fid);
    
    fid = fopen(sprintf('%s/%s_dT.txt', out_folder, tag), 'w');
    fprintf(fid, '%f', dT);
    fclose(fid);
end

%% make some windowed sincs with equal bandwidth and different TBW
%max_B1_for_fa_90_TBW_2_G = sum(msinc(N, 2)) * dT * 4.258 * 90/360; % should be close to 0.21
%max_B1_for_fa_90_TBW_2_G = 0.237;
max_B1_for_fa_90_TBW_2_G = 0.237 / 2;
TBWs = [2, 4, 8];

for TBW = TBWs
    B1_t = msinc(N, TBW / 2).' * max_B1_for_fa_90_TBW_2_G;
    
    B1_reim = cat(2, real(B1_t), imag(B1_t));
    
    tag = sprintf('msinc_tbw_%d', TBW);
    
    fid = fopen(sprintf('%s/%s_pulse.txt', out_folder, tag), 'w');
    fprintf(fid, '%f %f\n', B1_reim.');
    fclose(fid);
    
    fid = fopen(sprintf('%s/%s_dT.txt', out_folder, tag), 'w');
    fprintf(fid, '%f', dT  * TBW / 2);
    fclose(fid);
end

%% truncated sinc
TBW = 4;
B1_t = msinc(N, TBW / 2).' * max_B1_for_fa_90_TBW_2_G;
% circshift and zero pad so that get same total duration so animations are synchronized
B1_t = circshift(B1_t, [250]);
B1_t(1:250) = 0;
B1_reim = cat(2, real(B1_t), imag(B1_t));

tag = sprintf('msinc_truncated_tbw_%d', TBW);

fid = fopen(sprintf('%s/%s_pulse.txt', out_folder, tag), 'w');
fprintf(fid, '%f %f\n', B1_reim.');
fclose(fid);

fid = fopen(sprintf('%s/%s_dT.txt', out_folder, tag), 'w');
fprintf(fid, '%f', dT  * TBW / 2);
fclose(fid);
%% make some multiband pulses
TBW = 4; 
N_bands = [1 2 4];
band_cycles_1 = 2 * 2 * pi;
band_cycles_2 = 3 * band_cycles_1;
for N_band = N_bands
    B1_t = msinc(N, TBW / 2) .' * max_B1_for_fa_90_TBW_2_G;
    if N_band == 2
        modulation = (cos(linspace(-band_cycles_1, band_cycles_1, N)))';
        modulation = modulation ./ max(abs(modulation(:)));
        B1_t = B1_t .* modulation;
    elseif N_band == 4
        modulation = (cos(linspace(-band_cycles_1, band_cycles_1, N)) + ...
            cos(linspace(-band_cycles_2, band_cycles_2, N)))';
        modulation = modulation ./ max(abs(modulation(:)));
        B1_t = B1_t .* modulation;
    end

    B1_reim = cat(2, real(B1_t), imag(B1_t));


    tag = sprintf('multiband_%d', N_band);
    
    fid = fopen(sprintf('%s/%s_pulse.txt', out_folder, tag), 'w');
    fprintf(fid, '%f %f\n', B1_reim.');
    fclose(fid);
    
    fid = fopen(sprintf('%s/%s_dT.txt', out_folder, tag), 'w');
    fprintf(fid, '%f', dT  * TBW / 2);
    fclose(fid);
end

%% msinc_highres
TBW = 8;
B1_t = msinc(4000, TBW / 2).' * max_B1_for_fa_90_TBW_2_G / 2;
B1_reim = cat(2, real(B1_t), imag(B1_t));

tag = sprintf('msinc_highres_tbw_%d', TBW);

fid = fopen(sprintf('%s/%s_pulse.txt', out_folder, tag), 'w');
fprintf(fid, '%f %f\n', B1_reim.');
fclose(fid);

fid = fopen(sprintf('%s/%s_dT.txt', out_folder, tag), 'w');
fprintf(fid, '%f', dT  * TBW / 2);
fclose(fid);

%% make some windowed sincs with different FA 
target_fas = [45, 90, 135, 180];
TBW = 8;
for target_fa = target_fas
    B1_t = msinc(N, TBW / 2).' * max_B1_for_fa_90_TBW_2_G * target_fa / 90;
    
    B1_reim = cat(2, real(B1_t), imag(B1_t));
    
    tag = sprintf('msinc_tbw_%d_fa_%d', TBW, target_fa);
    
    fid = fopen(sprintf('%s/%s_pulse.txt', out_folder, tag), 'w');
    fprintf(fid, '%f %f\n', B1_reim.');
    fclose(fid);
    
    fid = fopen(sprintf('%s/%s_dT.txt', out_folder, tag), 'w');
    fprintf(fid, '%f', dT  * TBW / 2);
    fclose(fid);
end

%% make some adiabatic inversions

b1s_adiabatic = [.06, .1, .2, .4];

for b1_max = b1s_adiabatic
    
    T_adiabatic = 10;
    f_max = 4;
    N = 1000;
    dT = T_adiabatic / N;
    [mag, omega, phi] = get_wurst_pulse(b1_max, N, T_adiabatic, f_max, 40);
    
    B1_t = mag.' .* exp(1i * phi.');
    
    tag = sprintf('wurst_b1_%.2f', b1_max);
    
    B1_reim = cat(2, real(B1_t), imag(B1_t));
    
    fid = fopen(sprintf('%s/%s_pulse.txt', out_folder, tag), 'w');
    fprintf(fid, '%f %f\n', B1_reim.');
    fclose(fid);
    
    fid = fopen(sprintf('%s/%s_dT.txt', out_folder, tag), 'w');
    fprintf(fid, '%f', dT);
    fclose(fid);
end