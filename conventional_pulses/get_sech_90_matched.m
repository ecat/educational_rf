function [mag, omega, phi] = get_sech_90_matched(scheme, b1_max, nt, Tp, f_max, beta)
% input parameters that would give a nice 180
% schemes follow https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2688743/

if(scheme == 1)
    [mag, omega, phi] = get_sech_pulse(b1_max, nt, Tp * 2, f_max, beta);
elseif(scheme == 2)
    error('cant change slice select')
elseif (scheme == 3)
    [mag, omega, phi] = get_sech_pulse(b1_max, nt, Tp, f_max, beta/2);
end

end