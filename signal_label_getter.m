function [label] = signal_label_getter(signal)

signal_parsed = signal;
if(abs(real(signal_parsed)) < 0.01)
    % add special conditions so that avoids toggling between - and + in
    % animation
    signal_parsed = abs(real(signal_parsed)) + 1i * imag(signal_parsed);
end

if(abs(imag(signal_parsed)) < 0.01)
    signal_parsed = real(signal_parsed) + 1i * abs(imag(signal_parsed));
end



if(real(signal_parsed) > 0 && imag(signal_parsed) > 0)
    label = sprintf('Total Signal $=%.2f + %.2fj$', abs(real(signal_parsed)), abs(imag(signal_parsed)));
elseif(real(signal_parsed) < 0 && imag(signal_parsed) > 0)
    label = sprintf('Total Signal $=-%.2f + %.2fj$', abs(real(signal_parsed)), abs(imag(signal_parsed)));
elseif(real(signal_parsed) < 0 && imag(signal_parsed) < 0)
    label = sprintf('Total Signal $=-%.2f - %.2fj$', abs(real(signal_parsed)), abs(imag(signal_parsed)));
else
    label = sprintf('Total Signal $=%.2f - %.2fj$', abs(real(signal_parsed)), abs(imag(signal_parsed)));
end

end

