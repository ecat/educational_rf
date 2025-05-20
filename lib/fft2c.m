function y = fft2c(x)
y = cfft(cfft(x, 2), 1);