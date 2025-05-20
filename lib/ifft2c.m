function y = ifft2c(x)
    y = cifft(cifft(x, 2), 1);