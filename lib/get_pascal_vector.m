function v = get_pascal_vector(N)


p = pascal(N);


v = zeros(N, 1);

for nn = 1:N
    v(nn) = p(end - nn + 1, nn);
end