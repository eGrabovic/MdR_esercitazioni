function pk = DFT_vec(p, m)
% 
% p = discrete points
% m = size of series coefficients
N = size(p, 2);
pk = zeros(2, 2*m+1);

for k = -m:1:m
    id = k+m+1;

    for j = 1:N
        pk(:,id) = pk(:,id) + rotZ2D(-2*pi*k*(j)./N)*p(:, j);
    end
end
pk = pk./N;
end