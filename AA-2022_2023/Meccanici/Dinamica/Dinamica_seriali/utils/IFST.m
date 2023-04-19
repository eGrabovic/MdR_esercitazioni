function [p, p_i, p_ip] = IFST(x, pk, m, order)

p = zeros(2,1, class(x));
p_i = zeros(2, 2*m+1, class(x));
p_ip = zeros(2, 2*m+1, class(x));
counter = 1;
% for k = -m:1:m
%     id = k+m+1;
%     p_ip(:, counter) = p;
%     p = p + rotZ2D(k.*x)*pk(:, id);
%     p_i(:, counter) = p;
%     counter = counter+1;
% end

for id = order
    k = id - m - 1;
    p_ip(:, counter) = p;
    p = p + rotZ2D(k.*x)*pk(:, id);
    p_i(:, counter) = p;
    counter = counter+1;
end
