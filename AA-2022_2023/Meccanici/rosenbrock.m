function out = rosenbrock(x)
x = x(:);
out = sum(100.*(x(2:end) - x(1:end-1).^2).^2 + (1 - x(1:end-1)).^2);
end