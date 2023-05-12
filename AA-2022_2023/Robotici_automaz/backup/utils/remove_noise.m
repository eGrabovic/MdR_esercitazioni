function [X, Y, k] = remove_noise(X, Y)

d = sqrt(diff(X).^2 + diff(Y).^2); % contiguous points distance
k = find(d < 5*mean(d)); % find indices of points sufficiently close to each other
X = X(k);                % replace X and Y arrays with said points
Y = Y(k);
d = sqrt(diff(X).^2 + diff(Y).^2);
k = find(d > 5*mean(d));

for j = 1:length(k)-1
    if k(j+1) - k(j) < 5
        X(k(j)+1:k(j+1)) = [];
        Y(k(j)+1:k(j+1)) = [];
    end
end

% repeat process to remove any possible outliers
d = sqrt(diff(X).^2 + diff(Y).^2); % contiguous points distance
k = find(d < 5*mean(d)); % find indices of points sufficiently close to each other
X = X(k);                % replace X and Y arrays with said points
Y = Y(k);
d = sqrt(diff(X).^2 + diff(Y).^2);
k = find(d > 5*mean(d));

end