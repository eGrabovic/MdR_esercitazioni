%% creazione punti scritta
% read image and convert to grayscale matrix
I = rgb2gray(imread('Scritta.png'));
% h = ones(3,3)/10;
% I = imfilter(I,h); % filter image for smoother contour

% extract contour from image
[C,h] = imcontour(flipud(I), 1);
X = C(1,:);
Y = C(2,:);

figure; hold on; axis equal;
plot(X, Y);

% remove contiguous points that are too far apart from each other (contour noise)
[X, Y, k] = remove_noise(X, Y);


% center image around center of gravity
xm = mean(X); ym = mean(Y);
X = X-xm;
Y = Y-ym;

% close the loop (last point = initial point)
X(end+1) = X(1);
Y(end+1) = Y(1);

% interpolation and equispacing of points 
P =  interparc(2000, X, Y, 'linear').';
X = P(1,2:end);
Y = P(2,2:end);
figure; hold on; axis equal;
plot(X, Y);

% save coordinates to avoid recalculating 
save('coordinateX_scritta.mat', 'X')
save('coordinateY_scritta.mat', 'Y')
