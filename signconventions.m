function signconventions()
% signconventions   Tests relationships in sighs of traction, slip, and stress/strain.
% 

% Establish fault geometry
[bx, bz] = meshgrid(0:2.5:50, 0:2.5:25);
b.c = [bx(:), -bz(:)*cosd(60), -bz(:)]; % 60 degree south dip
b.nc = length(b.c);
b.v = delaunay(bx(:), bz(:)); 
b.v = fliplr(b.v); % Enforce counterclockwise circulation, upward normal vector
b.nEl = length(b.v);

% Data and boundary conditions arrays (zero traction everywhere)
d = repmat([0 0 0], size(b.v, 1), 1); % All zero values
bc = repmat([1 1 0], size(b.v, 1), 1); % All traction conditions

% Observation grid
[x, y] = meshgrid(-10.1:1:60, -15.1:1:15);
obs.x = x(:); obs.y = y(:); obs.z = 0*ones(size(obs.x)); 
% Calculate stress, displacement, and strain
obs.v = 'sde';
keyboard
% Run tribem with remote stress tensor of 1 unit in the s_yy direction
% (Tension perpendicular to fault)
[slip, trac, out, G] = tribemx(b, d, bc, [0 1 0 0 0 0], obs);

% Plot slips
meshview(b.c, b.v, slip(:, 1)); title('Strike slip (sinistral +)')
meshview(b.c, b.v, slip(:, 2)); title('Dip slip (reverse +)')
meshview(b.c, b.v, slip(:, 3)); title('Normal slip (opening +)')

% Plot tractions
meshview(b.c, b.v, trac(:, 1)); title('Strike traction (sinistral +)')
meshview(b.c, b.v, trac(:, 2)); title('Dip traction (reverse +)')
meshview(b.c, b.v, trac(:, 3)); title('Normal traction (tension +)')

% Plot displacements as vectors
figure
quiver3(obs.x, obs.y, obs.z, out.u(:, 1), out.u(:, 2), out.u(:, 3));
title('Displacements due to fault slip')

% Plot stress components as pcolors
figure; pcolor(x, y, reshape(out.s(:, 1), size(x))); colormap(bluewhitered); title(['\sigma_{xx}'])
figure; pcolor(x, y, reshape(out.s(:, 2), size(x))); colormap(bluewhitered); title(['\sigma_{yy}'])
figure; pcolor(x, y, reshape(out.s(:, 3), size(x))); colormap(bluewhitered); title(['\sigma_{zz}'])
figure; pcolor(x, y, reshape(out.s(:, 4), size(x))); colormap(bluewhitered); title(['\sigma_{xy}'])
figure; pcolor(x, y, reshape(out.s(:, 5), size(x))); colormap(bluewhitered); title(['\sigma_{xz}'])
figure; pcolor(x, y, reshape(out.s(:, 6), size(x))); colormap(bluewhitered); title(['\sigma_{yz}'])

% Check stress component signs by showing displacement gradients
figure; pcolor(x(:, 1:end-1), y(:, 1:end-1), diff(reshape(out.u(:, 1), size(x)), 1, 2)); colormap(bluewhitered); title(['\delta{u}_{x}/\delta{x}'])
figure; pcolor(x(1:end-1, :), y(1:end-1, :), diff(reshape(out.u(:, 2), size(x)))); colormap(bluewhitered); title(['\delta{u}_{y}/\delta{y}'])
dxy = diff(reshape(out.u(:, 1), size(x)));
dyx = diff(reshape(out.u(:, 2), size(x)), 1, 2);
exy = 0.5*(dxy(:, 1:end-1) + dyx(1:end-1, :));
figure; pcolor(x(1:end-1, 1:end-1), y(1:end-1, 1:end-1), exy); colormap(bluewhitered); title(['\delta{u}_{x}/\delta{y}'])