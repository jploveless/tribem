% tribem simple remote stress example

% Establish fault geometry
[bx, bz] = meshgrid(0:2.5:50, 0:2.5:25);
b.c = [bx(:), zeros(size(bx(:))), -bz(:)];
b.nc = length(b.c);
b.v = delaunay(bx(:), bz(:));
b.nEl = length(b.v);

% Data and boundary conditions arrays (zero traction everywhere)
d = zeros(size(b.v)); % All zero values
bc = ones(size(b.v)); % All traction conditions

% Observation grid
[x, y] = meshgrid(-10.1:1:60, -15.1:1:15);
obs.x = x(:); obs.y = y(:); obs.z = -15*ones(size(obs.x)); 
% Calculate stress, displacement, and strain
obs.v = 'sde';

% Run tribem with remote stress tensor of 490 units in the s_yy direction
% (Tension perpendicular to fault)
[slip, trac, out, G] = tribemx(b, d, bc, [0 490 0 0 0 0], obs);

% Visualize traction and slip in the element-normal direction
meshview(b.c, b.v, trac(:, 3)); 
view([0 0]); axis square; axis equal; title('Element-normal traction')
meshview(b.c, b.v, slip(:, 3)); 
view([0 0]); axis square; axis equal; title('Element-normal slip')
ca = caxis;


% Run tribem with remote stress tensor of 490 units in the s_yy direction
% and lithostatic loading (density of 2, corresponding to 490 units of 
% compression at bottom of fault (25 km depth), canceling out remote stress)
% (Tension perpendicular to fault)
[lslip, ltrac, lout] = tribemx(b, d, bc, [0 490 0 0 0 0], obs, G, 2);
meshview(b.c, b.v, ltrac(:, 3)); 
view([0 0]); axis square; axis equal; title('Element-normal traction, lithostatic')
meshview(b.c, b.v, lslip(:, 3));  
view([0 0]); axis square; axis equal; title('Element-normal slip, lithostatic')
caxis(ca); colormap(bluewhitered); % Set color scale equal to other slip figure
