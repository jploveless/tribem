function [C, angrange, Rbar, G] = triregstress(ds, patch, varargin)
% triregstress   Evaluate fit of regional stress to stress data
%   [C, shaz, Rbar] = triregstress(ds, patch) evaluates the fit of trial
%   regional stress tensors to stress indicator data in the presence of
%   triangular dislocation elements that represent faults. The input data
%   represent stress indicator data (fractures) and are specified in 
%   structure ds, containing fields: 
%      x: x-coordinate of stress observations
%      y: y-coordinate
%      z: z-coordinate
%      az: azimuth of fracture
%      ft: fracture type (1 = closing, 2 = shear, 3 = opening)
%   Input argument patch is a structure containing fault geometry data, 
%   as read from tridisl/ReadPatches.m. 
%
%   A suite of stress tensors is generated, with a ratio of principal 
%   stresses from 0-3, after Maerten et al., JSG, 2016. For each ratio, 
%   maximum horizontal stress azimuths between 0-180 degrees are tested. 
%   The total stress calculated at the observation coordinates is the sum
%   of the remote stress and perturbations to the stress field arising from
%   fault slip. Fault slip is assumed to take place to completely relieve
%   the shear traction imposed by regional stress on the fault elements.
%
%   The misfit of each stress field trial is returned to the cost matrix, C.
%   The range of maximum horizontal stress azimuths is returned to vector 
%   shaz, and the range of principal stress ratios is returned to vector 
%   Rbar. A plot of cost as a function of horizontal stress azimuth and 
%   principal stress ratio is generated, with the minimum misfit identified
%   at the intersection of the two gray lines.
%   



% Parse optional inputs
noa = length(varargin);
if noa > 0
   for i = 1:noa % Loop through optional arguments
      if isstruct(varargin{i}) % If it's a structure, it's partials
         G = varargin{i}; % It's partials
      else
         if length(varargin{i}) == 1 % If a single value was entered, it's either chop flag or Poisson's ratio
            if varargin{i} <= 0.5 % If it's less than or equal to 0.5, it's PR
               pr = varargin{i};
            end
         elseif length(varargin{i}) == 2 % If a 2-element vector was specified,
            lame = varargin{i}; % Lame parameters were specified            
         end
      end
   end
end

% Read patch file, if necessary
if ischar(patch)
   patch = ReadPatches(patch);
   patch = PatchCoordsx(patch);
end

% Calculate patch coordinates, if necessary
if ~isfield(patch, 'xc')
   patch = PatchCoordsx(patch);
end

% Total number of elements
tne = sum(patch.nEl);

% Number of data points
nds = numel(ds.x);

% Evaluate stress indicator type and calculate normals
if ~isfield(ds, 'nv')
   [dx, dy]= sph2cart(deg2rad(ds.az + 90), 0, 1);
   ds.nv = [dy, dx, 0*dx];
end

% Check for existence of specified Lame parameters
if ~exist('lame', 'var')
   mu = 3e10; lambda = 3e10; % Defaults
else
   mu = lame(2); lambda = lame(1);
end

% Make a structure containing perturbed element centroid coordinates
% Centroids are perturbed toward the half-space surface along the normal vector
psc = 1e-5; % Centroid perturbation scaling factor
[cent.x, cent.y, cent.z] = deal(patch.xc+psc.*patch.nv(:, 1), patch.yc+psc.*patch.nv(:, 2), patch.zc+psc.*patch.nv(:, 3));

% Augmented calculation coordinates (element centroids and observation coordinates)
cc.x = [cent.x; ds.x(:)];
cc.y = [cent.y; ds.y(:)];
cc.z = [cent.z; ds.z(:)];

% Define combined stress-slip-traction-remote stress partial derivatives
if ~exist('G', 'var')
   % Check for existence of specified Poisson's ratio
   if ~exist('pr', 'var')
      pr = lambda./(2*(lambda + mu));
   end
   % Full matrices of partial derivatives
   [~, G.e, G.tz] = GetTriCombinedPartialsx(patch, cc, [0 1], pr);
   % Convert strain partials to stress partials
   G.s = StrainToStressComp(G.e', mu, lambda)';
   % Project stress partials to traction partials
   G.sp = ProjectStrainPartialsMats(G.s(1:6*tne, :), patch.strike, patch.dip); % Isolating rows corresponding to elements
   g = G.sp; seye = eye(size(G.sp, 1)); % Set up design matrix
   g(:, 3:3:end) = seye(:, 3:3:end); % Substitute identity matrix for tensile slip columns
   % Extract stress partials for data points
   G.sd = G.s(6*tne+1:end, :); 
   % Calculate regional stress projection partials
   G.rp = GetProjectedRegionalTensorPartials(patch.strike, patch.dip);

   % Combine observation and slip stress partials with regional tensor projection partials
   G.spi = g\eye(size(g, 1)); % Invert traction partials
   G.spi(3:3:end, :) = 0; % Zero out tensile slip components
   G.tot = G.spi*-G.rp; % Traction partials times projection partials
   G.tot = G.sd*G.tot; % Observation stress partials times combined
end
%
% Define stress field test ranges
%
angrange = 0:1:179; % Azimuth of SHmax
naz = length(angrange);

sinc = 0.5; % Principal stress magnitude increment
[s1, s2, s3] = meshgrid(-1, -1:sinc:1, -1:sinc:1); % Nominal principal stress value ranges
sall = [s1(:), s2(:), s3(:)]; % Vectorize normal stress values
sall = sort(sall, 2); % Sort principal stress values
keep = sall(:, 1) ~= sall(:, 2) | sall(:, 2) ~= sall(:, 3); % Keep only those with non-all-equal stresses
sall = sall(keep, :); 
R = (sall(:, 2) - sall(:, 1))./(sall(:, 3) - sall(:, 1)); % Calculate principal stress ratios
[R, ui] = unique(R); % Find unique values of R
sall = sall(ui, :); % Retain smallest set of unique principal stresses
nstr = size(sall, 1);

%
% Begin stress field trials
%
cols = [2 3 1; 1 3 2; 1 2 3]; % Column indices to map principal stresses to Cartesian axes
% Allocate space for arrays
c = zeros(nds, 1);
C = zeros(naz, 3*nstr);
Rbar = zeros(3*nstr, 1);

% For each stress domain (normal, wrench, reverse),
for h = 1:3
   % For each set of principal stress magnitudes,
   for i = 1:nstr 
      ns = diag(sall(i, cols(h, :))); % Create principal stress matrix
      Rbar(nstr*(h-1)+i) = R(i) + (h-1); % Place Rbar value
      % For each horizontal stress azimuth,
      for j = 1:naz 
         % Define azimuthal rotation matrix; these are principal stress eigenvectors
         rot = [sind(angrange(j)), cosd(angrange(j)), 0; cosd(angrange(j)), -sind(angrange(j)), 0; 0 0 1];
         % Rotate stress components
         smat = rot'*ns*rot;
   
         % Multiply remote stress by fault-related partials to give stress at observation coordinates
         % Add remote stress to give total stress at each point
         stot = G.tot*smat([1 5 9 2 3 6])' + repmat(smat([1 5 9 2 3 6])', nds, 1);
         stotm = stressrow2mat(unstack6(stot)); % Create stacked 3x3 tensors for each observation 
         % For each observation,
         for k = 1:nds
            % Calculate principal stresses orientations
            [psv, ~] = eig(stotm(3*k-2:3*k, :), 'vector');
            % Calculate misfit: 1 minus dot product of sigma_3 and element normal
            c(k) = 1 - dot(psv(:, ds.ft(k)), ds.nv(k, :)').^2;
         end % End observation loop
         % Cost for this stress field trial is mean of observation costs
         C(j, (h-1)*nstr + i) = mean(c);
      end % End azimiuth loop
   end % End stress magnitude loop
end % End stress domain loop

% Create cost figure
figure
pcolor(Rbar, angrange, C); shading flat
xlabel('Stress ratio (R)');
ylabel('\sigma_H orientation (\theta)');
[minaz, minr] = find(C == min(C(:)));
hold on
aa = axis;
line([1 2; 1 2], [aa(3:4); aa(3:4)], 'color', 'k', 'linewidth', 1.5);
line([Rbar(minr); Rbar(minr)], aa(3:4)', 'color', 0.5*[1 1 1], 'linewidth', 1);
line(aa(1:2)', [angrange(minaz); angrange(minaz)], 'color', 0.5*[1 1 1], 'linewidth', 1);
c = colorbar;
c.Label.String = 'Cost';
fs = 12; % font size
lw = 1; % line width
tl = [0.008 0.008]; % tick length
set(gca, 'fontsize', fs, 'linewidth', lw, 'ticklength', tl, 'layer', 'top', 'box', 'on');