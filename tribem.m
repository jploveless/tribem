function [slip, trac, o, G] = tribem(patch, d, bc, varargin)
% TRIBEM  Triangular element boundary element solver using geographic coordinates.
%   [SLIP, TRAC] = TRIBEM(PATCH, D, BC) estimates slip or traction on triangular
%   elements given a mesh in structure PATCH, a set of slip and traction
%   boundary conditions D and a boundary condition identifer matrix, BC.
%
%   PATCH can either be a string giving the full patch to a Gmsh .msh file 
%   or a structure given as PATCH = READPATCHES(PATCHFILE).  Given NEL 
%   triangular elements, D should be an NEL-by-3 array of boundary conditions, 
%   giving either the slip or traction conditions on each element.  For elements
%   with slip prescribed, the corresponding line of D should be:
%
%   [strike slip, dip slip, tensile slip]
% 
%   and for elements whose traction conditions are specified, D should read:
%
%   [shear stress (strike), shear stress (dip), normal stress].
%
%   The NEL-by-3 array BC specifies the type of boundary conditions for each 
%   element, given as 0 for slip conditions and 1 for traction conditions.
%
%   [SLIP, TRAC, O] = TRIBEM(PATCH, D, BC, OBS) will also calculate values at a set of 
%   observation coordinates defined by input argument OBS.  OBS can be a structure
%   containing fields lon, lat, dep, and v, where lon, lat, and dep are N-by-1 arrays 
%   giving the coordinates of N observation points (dep in km), and v is an integer 
%   specifying what will be calculated at these points:
%      v = 1: Calculate displacement only (returned to O.u as a 3N-by-1 vector)
%      v = 2: Calculate stress only (returned to O.s as a 6N-by-1 vector)
%      v = 3: Calculate displacement and stress
%
%   [SLIP, TRAC, G] = TRIBEM(PATCH, D, BC) will output the structure G containing 
%   fields of Green's functions relating slip to displacement and/or traction. This 
%   structure can be reused as an input argument, speeding computation when only the
%   boundary conditions and not the problem geometry have changed. 
%
%   [SLIP, TRAC, O, G] = TRIBEM(PATCH, D, BC, OBS) will output both O and G, when OBS
%   is specified as an input argument. 
%
%   [...] = TRIBEM(PATCH, D, BC, REMS) allows specification of a remote stress tensor,
%   REMS. The tensor components, assuming Cartesian coordinates, should be given as
% 
%   REMS = [sigma_xx, sigma_yy, 0, sigma_xy, 0, 0]'; % 0 values for components involving z
%
%   [...] = TRIBEM(PATCH, D, BC, G) uses a structure G of already-calculated
%   traction (G.sp) and displacement partials (G.u).  No checking is done to assure
%   that the sizes of G fields are appropriate for the problem; this must be done
%   externally.  Therefore, this option is best used when only the boundary conditions
%   change from one model run to the next, i.e., all of the problem geometry remains
%   consistent.
%
%   [...] = TRIBEM(PATCH, D, BC, CHOP) accepts CHOP as a flag to indicate whether
%   or not the components of slip should be chopped according to element geometry.  If
%   CHOP > 0,  only strike-slip and dip-slip will be considered for dipping faults, and
%   only strike-slip and tensile motion will be considered for vertical faults.
%
%   [...] = TRIBEM(PATCH, D, BC, G, CHOP) accepts both G and CHOP as input arguments.
%   If the input G fields are already chopped, no additional trimming will be done.
%
%
% SIGN CONVENTIONS:
%------------------
%
% SLIP:
%                +      |     - 
% -------------------------------
% 1.  Strike: Sinistral | Dextral
% 2.     Dip: Reverse   | Normal
% 3. Tensile: Opening   | Closing
%
% TRACTION:
%                     +    |     - 
% ----------------------------------------
% 1. Strike shear: Strike  | Strike + 180
% 2.    Dip shear: Updip   | Downdip
% 3.       Normal: Tensile | Compressional
%
% For comparison, Poly3D's sign conventions are:
%
% SLIP:
%                +      |     - 
% -------------------------------
% 1.     Dip: Normal    | Reverse
% 2.  Strike: Sinistral | Dextral
% 3. Tensile: Opening   | Closing
%
% TRACTION:
%                      +   |     - 
% ----------------------------------------
% 1.    Dip shear: Downdip | Updip
% 2. Strike shear: Strike  | Strike + 180 
% 3.       Normal: Tensile | Compressional
%

chop = 0;

% Parse optional inputs
noa = length(varargin);
if noa > 0
   for i = 1:noa % Loop through optional arguments
      if isstruct(varargin{i}) % If it's a structure, it's either observations or partials
         if isfield(varargin{i}, 'lon') % If field lon exists, it's observations
            obs = varargin{i};
         else % If not, 
            G = varargin{i}; % It's partials
         end
      else
         if length(varargin{i}) == 1 % If a single value was entered, it's either chop flag or Poisson's ratio
            if varargin{i} <= 0.5 % If it's less than or equal to 0.5, it's PR
               pr = varargin{i};
            else
               chop = varargin{i}; % If not, it's chop flag
            end
         elseif length(varargin{i}) == 2 % If a 2-element vector was specified,
            lame = varargin{i}; % Lame parameters were specified            
         elseif length(varargin{i}) == 6
            rems = varargin{i}; % Remote stress tensor
         end
      end
   end
end

% Read patch file, if necessary
if ischar(patch)
   patch = ReadPatches(patch);
   patch = PatchCoords(patch);
end

% Calculate patch coordinates, if necessary
if ~isfield(patch, 'lonc')
   patch = PatchCoords(patch);
end

% Total number of elements
tne = sum(patch.nEl);

% Project remote stress tensor onto elements
if exist('rems', 'var')
   rems([3 5 6]) = 0; % Make sure traction-free surface of half space condition is met
   remsp = ProjectStrainPartialsMats(repmat(rems(:), sum(patch.nEl), 1), patch.strike, patch.dip);
end
   
% Make a structure containing perturbed element centroid coordinates
% Centroids are perturbed toward the half-space surface along the normal vector
psc = 1e3; % Centroid perturbation scaling factor
[cent.lon, cent.lat, cent.dep] = deal(patch.lonc+psc.*patch.nv(:, 1), patch.latc+psc.*patch.nv(:, 2), patch.zc+psc.*patch.nv(:, 3));

% Augment this structure with any observation coordinates
if exist('obs', 'var')
   obs.lon = obs.lon(:);
   obs.lat = obs.lat(:);
   obs.dep = obs.dep(:);
   if ~isfield(obs, 'v') % Default behavior is to calculate displacements only at observation coordinates
      obs.v = 1;
   end
   if obs.v == 1 % Calculate displacements at observation coordinates
      opt = repmat([1 0], length(obs.lon), 1);
   elseif obs.v == 2 % Calculate stresses at observation coordinates
      opt = repmat([0 1], length(obs.lon), 1);
   elseif obs.v == 3 % Calculate displacements and stresses at observation coordinates
      opt = repmat([1 1], length(obs.lon), 1);
   else
      opt = zeros(0, 2);
   end
else
   [obs.lon, obs.lat, obs.dep] = deal([]);
   obs.v = 999;
   opt = zeros(0, 2);
end

% Augmented calculation coordinates
cc.lon = [cent.lon; obs.lon];
cc.lat = [cent.lat; obs.lat];
cc.dep = [cent.dep; obs.dep];
% Set partials calculation option; this will be adjusted if need be when observations are considered
opt = logical([repmat([0 1], tne, 1); opt]); 

% Call GetTriCombinedPartials for elements and coordinates
if ~exist('G', 'var')
   % Check for existence of specified Lame parameters
   if ~exist('lame', 'var')
      mu = 3e10; lambda = 3e10; % Defaults
   else
      mu = lame(1); lambda = lame(2);
   end
   % Check for existence of specified Poisson's ratio
   if ~exist('pr', 'var')
      pr = 0.25;
   end
   [G.u, G.e, G.tz] = GetTriCombinedPartials(patch, cc, opt, pr);
   G.e = G.e./1e6; % Internal conversion assuming mm is unit of slip
   G.e = StrainToStressComp(G.e', mu, lambda)'; 
   G.sp = ProjectStrainPartialsMats(G.e(1:6*tne, :), patch.strike, patch.dip);
end

% Trim arrays for slip components, if requested and if not already done
if chop > 0
   triD = find(G.tz == 2);
   triT = find(G.tz == 3);
   colkeep = setdiff(1:3*tne, [3*triD-0; 3*triT-1]); % Define columns to keep
   if size(G.sp, 2) == numel(patch.v) % If the columns have not already been trimmed...
      G.e = G.e(:, colkeep); % ...trim them
      G.sp = G.sp(:, colkeep);
      if size(G.u, 2) > size(G.e, 2)
         G.u = G.u(:, colkeep);
      end   
   end
else
   colkeep = 1:3*tne; % All columns retained if no chopping is requested
end

% Determine which elements need slip estimated, and which have slip prescribed
c3 = 0; % Flag to indicate that the arrays came in as 3 column arrays 
if max([size(bc, 2), size(d, 2)]) == 3
   c3 = 1; 
   d = reshape(d', numel(d), 1);
   bc = reshape(bc', numel(bc), 1);
end
bc = logical(bc(:));
islip = find(~bc(colkeep)); % Indices of prescribed slips
pslip = d(colkeep); pslip = pslip(islip); % Actual prescribed slips
itrac = find(bc); % Indices of prescibed tractions
ptrac = d(itrac); % Actual prescribed tractions
esidx = find(bc(colkeep)); % Indices of the slip components that will be estimated

% Determine traction induced by any input displacement boundary conditions
dstress = G.sp(:, islip)*pslip;

if ~exist('remsp', 'var') % If there is no remote stress
   % Corrected traction is prescribed traction minus contribution from displacement b.c.
   elstress = ptrac - dstress(itrac);
   % Estimate slips that produce corrected traction
   m = G.sp(itrac, esidx)\elstress;
   % Assemble full slip vector (prescribed and estimated)
   slip = zeros(3*tne, 1);
   slip(islip) = pslip;
   slip(itrac) = m;
   % Assemble full traction vector (prescribed and estimated)
   trac = zeros(3*tne, 1);
   trac(islip) = dstress(islip);
   trac(itrac) = ptrac;
else
   elstress = 0*d; % Full element stress vector
   elstress(itrac) = ptrac; % Corrected traction is prescribed traction minus contribution from displacement b.c. 
   elstress = elstress - dstress - remsp; % Also corrected for remote stress projected onto elements
   g = G.sp; seye = eye(size(G.sp, 1));
   g(:, islip) = seye(:, islip);
   m = g\elstress;
   % Assemble full slip vector (prescribed and estimated)
   slip = m;
   slip(islip) = pslip;
   % Assemble full traction vector (prescribed and estimated)
   trac = elstress;
end

% Forward calculation for observations
if obs.v == 1 | obs.v == 3
   o.u = G.u*slip(:);
else
   o.u = [];
end
if obs.v == 2 | obs.v == 3
   o.s = G.e(6*tne+1:end, :)*slip(:);
else
   o.s = [];
end

% Reshape if original boundary condition array was n-by-3
if c3 == 1
   slip = reshape(slip, 3, tne)';
   trac = reshape(trac, 3, tne)';
end

% Process optional output arguments
if nargout == 3 % 3rd argument could be o or G
   if obs.v == 999 % If no observation points are specified, it must be G
      varargout{:} = G;
   else
      varargout{:} = o;
   end
end

if nargout == 4
   varargout{1} = o;
   varargout{2} = G;
end
