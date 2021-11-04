function [p, d, bc, rems, obs] = poly3d2tribem(file)



% Read .in file as text
a = opentxt(file);

% Find indices of objects, nodes, and elements
objects = findstr(a(:, 1)', 'o');
nodes = findstr(a(:, 1)', 'v');
els1 = findstr(a(:, 1)', 'e');
els2 = findstr(a(:, 2)', ' ');
els = intersect(els1, els2);

% Index differences
dnode = diff(nodes);
dels = diff(els);

% Beginning and end of node, element series
nbegs = [nodes(1), nodes(find(dnode > 1) + 1)];
nends = [nodes(find(dnode > 1)), nodes(end)];

ebegs = [els(1), els(find(dels > 1) + 1)];
eends = [els(find(dels > 1)), els(end)];

% Define number of elements 
no = size(objects, 2);
p.nc = (nends - nbegs + 1)';
p.nEl = (eends - ebegs +1)';

% Gather nodes, elements, and boundary conditions
p.c = zeros(sum(p.nc), 3);
p.v = zeros(sum(p.nEl), 3);
d = p.v;
bc = p.v;

elcount = 1;
for i = 1:length(nbegs)
   for j = nbegs(i):nends(i)
      nline = textscan(a(j, :), 'v %f global %f %f %f\r');
      p.c(nline{1}, :) = [nline{2}, nline{3}, nline{4}];
   end
   for j = ebegs(i):eends(i)
      eline = textscan(a(j, :), 'e 3 elocal %s %f %f %f %f %f %f\r');
      p.v(elcount, :) = [eline{5}, eline{6}, eline{7}];
      bctext = char(eline{1}); bctext1 = bctext; bctext(1) = bctext1(2); bctext(2) = bctext1(1);
      bct = strfind(bctext, 't');
      bc(elcount, bct) = 1;
      d(elcount, :) = [eline{3}, -eline{2}, eline{4}]; % Reverse order of strike and dip directions
      elcount = elcount + 1;
   end
end