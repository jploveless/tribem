function ec = PatchEdges(p)
%  PatchEdges  Determines coordinates of fault edges in patch structure.
%    EC = PatchEdges(P) determines the 3-D coordinates of the outlines of
%    faults as specified in structure P. The output EC is 3 column array
%    containing [x, y, z] (or [lon, lat, z]) coordinates of the edges, with
%    [NaN NaN NaN] rows separating edges of distinct entities so that all
%    edges can be plotted with a simple plot or plot3 command.
%

ends = cumsum(p.nEl);
begs = [1; ends(1:end-1)+1];
cc = [0; cumsum(p.nc)];
ec = [];

for i = 1:length(p.nEl) % For each entity, 
   oe = OrderedEdges(p.c, p.v(begs(i):ends(i), :));
   ec = [ec; p.c([oe(1, :), oe(end)], :); nan(1, 3)];
end