function [ totalelements ] = rake(totalelements, slip)  %ppositive
%calculates rake on each element of the mesh, given mesh and slip inputs
totalelements.sliprake=totalelements.dip*0;
totalelements.sliprake=atan2d(-slip(:, 2), slip(:, 1));


