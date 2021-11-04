function [faults, faultsHotspot] = meshFaults(filename, hotspot, dips, path)
%meshFaults takes in fault mesh and spherical hotspot mesh and outputs
%patch for running in tribemx
faults=load(filename);
coords=utmconv(faults(:,[1 2]), 0);
faults=[coords, faults(:,3)];
save ('faultsutm.txt', 'faults', '-ascii')
maketopbotcoords('faultsutm.txt', 10e3);
rotatefaultcoords_dir('faultsutm_botcoords', dips);
addpath(strcat(erase(filename, '.txt'), 'utm_botcoords'), strcat(erase(filename, '.txt'), 'utm_botcoords_truedip'), strcat(erase(filename, '.txt'), 'utm_topcoords'));
faults=gmshfaults_Windows('faultsutm_topcoords', 'faultsutm_botcoords_truedip', 7500, path);
faults = ReadPatches(faults);
faults = PatchCoordsx(faults);
faultsHotspot=patchcat(faults, hotspot);
end