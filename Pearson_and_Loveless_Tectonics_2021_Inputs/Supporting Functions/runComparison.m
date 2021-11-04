function [rakes, results, anglesbetween, misfit] = runComparison(faults, faultsHotspot, hotspotContribution, azimuth, faultOpening, comparison, slipData, stressData, eqData, eqDataStress, GPSData)
%prepare inputs for tribemx
d=zeros(sum(faultsHotspot.nEl), 3);
d(sum(faults.nEl)+1:end, 3)=hotspotContribution;
bc=ones(sum(faultsHotspot.nEl), 3);
bc(:, 3)=0;

if strcmp(faultOpening, 'None')==1
    %standard, do not allow fault opening
    bc(:, 3)=0;
elseif strcmp(faultOpening,'All')==1
    %let all faults open
    bc(1:sum(faults.nEl), 3)=1;
else
    %let poorly fit faults open
    for i=1:length(faultOpening)
        bc(begs(i):ends(i), 3)=1;
    end
end


[REMSrotate, ~]=calculateREMS(azimuth);
if comparison==0 
    [slip, trac, G] = tribemx(faultsHotspot, d, bc, REMSrotate);
    faultsHotspot=rake(faultsHotspot, slip);
    rakes=faultsHotspot.sliprake;
end




if comparison==0  %fault-slip derived, moment tensor axes
    [~, ~, misfit, anglebetweenP, anglebetweenT, pressure, tension]=momentTensorsFaults(faults, faultsHotspot, d, bc, azimuth, slipData(:, 1:2), slipData(:, 3:4), G);
    anglesbetween=(anglebetweenP(:, 1)+anglebetweenT(:,1))/2;
    anglesbetween=nonzeros(unique(anglesbetween, 'rows', 'stable'));
    results=[pressure, tension];
    misfit=mean(mean(anglesbetween));
elseif comparison==1 %fault-slip derived, stress tensor axes
    [misfit, anglesbetween, stressAxes, stressRatio, flips, ~, rakes]=stressTensorsFaults(faults, faultsHotspot, d, bc, azimuth,stressData);
    results=stressAxes;
    anglesbetween=unique(anglesbetween, 'rows', 'stable');
    anglesbetween = anglesbetween(any(anglesbetween,2),:);
elseif comparison==2  %earthquake derived, moment tensor axes
    [momAxes, misfit, anglesbetween, rakes] = momentTensorsCoordinates(faults, faultsHotspot, d, bc, azimuth,eqData);
    results=momAxes;
elseif comparison==3  %earthquake-derived, stress tensor axes
    [stressAxes, anglesbetween, misfit, stressRatio, flips, sigma, rakes] = stressTensorsCoordinates(faults, faultsHotspot, d, bc, azimuth,eqDataStress);
    results=stressAxes;
elseif comparison==4 %GPS axes
    [anglesbetween, misfit, results, rakes] = GPSCoordinates(faults, faultsHotspot, d, bc, azimuth,hotspotContribution, GPSData);
end

end