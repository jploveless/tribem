function [ fitBoth, fit1, degreeOfMisfit, anglebetweenPnew, anglebetweenTnew, pressure, tension] = momentTensorsFaults( faults, faultsHotspot, d, bc, alpha, toCompareP, toCompareT, G, varargin)
%momentTensorsFaults solves for angular misfit with comparison moment tensor axes
%defined at a series of faults

%Inputs: 
% faults: meshed fault map       faultsHotspot: meshed fault map and spherical mesh hotspot
% d: nelx3 array of boundary conditions & bc: NEL-by-3 array specifying type
% of boundary conditions, see tribemx for specifications 
% alpha: angle of uniaxial remote stress, defined clockwise from due east
% G: pre-calculated Green's functions, derived from tribemx
% toCompareP/toCompareT: pressure and tension axes, respectively, from
% comparison data
% varargin: 

%Outputs: 
%fitBoth/fit1: the number of regions in which both moment tensor axes are fit well/one is fit well
%degreeOfMisfit: average overall angular misfit
%anglebetweenPnew/anglebetweenTnew: angular misfits between pressure and tension axes for each
%comparison data point
%pressure/tension: model output of pressure and tension axes


%run slip, rake, and slip magnitude calculations
[REMSrotate, ~]=calculateREMS(alpha);
numvarargs=length(varargin);  %slip is pre-calculated when running for a gradient of remote stress
if numvarargs==1
    slip=varargin;
else
    [slip, ~] = tribemx(faultsHotspot, d, bc, REMSrotate, G);
end
faultsHotspot=rake(faultsHotspot, slip);
rakes=faultsHotspot.sliprake(1:size(faults.v));
dips=faultsHotspot.dip(1:size(faults.v));
strikes=faultsHotspot.strike(1:size(faults.v));
strikesAndDips=[deg2rad(strikes), deg2rad(dips)];
slipmag=zeros(sum(faults.nEl), 1);
for i=1:sum(faults.nEl)
    slipmag(i)=sqrt(slip(i, 1)^2+slip(i, 2)^2+slip(i, 3)^2);
end


%convert rakes on each element to plunges and trends
[trends, plunges]=raketoTrendPlunge(strikes, dips, rakes);

%run PTAxes to determine trend and plunge of 
%pressure and tension for each element, then adjust for multiple faults in
%a single region
ends=cumsum(faults.nEl);
begs=[1; ends(1:length(ends)-1)+1];
pressures=zeros(length(ends), 2);
tensions=zeros(length(ends), 2);

for i=1:length(ends)
    [pressures(i, :), tensions(i, :)]=MomTens(strikesAndDips(begs(i):ends(i), 1:2), [deg2rad(trends(begs(i):ends(i))), deg2rad(plunges(begs(i):ends(i)))], slipmag(begs(i):ends(i)));
end
%adjust for regions with multiple faults
numbers=[25, 26, 0, 0; 27, 28, 0, 0; 29, 30, 0, 0; 5, 34, 35, 37]; %defined numbers of faults where multiple are in a single region 
for i=1:size(numbers, 1)
    [pressures, tensions]=multipleFaults(pressures, tensions, begs, ends, strikesAndDips, trends, plunges, slipmag, nonzeros(numbers(i, :)).');
end


%adjust formatting for plunges (must be positive in order to come out
%correctly), trends must then be switched by 180' to remain consistent
%with which side of plunge is being considered
for i=1:length(pressures)
    if pressures(i, 2)<0
        if pressures(i, 1)<pi
            pressures (i, 1)=pressures(i, 1)+pi;
        else
            pressures(i, 1)=pressures(i, 1)-pi;
        end
        pressures(i, 2)=-pressures(i, 2);
    end
end
 for i=1:length(pressures)
     if tensions(i, 2)<0
         if tensions(i, 1)<pi
              tensions (i, 1)=tensions(i, 1)+pi;
          else
              tensions(i, 1)=tensions(i, 1)-pi;
          end
         tensions(i, 2)=-tensions(i, 2);
     end
 end

pressure=rad2deg(pressures);
tension=rad2deg(tensions);
%--------------------------------------------------------------------------
%calculate angles between calculated P and T axes and the axes chosen for
%comparison
%fix trend angle to be in correct format (0'=due east)
%for conversion to cartesian coordinates
pressureTrendvector=mod(pressure(:, 1)+90, 360);
for i=1:length(pressureTrendvector)
    if pressureTrendvector(i, 1)<180
        pressureTrendvector(i, 1)=180-pressureTrendvector(i, 1);
    else
        pressureTrendvector(i, 1)=540-pressureTrendvector(i, 1);
    end
end

tensionTrendvector=mod(tension(:, 1)+90, 360);
for i=1:length(tensionTrendvector)
    if tensionTrendvector(i, 1)<180
        tensionTrendvector(i, 1)=180-tensionTrendvector(i, 1);
    else
        tensionTrendvector(i, 1)=540-tensionTrendvector(i, 1);
    end
end

toComparePTrend=mod(toCompareP(:, 1)+90, 360);
for i=1:length(pressureTrendvector)
    if toComparePTrend(i, 1)<180
        toComparePTrend(i, 1)=180-toComparePTrend(i, 1);
    else
        toComparePTrend(i, 1)=540-toComparePTrend(i, 1);
    end
end

toCompareTTrend=mod(toCompareT(:, 1)+90, 360);
for i=1:length(pressureTrendvector)
    if toCompareTTrend(i, 1)<180
        toCompareTTrend(i, 1)=180-toCompareTTrend(i, 1);
    else
        toCompareTTrend(i, 1)=540-toCompareTTrend(i, 1);
    end
end

%calculate alternative trends for flipped axes (used to determine the
%minimum angle between the axes)
pressureTrendFlip=pressureTrendvector-180;
tensionTrendFlip=tensionTrendvector-180;
toComparePTrendFlip=toComparePTrend-180;
toCompareTTrendFlip=toCompareTTrend-180;
for i=1:length(pressureTrendFlip)
    if pressureTrendFlip(i, 1)<0
        pressureTrendFlip(i, 1)=360+pressureTrendFlip(i, 1);
    end
    if tensionTrendFlip(i, 1)<0
        tensionTrendFlip(i, 1)=360+tensionTrendFlip(i, 1);
    end
    if toComparePTrendFlip(i, 1)<0
        toComparePTrendFlip(i, 1)=360+toComparePTrendFlip(i, 1);
    end
    if toCompareTTrendFlip(i, 1)<0
        toCompareTTrendFlip(i, 1)=360+toCompareTTrendFlip(i, 1);
    end
end

%calculate angle between
vectorP=[(sind(pressure(:, 2)+90).* cosd(pressureTrendvector(:, 1))), (sind(pressure(:, 2)+90).* sind(pressureTrendvector(:, 1))), cosd(pressure(:, 2)+90)];
vectorT=[(sind(tension(:, 2)+90).* cosd(tensionTrendvector(:, 1))), (sind(tension(:, 2)+90).* sind(tensionTrendvector(:, 1))), cosd(tension(:, 2)+90)];
altvectorP=[(sind(-pressure(:, 2)+90).* cosd(pressureTrendFlip(:, 1))), (sind(-pressure(:, 2)+90).* sind(pressureTrendFlip(:, 1))), cosd(-pressure(:, 2)+90)];
altvectorT=[(sind(-tension(:, 2)+90).* cosd(tensionTrendFlip(:, 1))), (sind(-tension(:, 2)+90).* sind(tensionTrendFlip(:, 1))), cosd(-tension(:, 2)+90)];

toComparePvector=[(sind(toCompareP(:, 2)+90).* cosd(toComparePTrend(:, 1))), (sind(toCompareP(:, 2)+90).* sind(toComparePTrend(:, 1))), cosd(toCompareP(:, 2)+90)];
toCompareTvector=[(sind(toCompareT(:, 2)+90).* cosd(toCompareTTrend(:, 1))), (sind(toCompareT(:, 2)+90).* sind(toCompareTTrend(:, 1))), cosd(toCompareT(:, 2)+90)];
altTocompareP=[(sind(-toCompareP(:, 2)+90).* cosd(toComparePTrendFlip(:, 1))), (sind(-toCompareP(:, 2)+90).* sind(toComparePTrendFlip(:, 1))), cosd(-toCompareP(:, 2)+90)];
altTocompareT=[(sind(-toCompareT(:, 2)+90).* cosd(toCompareTTrendFlip(:, 1))), (sind(-toCompareT(:, 2)+90).* sind(toCompareTTrendFlip(:, 1))), cosd(-toCompareT(:, 2)+90)];
anglebetweenP=zeros(length(toCompareP), 1); anglebetweenPalt1=anglebetweenP; anglebetweenPalt2=anglebetweenP; anglebetweenPalt3=anglebetweenP;
anglebetweenT=zeros(length(toCompareP), 1); anglebetweenTalt1=anglebetweenT; anglebetweenTalt2=anglebetweenT; anglebetweenTalt3=anglebetweenT;
for i=1:length(toCompareP)
    if toCompareP(i, 1)~=0
        anglebetweenP(i, 1) = atan2d(norm(cross(vectorP(i, :), toComparePvector(i, :))),dot(vectorP(i, :), toComparePvector(i, :)));
        anglebetweenPalt1(i, 1) = atan2d(norm(cross(vectorP(i, :), altTocompareP(i, :))),dot(vectorP(i, :), altTocompareP(i, :)));
        anglebetweenPalt2(i, 1) = atan2d(norm(cross(altvectorP(i, :), toComparePvector(i, :))),dot(altvectorP(i, :), toComparePvector(i, :)));
        anglebetweenPalt3(i, 1) = atan2d(norm(cross(altvectorP(i, :), altTocompareP(i, :))),dot(altvectorP(i, :), altTocompareP(i, :)));
    end
end

for i=1:length(toCompareT)
    if toCompareT(i, 1)~=0
        anglebetweenT(i, 1) = atan2d(norm(cross(vectorT(i, :), toCompareTvector(i, :))), dot(vectorT(i, :), toCompareTvector(i, :)));
        anglebetweenTalt1(i, 1) = atan2d(norm(cross(vectorT(i, :), altTocompareT(i, :))), dot(vectorT(i, :), altTocompareT(i, :)));
        anglebetweenTalt2(i, 1) = atan2d(norm(cross(altvectorT(i, :), toCompareTvector(i, :))), dot(altvectorT(i, :), toCompareTvector(i, :)));
        anglebetweenTalt3(i, 1) = atan2d(norm(cross(altvectorT(i, :), altTocompareT(i, :))), dot(altvectorT(i, :), altTocompareT(i, :)));
    end
end

%define misfit angles as the minimum angle between axes
for i=1:length(anglebetweenP)
    anglebetweenPnew(i, 1)=min([anglebetweenP(i, 1), anglebetweenPalt1(i, 1), anglebetweenPalt2(i, 1), anglebetweenPalt3(i, 1)]);
    anglebetweenTnew(i, 1)=min([anglebetweenT(i, 1), anglebetweenTalt1(i, 1), anglebetweenTalt2(i, 1), anglebetweenTalt3(i, 1)]);
end
%-------------------------------------------------------------------------
%calculate degree of similarity between the model and Karson data

%degree of fit as the number of faults that are relatively accurately
%modeled for both P and T (fitBoth) and for only 1 of P or T (fit1)
anglesbetween=[anglebetweenPnew, anglebetweenTnew];
count=0; count2=0; 
for i=1:length(ends)
    if anglebetweenPnew(i, 1)>0
        angles=anglesbetween(i, :);
        if length(nonzeros(logical(angles<35)))==1
            count=count+1;
        end
        if length(nonzeros(logical(angles<35)))==2
            count2=count2+1;
        end
    end
end
fit1=count;
fitBoth=count2;

%subtract overcount due to multiple faults representing one region
numbers=[25, 26, 0, 0; 27, 28, 0, 0; 29, 30, 0, 0; 5, 34, 35, 37]; %defined numbers of faults where multiple are in a single region 
for i=1:size(numbers, 1)
    [fitBoth, fit1]=overcount(nonzeros(numbers(i, :)).', anglebetweenPnew, anglebetweenTnew, fitBoth, fit1);
end


%degree of misfit as the average angle between Ps and Ts for the two vectors
count=0;
count2=0;
for i=1:length(ends)
    if toCompareP(i, 1)~=0
        count=count+(anglebetweenPnew(i, 1)+anglebetweenTnew(i, 1));
        count2=count2+2;
    end
end
degreeOfMisfit=count/(count2);
end

