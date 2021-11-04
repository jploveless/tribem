function [momAxes, degreeOfMisfit, anglesbetween, rakes]= momentTensorsCoordinates( faults, faultsHotspot, d, bc, alpha, momentTensorGrid)
%momentTensorsCoordinates solves for angular misfit with comparison moment tensor axes
%defined at a series of observation coordinates

%Inputs: 
% faults: meshed fault map       faultsHotspot: meshed fault map and spherical mesh hotspot
% d: nelx3 array of boundary conditions & bc: NEL-by-3 array specifying type
% of boundary conditions, see tribemx for specifications 
% alpha: angle of uniaxial remote stress, defined clockwise from due east
% momentTensorGrid: comparison data: longitude/latitude (columns 1-2),
% depth (column 3), pressure axes (columns 4:5), tension axes (columns 6:7,
% moment magnitude (column 8)


%Outputs: 
%momAxes: moment tensor axes for each observation coordinate defined for
%the model
%degreeOfMisfit: average overall angular misfit
%anglesbetween: angular misfits between pressure and tension axes for each
%comparison data point
%rakes: overall rake calculated for each fault from the model 


%calculate starting inputs
rotation=zeros(2,2);
rotation(1,1)=cosd(alpha);
rotation(2,2)=cosd(alpha);
rotation(2,1)=sind(alpha);
rotation(1,2)=-sind(alpha);
a=[6e6 0; 0 0];
remsinitial=rotation'*a*rotation;
REMS=zeros(1,6);
REMS(1)=remsinitial(1,1);
REMS(2)=remsinitial(2,2);
REMS(4)=-remsinitial(1,2);
ends=cumsum(faults.nEl); begs=[1; ends(1:end-1)+1]; 
momComparison=momentTensorGrid(:,4:7);

%make grid of observation coodinates
momGrid=[utmconv(momentTensorGrid(:, 1:2), 0, '27W'), momentTensorGrid(:,3)];
obs.x=momGrid(:,1); obs.y=momGrid(:,2); 
obs.z=momGrid(:,3)*-1000;
obs.v = 'e'; 
[slip, trac, out] = tribemx(faultsHotspot, d, bc, REMS, obs);
faultsHotspot=rake(faultsHotspot, slip);
rakes=faultsHotspot.sliprake;
%calculate moment tensor axes
maxes=zeros(length(obs.x), 6); paxes2=[]; magnitudes=[];
for i=1:length(obs.x)
    maxes(i, :)=-out.e(i, 1:6);
    paxes = PrincipalStress([maxes(i, 1), maxes(i, 4), maxes(i, 5); maxes(i, 4), maxes(i, 2), maxes(i, 6); maxes(i, 5), maxes(i, 6), maxes(i, 3)], pi/2, 0, pi/2);
    paxes2=[paxes2; paxes(:, 2:3)];
    magnitudes=[magnitudes; paxes(:, 1)];
end

%adjust formatting for plunges (must be positive in order to come out
%correctly), trends must then be switched by 180' to remain consistent
%with which side of plunge is being considered
momAxes=[];
for j=1:3
    k=paxes2(j:3:end, :); %take either P-axis, B-axis, or T-axis
    for i=1:length(k)
        if k(i, 2)<0
            k(i, 2)=-k(i, 2);
            if k(i, 1)<pi
                k (i, 1)=k(i, 1)+pi;
            else
                k(i, 1)=k(i, 1)-pi;
            end
        end
    end
    momAxes=[momAxes, k];
end

%extract individual moment tensor axes and convert to degrees for ease of comparison
momAxes=rad2deg(momAxes);
pressure=momAxes(:, 1:2);
tension=momAxes(:, 5:6);

momAxes=[momAxes(:, 1:2), momAxes(:,5:6)];
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

toComparePTrend=mod(momComparison(:, 1)+90, 360);
for i=1:length(pressureTrendvector)
    if toComparePTrend(i, 1)<180
        toComparePTrend(i, 1)=180-toComparePTrend(i, 1);
    else
        toComparePTrend(i, 1)=540-toComparePTrend(i, 1);
    end
end

toCompareTTrend=mod(momComparison(:, 3)+90, 360);
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

toComparePvector=[(sind(momComparison(:, 2)+90).* cosd(toComparePTrend(:, 1))), (sind(momComparison(:, 2)+90).* sind(toComparePTrend(:, 1))), cosd(momComparison(:, 2)+90)];
toCompareTvector=[(sind(momComparison(:, 4)+90).* cosd(toCompareTTrend(:, 1))), (sind(momComparison(:, 4)+90).* sind(toCompareTTrend(:, 1))), cosd(momComparison(:, 4)+90)];
altTocompareP=[(sind(-momComparison(:, 2)+90).* cosd(toComparePTrendFlip(:, 1))), (sind(-momComparison(:, 2)+90).* sind(toComparePTrendFlip(:, 1))), cosd(-momComparison(:, 2)+90)];
altTocompareT=[(sind(-momComparison(:, 4)+90).* cosd(toCompareTTrendFlip(:, 1))), (sind(-momComparison(:, 4)+90).* sind(toCompareTTrendFlip(:, 1))), cosd(-momComparison(:, 4)+90)];
anglebetweenP=zeros(length(momComparison), 1); anglebetweenPalt1=anglebetweenP; anglebetweenPalt2=anglebetweenP; anglebetweenPalt3=anglebetweenP;
anglebetweenT=zeros(length(momComparison), 1); anglebetweenTalt1=anglebetweenT; anglebetweenTalt2=anglebetweenT; anglebetweenTalt3=anglebetweenT;
for i=1:length(momComparison)
    if momComparison(i, 1)~=0
        anglebetweenP(i, 1) = atan2d(norm(cross(vectorP(i, :), toComparePvector(i, :))),dot(vectorP(i, :), toComparePvector(i, :)));
        anglebetweenPalt1(i, 1) = atan2d(norm(cross(vectorP(i, :), altTocompareP(i, :))),dot(vectorP(i, :), altTocompareP(i, :)));
        anglebetweenPalt2(i, 1) = atan2d(norm(cross(altvectorP(i, :), toComparePvector(i, :))),dot(altvectorP(i, :), toComparePvector(i, :)));
        anglebetweenPalt3(i, 1) = atan2d(norm(cross(altvectorP(i, :), altTocompareP(i, :))),dot(altvectorP(i, :), altTocompareP(i, :)));
    end
end

for i=1:length(momComparison)
    if momComparison(i, 1)~=0
        anglebetweenT(i, 1) = atan2d(norm(cross(vectorT(i, :), toCompareTvector(i, :))), dot(vectorT(i, :), toCompareTvector(i, :)));
        anglebetweenTalt1(i, 1) = atan2d(norm(cross(vectorT(i, :), altTocompareT(i, :))), dot(vectorT(i, :), altTocompareT(i, :)));
        anglebetweenTalt2(i, 1) = atan2d(norm(cross(altvectorT(i, :), toCompareTvector(i, :))), dot(altvectorT(i, :), toCompareTvector(i, :)));
        anglebetweenTalt3(i, 1) = atan2d(norm(cross(altvectorT(i, :), altTocompareT(i, :))), dot(altvectorT(i, :), altTocompareT(i, :)));
    end
end

for i=1:length(anglebetweenP)
    anglebetweenPnew(i, 1)=min([anglebetweenP(i, 1), anglebetweenPalt1(i, 1), anglebetweenPalt2(i, 1), anglebetweenPalt3(i, 1)]);
    anglebetweenTnew(i, 1)=min([anglebetweenT(i, 1), anglebetweenTalt1(i, 1), anglebetweenTalt2(i, 1), anglebetweenTalt3(i, 1)]);
end
anglesbetween=[anglebetweenPnew, anglebetweenTnew];
%-------------------------------------------------------------------------
%calculate degree of similarity between the model and comparison data 
%degree of fit as the number of faults that are relatively accurately
%modeled for both P and T (fitBoth) and for only 1 of P or T (fit1)
anglesbetween=[anglebetweenPnew, anglebetweenTnew];
count=0; count2=0; 
for i=1:length(anglesbetween)
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

%degree of misfit as the average angle between Ps and Ts for the two vectors
anglesbetween2=(anglesbetween(:,1)+anglesbetween(:,2))/2;
for i=1:length(momComparison)
    if momComparison(i, 1)~=0
        count(i,1)=anglesbetween2(i,1)*momentTensorGrid(i,8)/sum(momentTensorGrid(:,8));
    end
end
degreeOfMisfit=sum(count);
end

