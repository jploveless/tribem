function[stressAxes, anglesbetween, misfit, stressRatio, flips, sigmaflipped, rakes] = stressTensorsCoordinates( faults, faultsHotspot, d, bc, alpha, stressGrid)
%stressTensorsCoordinates solves for angular misfit with comparison stress tensor axes
%defined at a series of observation coordinates

%Inputs: 
% faults: meshed fault map       faultsHotspot: meshed fault map and spherical mesh hotspot
% d: nelx3 array of boundary conditions & bc: NEL-by-3 array specifying type
% of boundary conditions, see tribemx for specifications 
% alpha: angle of uniaxial remote stress, defined clockwise from due east
% stressGrid: comparison data: longitude/latitude (columns 1-2),
% depth (column 3), sigma 1 (columns 4:5), sigma 2 (columns 6:7),
% sigma 3 (columns 8:9)


%Outputs: 
%stressAxes: stress tensor axes for each observation coordinate defined for
%the model
%anglesbetween: angular misfits between pressure and tension axes for each
%comparison data point
%misfit: average overall angular misfit
%stressRatio: calculated deviatoric stress ratios
%flips: text marking down identified possible stress permutations
%sigmaflipped: stress tensor axes adjusted for possible stress permutations
%rakes: overall rake calculated for each fault from the model 


%calculate REMS based on desired angle clockwise from due east
rotation=zeros(2,2);
rotation(1,1)=cosd(alpha); rotation(2,2)=cosd(alpha);
rotation(2,1)=sind(alpha); rotation(1,2)=-sind(alpha);
a=[6e6 0; 0 0];
remsinitial=rotation'*a*rotation;
REMS=zeros(1,6);
REMS(1)=remsinitial(1,1); REMS(2)=remsinitial(2,2); 
REMS(4)=-remsinitial(1,2); 

%calculate starting inputs
ends=cumsum(faults.nEl); begs=[1; ends(1:end-1)+1]; 
stressesComparison=stressGrid(:,4:9);

%calculate grid of observation coordinates
stressGrid2=[utmconv(stressGrid(:, 1:2), 0), stressGrid(:,3)];
obs.x=stressGrid2(:,1); obs.y=stressGrid2(:,2); 
obs.z=stressGrid2(:,3)*-1000; 
obs.v = 's'; 
[slip, trac, out] = tribemx(faultsHotspot, d, bc, REMS, obs);
faultsHotspot=rake(faultsHotspot, slip);
rakes=faultsHotspot.sliprake;
%calculate stress tensor axes
mstress=zeros(length(obs.x), 6); pstresses=[]; magnitudes=[];
for i=1:length(obs.x)
    mstress(i, :)=-out.s(i, 1:6);
    pstress = PrincipalStress([mstress(i, 1), mstress(i, 4), mstress(i, 5); mstress(i, 4), mstress(i, 2), mstress(i, 6); mstress(i, 5), mstress(i, 6), mstress(i, 3)], pi/2, 0, pi/2);
    pstresses=[pstresses; pstress(:, 2:3)];
    magnitudes=[magnitudes; pstress(:, 1)];
end

%calculate stress ratio
magnitudes=reshape(magnitudes, [3, length(magnitudes)/3]).';
stressRatio=zeros(length(magnitudes), 1);
deviatoricStress=zeros(length(magnitudes), 3);
for i=1:length(magnitudes)
    for j=1:3
        meanMag=mean(magnitudes(i, :));
        deviatoricStress(i, j)=magnitudes(i, j)-meanMag;
    end
end
for i=1:length(magnitudes)
    stressRatio(i, 1)=(deviatoricStress(i, 2)-deviatoricStress(i, 3))/(deviatoricStress(i, 1)-deviatoricStress(i, 3));
end

%adjust formatting for plunges (must be positive in order to come out
%correctly), trends must then be switched by 180' to remain consistent
%with which side of plunge is being considered
stressAxes=[];
for j=1:3
    k=pstresses(j:3:end, :); %take either sigma1, sigma2, or sigma3
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
    stressAxes=[stressAxes, k];
end
%extract individual stress tensor axes and convert to degrees for ease of comparison
stressAxes=rad2deg(stressAxes);
sigma1=stressAxes(:, 1:2);
sigma2=stressAxes(:, 3:4);
sigma3=stressAxes(:, 5:6);
%-------------------------------------------------------------------------
%calculate angles between comparison data and principal stresses calculated
%above
%fix trend angle to be in correct format (0'=due east)
%for conversion to cartesian coordinates
sigmaTrends=mod(stressAxes(:, 1:2:5)+90, 360);

for j=1:3 %repeat for each component of principal stress
    for i=1:length(sigmaTrends)
        if sigmaTrends(i, j)<180
            sigmaTrends(i, j)=180-sigmaTrends(i, j);
        else
            sigmaTrends(i, j)=540-sigmaTrends(i, j);
        end
    end
end

%convert modeled stress axis vectors to Cartesian format
vector1=[(sind(sigma1(:, 2)+90).* cosd(sigmaTrends(:, 1))), (sind(sigma1(:, 2)+90).* sind(sigmaTrends(:, 1))), cosd(sigma1(:, 2)+90)];
vector2=[(sind(sigma2(:, 2)+90).* cosd(sigmaTrends(:, 2))), (sind(sigma2(:, 2)+90).* sind(sigmaTrends(:, 2))), cosd(sigma2(:, 2)+90)];
vector3=[(sind(sigma3(:, 2)+90).* cosd(sigmaTrends(:, 3))), (sind(sigma3(:, 2)+90).* sind(sigmaTrends(:, 3))), cosd(sigma3(:, 2)+90)];
vectors=[vector1, vector2, vector3];

%do the same for comparison data
stressesTrends=[];
for k=1:size(stressesComparison, 2)/6
    stressesTrends=[stressesTrends, mod(stressesComparison(:, (6*k)-5:2:(6*k)-1)+90, 360)];
end
for j=1:size(stressesComparison, 2)/2
    for i=1:length(stressesTrends)
        if stressesTrends(i, j)<180
            stressesTrends(i, j)=180-stressesTrends(i, j);
        else
            stressesTrends(i, j)=540-stressesTrends(i, j);
        end
    end
end


%calculate angle between them
[anglesbetween, flips, sigmaflipped]=permutationsCheck(vectors, stressesComparison(:, 1:6), stressesTrends(:, 1:3), sigma1, sigma2, sigma3);



%adjust to consider the smallest angle between axes
for i=1:length(stressesComparison)
    for j=1:3
        if anglesbetween(i, j)>90 
            anglesbetween(i, j)=180-anglesbetween(i, j);
        end
    end
end
%-------------------------------------------------------------------------
%calculate degree of similarity with comparison data
%degree of fit as the number of faults that fit well for 1 principal stress
count=0; counts=[];
for j=1
    for i=1:length(ends)
        if anglesbetween(i, (3*j)-2)>0
            angles=anglesbetween(i, (3*j)-2:(3*j));
            if length(nonzeros(logical(angles<35)))==1
                count=count+1;
            end
        end
    end
    counts=[counts, count];
end
fit1=counts;

%degree of fit as the number of faults that fit well for 2 of 3 principal stresses
count=0; counts=[];
for j=1
    for i=1:length(ends)
        if anglesbetween(i, (3*j)-2)>0
            angles=anglesbetween(i, (3*j)-2:(3*j));
            if length(nonzeros(logical(angles<35)))==2
                count=count+1;
            end
        end
    end
    counts=[counts, count];
end
fit2=counts;

%degree of fit as the number of faults that fit well for all 3 principal stresses
count=0; counts=[];
for j=1
    for i=1:length(ends)
        if anglesbetween(i, (3*j)-2)>0 && anglesbetween(i, (3*j)-2)<35 && anglesbetween(i, (3*j)-1)<35 && anglesbetween(i, (3*j))<35
            count=count+1;
        end
    end
    counts=[counts, count];
end
fit3=counts;
fit=[fit1; fit2; fit3];


%degree of misfit as the average angle between principal stresses for the 3
%vectors
count=0;
count2=0; misfits=[];
for j=1
    for i=1:length(ends)
        if anglesbetween(i, (3*j)-2)~=0
            count=count+anglesbetween(i, (3*j)-2)+anglesbetween(i, (3*j)-1);
            count2=count2+2;
        end
    end
    misfit=count/(count2);
    misfits=[misfits, misfit];
    count=0;  count2=0;
end
misfit=mean(misfits);
anglesbetween=[anglesbetween(:,1), anglesbetween(:,3)]; 
end