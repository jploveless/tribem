function[misfit, anglesbetween, stressAxes, stressRatio, flips, sigmaFlipped, rakes] = stressTensorsFaults( faults, faultsHotspot, d, bc, alpha, stressData)
%stressTensorsFaults solves for angular misfit with comparison stress tensor 
%defined at a series of faults

%Inputs: 
% faults: meshed fault map       faultsHotspot: meshed fault map and spherical mesh hotspot
% d: nelx3 array of boundary conditions & bc: NEL-by-3 array specifying type
% of boundary conditions, see tribemx for specifications 
% alpha: angle of uniaxial remote stress, defined clockwise from due east
% G: pre-calculated Green's functions, derived from tribemx
% stressData: comparison data: sigma 1 (columns 1:2), sigma 2(columns 3:4),
% sigma 3(columns 5:6), if data contains multiple comparison points for 1
% fault: pattern continues, with sigma 1 (columns 7:8), etc. 


%Outputs: 
%misfit: average overall angular misfit
%fit: number of faults for which 1, 2, and all of the principal stress axes
%fit well 
%anglesbetween: angular misfits between sigma 1, 2, and 3 axes for each
%comparison data point
%stressAxes: stress tensor axes for each fault defined for
%the model
%stressRatio: calculated deviatoric stress ratios
%flips: text marking down identified possible stress permutations
%sigmaflipped: stress tensor axes adjusted for possible stress permutations

%calculate REMS based on angle counting clockwise from due East
rotation=zeros(2,2);
rotation(1,1)=cosd(alpha); rotation(2,2)=cosd(alpha);
rotation(2,1)=sind(alpha); rotation(1,2)=-sind(alpha);
a=[6e6 0; 0 0];
remsinitial=rotation'*a*rotation;
REMS=zeros(1,6);
REMS(1)=remsinitial(1,1); REMS(2)=remsinitial(2,2); 
REMS(4)=-remsinitial(1,2); 

%calculate starting inputs
stressesComparison=stressData(:,5:end);
ends=cumsum(faults.nEl); begs=[1; ends(1:end-1)+1]; 

%calculate grid of observation coordinates
stressGrid2=[utmconv(stressData(:, 1:2), 0), stressData(:,3)];
obs.x=stressGrid2(:,1); obs.y=stressGrid2(:,2); 
obs.z=stressGrid2(:,3)*-1000;
obs.v = 's';
[slip, trac, out, G] = tribemx(faultsHotspot, d, bc, REMS, obs);
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

%calculate deviatoric stress ratio
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

%convert modeled principal stress axes to Cartesian format
vector1=[(sind(sigma1(:, 2)+90).* cosd(sigmaTrends(:, 1))), (sind(sigma1(:, 2)+90).* sind(sigmaTrends(:, 1))), cosd(sigma1(:, 2)+90)];
vector2=[(sind(sigma2(:, 2)+90).* cosd(sigmaTrends(:, 2))), (sind(sigma2(:, 2)+90).* sind(sigmaTrends(:, 2))), cosd(sigma2(:, 2)+90)];
vector3=[(sind(sigma3(:, 2)+90).* cosd(sigmaTrends(:, 3))), (sind(sigma3(:, 2)+90).* sind(sigmaTrends(:, 3))), cosd(sigma3(:, 2)+90)];
vectors=[vector1, vector2, vector3];

%do the same for comparison data
plateauxTrends=[];
for k=1:size(stressesComparison, 2)/6
    plateauxTrends=[plateauxTrends, mod(stressesComparison(:, (6*k)-5:2:(6*k)-1)+90, 360)];
end
for j=1:size(stressesComparison, 2)/2
    for i=1:length(plateauxTrends)
        if plateauxTrends(i, j)<180
            plateauxTrends(i, j)=180-plateauxTrends(i, j);
        else
            plateauxTrends(i, j)=540-plateauxTrends(i, j);
        end
    end
end


%calculate angle between them
[anglebetween1, flips1, sigma1flipped]=permutationsCheck(vectors, stressesComparison(:, 1:6), plateauxTrends(:, 1:3), sigma1, sigma2, sigma3);
[anglebetween2, flips2, sigma2flipped]=permutationsCheck(vectors, stressesComparison(:, 7:12), plateauxTrends(:, 4:6), sigma1, sigma2, sigma3);
[anglebetween3, flips3, sigma3flipped]=permutationsCheck(vectors, stressesComparison(:, 13:18), plateauxTrends(:, 7:9), sigma1, sigma2, sigma3);
if size(stressesComparison, 2)/6==5
    [anglebetween4, flips4, sigma4flipped]=permutationsCheck(vectors, stressesComparison(:, 19:24), plateauxTrends(:, 10:12), sigma1, sigma2, sigma3);
    [anglebetween5, flips5, sigma5flipped]=permutationsCheck(vectors, stressesComparison(:, 25:30), plateauxTrends(:, 13:15), sigma1, sigma2, sigma3);
end
flips={flips1; flips2; flips3};
sigmaFlipped={sigma1flipped, sigma2flipped, sigma3flipped};

if size(stressesComparison, 2)/6==5
    flips={flips1; flips2; flips3; flips4; flips5};
    sigmaFlipped={sigma1flipped, sigma2flipped, sigma3flipped, sigma4flipped, sigma5flipped};%, sigma6flipped, sigma7flipped};
end

%adjust to consider the smallest angle between axes
for i=1:length(stressesComparison)
    for j=1:3
        if anglebetween1(i, j)>90 
            anglebetween1(i, j)=180-anglebetween1(i, j);
        end
        if anglebetween2(i, j)>90 
            anglebetween2(i, j)=180-anglebetween2(i, j);
        end
        if anglebetween3(i, j)>90 
            anglebetween3(i, j)=180-anglebetween3(i, j);
        end
        if size(stressesComparison,2)/6==5
            if anglebetween4(i, j)>90 
                anglebetween4(i, j)=180-anglebetween4(i, j);
            end
            if anglebetween5(i, j)>90 
                anglebetween5(i, j)=180-anglebetween5(i, j);
            end
        end
    end
end

anglesbetween=[anglebetween1, anglebetween2, anglebetween3];
if size(stressesComparison, 2)/6==5
    anglesbetween=[anglebetween1, anglebetween2, anglebetween3, anglebetween4, anglebetween5];
end
%-------------------------------------------------------------------------
%calculate degree of similarity with comparison data
%degree of misfit as the average angle between principal stresses for the 2
%vectors
anglesbetween=anglesbetween(:,[1,3,4,6,7,9,10,12,13,15]);
anglesbetween2=anglesbetween;
for i=1:size(stressData,1)
    anglesbetween2(i,11)=sum(anglesbetween2(i,1:10))/length(nonzeros(anglesbetween2(i, 1:10)));
end
misfit=sum(anglesbetween2(:,11), 'omitnan')/size(stressData,1);
end