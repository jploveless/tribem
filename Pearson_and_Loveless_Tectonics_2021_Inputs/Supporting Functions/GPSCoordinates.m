function [anglesbetween, misfit, results, rakes] = GPSCoordinates(faults, faultsHotspot, d, bc, azimuth, hotspotContribution, GPSData)
%GPScomparison evaluates model results in comparison with strain axes
%derived from GPS velocity data 
%Inputs: 
% faults: meshed fault map       faultsHotspot: meshed fault map and spherical mesh hotspot
% d: nelx3 array of boundary conditions & bc: NEL-by-3 array specifying type
% of boundary conditions, see tribemx for specifications 
% azimuth: angle of uniaxial remote stress, defined clockwise from due east
% hotspot Contribution: amount of element-normal dislocation on the hotspot
% in mm/yr
% GPSData: comparison data, including: stationsTotal (UTM coordinates of GPS stations, column 1: 
% Eastings, column 2: Northings), GPS (trends of axes of maximum and minimum extensional strain calculated
% directly from GPS velocities, 1:max, 2:min), and stationsList(UTM coordinates, broken down by region in order:
% NW, CW, SW, NC, CC, SC, NE, CE, SE)


%Outputs: 
%anglesbetween: angular misfits between axis of maximum extension for each
%comparison data point
%misfit: overall average misfit for model 
%results:  axes of maximum extension for each observation coordinate defined for
%the model
%rakes: overall rake calculated for each fault from the model 

%set up grid of observation points at each GPS station location
stationsTotal=GPSData{1};
GPS=GPSData{2};
stationsList=GPSData{3};
xyutm=GPSData{4};
obs.x=stationsTotal(:, 1); obs.y=stationsTotal(:, 2); obs.z=0*obs.x;
obs.v="d";
obs.rc=[xyutm(1,1), xyutm(1,2), 0]; %reference coordinates, close to the plate boundary

%run tribemx to calculate displacements at each station
d(sum(faults.nEl)+1:end, 3)=hotspotContribution;
[REMSrotate, ~]=calculateREMS(azimuth);
[slip, trac, out] = tribemx(faultsHotspot, d, bc, REMSrotate, obs);
disp2=[out.u(:,1), out.u(:,2)];
faultsHotspot=rake(faultsHotspot, slip);
rakes=faultsHotspot.sliprake;


start=0; type=1; plotpar=0;
array=zeros(9, 1);  arrayAlt=zeros(9, 1); misfitAlt=zeros(9, 1);
misfit=zeros(9, 1); anglesbetween=zeros(9,1); results=zeros(9,1);
for k=1:9  %run for each of the 9 regions
    %calculate strains based on displacement
    stationsNew=stationsList{1, k};
    nStations=length(stationsNew);
    par = [500e3, nStations, 250e3];
    [~,~,~,pstrains,~] = GridStrain(stationsNew,disp2(start+1:start+length(stationsNew), :),type,par,plotpar);
    start=start+length(stationsNew);
    %rename
    e1mag = pstrains(1, 1, :); e1mag = e1mag(:);
    e1tre = pstrains(1, 2, :); e1tre = e1tre(:);
    e1plu = pstrains(1, 3, :); e1plu = e1plu(:);
    e2mag = pstrains(2, 1, :); e2mag = e2mag(:);
    e2tre = pstrains(2, 2, :); e2tre = e2tre(:);
    e2plu = pstrains(2, 3, :); e2plu = e2plu(:);
    e3mag = pstrains(3, 1, :); e3mag = e3mag(:);
    e3tre = pstrains(3, 2, :); e3tre = e3tre(:);
    e3plu = pstrains(3, 3, :); e3plu = e3plu(:);
    
    %permutations
    e1horizmag = e1mag; % Default assumption is that e1 horizontal is e1
    e1horiztre = e1tre;
    e1horizplu = e1plu;
    vert_e1 = e1plu == pi/2; % Logical index of any e1 plunges that are vertical
    e1horizmag(vert_e1) = e2mag(vert_e1); % Reassign those values to be e2 values
    e1horiztre(vert_e1) = e2tre(vert_e1);
    e1horizplu(vert_e1) = e2plu(vert_e1);
    e3horizmag = e3mag; % Default assumption is that e3 horizontal is e3
    e3horiztre = e3tre;
    e3horizplu = e3plu;
    vert_e3 = e3plu == pi/2; % Logical index of any e3 plunges that are vertical
    e3horizmag(vert_e3) = e2mag(vert_e3); % Reassign those values to be e2 values
    e3horiztre(vert_e3) = e2tre(vert_e3);
    e3horizplu(vert_e3) = e2plu(vert_e3);
    array(k, 1)=rad2deg(e1horiztre);
    arrayAlt(k, 1)=mod(rad2deg(e1horiztre)+180, 360);
    misfit(k, 1)=mod(abs(array(k, 1)-GPS(k, 1)), 180);
    misfitAlt(k, 1)=mod(abs(arrayAlt(k, 1)-GPS(k, 1)), 180);
    anglesbetween(k, 1)=min(misfit(k,1), misfitAlt(k,1));
    results(k,1)=array(k,1);
end
results(:,2)=zeros(9,1);
misfit=mean(anglesbetween);
end

