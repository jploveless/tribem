function [pressures, tensions]=multipleFaults(pressures, tensions, begs, ends, strikesAndDips, trends, plunges, slipmag, faultNumbers)
%calculate pressure and tension axes for region containing multiple faults

total=[];
for i=faultNumbers
    total=[total, begs(i):ends(i)];
end
[pressures(faultNumbers(1), :), tensions(faultNumbers(1), :)]=MomTens(strikesAndDips(total, 1:2), [deg2rad(trends(total)), deg2rad(plunges(total))], slipmag(total));

for i=faultNumbers(2:end)
    pressures(i, :)=pressures(faultNumbers(1), :);
    tensions(i, :)=tensions(faultNumbers(1), :);
end
