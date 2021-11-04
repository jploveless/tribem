function[anglesbetween, flips, sigma]=permutationsCheck(modeledVectors, stressAxes, stressAxesTrends, sigma1, sigma2, sigma3)
%checks principal stress axes for any possible stress permutations 

%prepare variables 
anglesbetween1=zeros(length(stressAxes), 3);
anglesbetween2=anglesbetween1; anglesbetween3=anglesbetween1;
anglesbetween=anglesbetween1;
toCompare=zeros(length(stressAxes), 9);
flips=[];

for i=1:length(stressAxes)
    if stressAxes(i, 1)~=0 %if there is comparison data for this fault
        %collect vectors for this fault
        vector1=modeledVectors(i, 1:3); 
        vector2=modeledVectors(i, 4:6); 
        vector3=modeledVectors(i, 7:9);
        
        %convert comparison data to Cartesian format
        toCompare1(i, 1:3)=[(sind(stressAxes(i, 2)+90).* cosd(stressAxesTrends(i, 1))), (sind(stressAxes(i, 2)+90).* sind(stressAxesTrends(i, 1))), cosd(stressAxes(i, 2)+90)];
        toCompare2(i, 1:3)=[(sind(stressAxes(i, 4)+90).* cosd(stressAxesTrends(i, 2))), (sind(stressAxes(i, 4)+90).* sind(stressAxesTrends(i, 2))), cosd(stressAxes(i, 4)+90)];
        toCompare3(i, 1:3)=[(sind(stressAxes(i, 6)+90).* cosd(stressAxesTrends(i, 3))), (sind(stressAxes(i, 6)+90).* sind(stressAxesTrends(i, 3))), cosd(stressAxes(i, 6)+90)];
        
        %calculate angles between each of the modeled principal stress axes
        %and the axes in the comparison data 
        anglesbetween1(i,1)=acosd(dot(vector1, toCompare1(i, :))/(norm(vector1)*norm(toCompare1(i, :))));
        anglesbetween1(i,2)=acosd(dot(vector1, toCompare2(i, :))/(norm(vector1)*norm(toCompare2(i, :))));
        anglesbetween1(i,3)=acosd(dot(vector1, toCompare3(i, :))/(norm(vector1)*norm(toCompare3(i, :))));
        anglesbetween2(i,1)=acosd(dot(vector2, toCompare1(i, :))/(norm(vector2)*norm(toCompare1(i, :))));
        anglesbetween2(i,2)=acosd(dot(vector2, toCompare2(i, :))/(norm(vector2)*norm(toCompare2(i, :))));
        anglesbetween2(i,3)=acosd(dot(vector2, toCompare3(i, :))/(norm(vector2)*norm(toCompare3(i, :))));
        anglesbetween3(i,1)=acosd(dot(vector3, toCompare1(i, :))/(norm(vector3)*norm(toCompare1(i, :))));
        anglesbetween3(i,2)=acosd(dot(vector3, toCompare2(i, :))/(norm(vector3)*norm(toCompare2(i, :))));
        anglesbetween3(i,3)=acosd(dot(vector3, toCompare3(i, :))/(norm(vector3)*norm(toCompare3(i, :))));
        
        %check that the smallest angle between the axes are being
        %considered
        for j=1:3
            if anglesbetween1(i, j)>90 
                anglesbetween1(i, j)=180-anglesbetween1(i, j);
            end
            if anglesbetween2(i, j)>90 
                anglesbetween2(i, j)=180-anglesbetween2(i, j);
            end
            if anglesbetween3(i, j)>90 
                anglesbetween3(i, j)=180-anglesbetween3(i, j);
            end
        end
        
        %check for existence of any stress permutations 
        %sigma1-sigma2
        if anglesbetween1(i, 1)>anglesbetween1(i, 2) && anglesbetween2(i, 2)>anglesbetween2(i, 1) && anglesbetween1(i, 2)<50 && anglesbetween2(i, 1)<50
            flips=[flips, sprintf('%d (sigma1-sigma2)', i)];
            [sigma1(i, :), sigma2(i, :)]=swap(sigma2(i, :), sigma1(i, :));
            anglesbetween(i, 1)=anglesbetween1(i, 2);
            anglesbetween(i, 2)=anglesbetween2(i, 1);
            anglesbetween(i, 3)=anglesbetween3(i, 3);
        %sigma2-sigma3
        elseif anglesbetween2(i, 2)>anglesbetween2(i, 3) && anglesbetween3(i, 3)>anglesbetween3(i, 2) && anglesbetween2(i, 3)<50 && anglesbetween3(i, 2)<50
            flips=[flips, sprintf('%d (sigma2-sigma3)', i)];
            [sigma2(i, :), sigma3(i, :)]=swap(sigma3(i, :), sigma2(i, :));
            anglesbetween(i, 2)=anglesbetween2(i, 3);
            anglesbetween(i, 3)=anglesbetween3(i, 2);
            anglesbetween(i, 1)=anglesbetween1(i, 1);
        else %no permutation 
            anglesbetween(i, 1)=anglesbetween1(i,1);
            anglesbetween(i, 2)=anglesbetween2(i, 2);
            anglesbetween(i, 3)=anglesbetween3(i, 3);
        end         
    end
end
sigma=[sigma1, sigma2, sigma3];
end
