function [] = makeStereonet(results, type)
%makeStereonet plots an equal area stereonet, plotting axes on top 
figure
Stereonet(0,90*pi/180,10*pi/180,1);
hold on
results=deg2rad(results);
for i=1:length(results)
    %Plot sigma1/P axis/maximum extension (black, filled circle)
    [xp,yp] = StCoordLine(results(i,1),results(i,2),1);
    h2=plot(xp,yp,'ko','MarkerFaceColor','k');
    %Plot T axis (red circle)
    if size(results, 2)==4 
        [xp,yp] = StCoordLine(results(i,3),results(i,4),1);
        h3=plot(xp,yp,'ro','MarkerFaceColor','r');
    elseif size(results, 2)==6
        [xp,yp] = StCoordLine(results(i,3),results(i,4),1);
        h3=plot(xp,yp,'bo','MarkerFaceColor','b');
        [xp,yp] = StCoordLine(results(i,5),results(i,6),1);
        h4=plot(xp,yp,'ro','MarkerFaceColor','r');
    end
end
legend('Location','northwest');
if type==0 || type==2
    legend([h2 h3],{'P', 'T'});
elseif type==1 || type==3
    legend([h2 h3 h4],{'σ1','σ2', 'σ3'});
else
    legend(h2,{'Extension'}); 
end
end