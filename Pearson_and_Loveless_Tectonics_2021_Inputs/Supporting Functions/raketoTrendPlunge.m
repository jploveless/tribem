function [ trends, plunges ] = raketoTrendPlunge( strikes, dips, rakes )
%takes series of faults with known strike, dip and rake and outputs the
%trend and plunge for each rake
tanR=zeros(size(rakes));
cosD=tanR; sinD=tanR; sinR=tanR; beta=tanR; betaTwo=tanR; rawtrend=tanR;
trends=zeros(size(rakes));
plunges=zeros(size(rakes));
for i=1:size(rakes)
    tanR(i, 1)=tan(degtorad(rakes(i, 1)));
    cosD(i, 1)=cos(degtorad(dips(i, 1)));
    sinD(i, 1)=sin(degtorad(dips(i, 1)));
    sinR(i, 1)=sin(degtorad(rakes(i, 1)));
    beta(i, 1)=abs(radtodeg(atan(tanR(i,1)*cosD(i,1))));
    betaTwo(i,1)=180-beta(i, 1);
    if rakes(i,1)>90
        rawtrend(i,1)= strikes(i, 1)+ betaTwo(i,1);
    else
        rawtrend(i,1)=strikes(i, 1)+ beta(i, 1);
    end
    if rawtrend(i,1)>360
        trends(i,1)=rawtrend(i,1)-360;
    else 
        trends(i,1)=rawtrend(i, 1);
    end
    plunges(i,1)=radtodeg(asin(sinD(i, 1)*sinR(i, 1)));
end

end

