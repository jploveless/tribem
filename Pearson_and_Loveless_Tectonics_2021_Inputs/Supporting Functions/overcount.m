function[degreeOfFitBoth, degreeOfFit1]=overcount(faultNumbers, anglebetweenP, anglebetweenT, degreeOfFitBoth, degreeOfFit1)
%subtract overcount due to multiple faults associated with 1 region in the
%fault-slip data 
if anglebetweenP(faultNumbers(1), 1)<35
    if anglebetweenT(faultNumbers(1), 1)<35
        degreeOfFitBoth=degreeOfFitBoth-length(faultNumbers)+1;
    end
end
if anglebetweenP(faultNumbers(1), 1)<35
    if anglebetweenT(faultNumbers(1), 1)>=35
        degreeOfFit1=degreeOfFit1-length(faultNumbers)+1;
    end
elseif anglebetweenT(faultNumbers(1), 1)<35
    degreeOfFit1=degreeOfFit1-length(faultNumbers)+1;
end

end

