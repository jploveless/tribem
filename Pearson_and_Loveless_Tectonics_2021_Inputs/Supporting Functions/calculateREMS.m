function[REMS, remsinitial]=calculateREMS(alpha)
%calculate uniaxial remote stress tensor at a given azimuth, given as an
%angle counterclockwise from due East
rotation=zeros(2,2);
rotation(1,1)=cosd(alpha); rotation(2,2)=cosd(alpha);
rotation(2,1)=sind(alpha); rotation(1,2)=-sind(alpha);
a=[6e6 0; 0 0];
remsinitial=rotation'*a*rotation;
REMS=zeros(1,6);
REMS(1)=remsinitial(1,1); REMS(2)=remsinitial(2,2); REMS(4)=remsinitial(1,2);
end