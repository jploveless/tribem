function[]=numberedFaults(all)
ends=cumsum(all.nEl);
begs=[1; ends(1:end-1)+1];
cols=zeros(size(all.v, 1), 1);
for i=1:length(ends)
    cols(begs(i):ends(i)) = i;
end
meshview(all.c, all.v);
text(all.x1, all.y1, num2str(cols));
end

