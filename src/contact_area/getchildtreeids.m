function ret=getchildtreeids(data,parentlist)
%Uses the 24-column-data matrix as analyzed from VAST color files
%Gets a list of the IDs of the segment's children's tree (if it exists)

ret=[];
pal=parentlist(:);
for p=1:1:size(pal,1)
  index=parentlist(p);
  if (data(index,15)>0)
    i=data(index,15);
    ret=[ret i getchildtreeids(data,i)]; %Add size of child tree
    while (data(i,17)>0) %add sizes of all nexts
      i=data(i,17);
      ret=[ret i getchildtreeids(data,i)];
    end;
  end;
end;