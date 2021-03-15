function ret=getparent(data,refid)
%Returnd id of parent of refid
%Uses the 24-column-data matrix as analyzed from VAST color files

ret=data(refid,14);
