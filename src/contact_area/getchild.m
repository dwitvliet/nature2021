function ret=getchild(data,refid)
%Returnd id of first child of refid
%Uses the 24-column-data matrix as analyzed from VAST color files

ret=data(refid,15);
