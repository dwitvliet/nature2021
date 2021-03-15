function ret=getnext(data,refid)
%Returnd id of next of refid
%Uses the 24-column-data matrix as analyzed from VAST color files

ret=data(refid,17);
