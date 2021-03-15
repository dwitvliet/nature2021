function ret=getprev(data,refid)
%Returnd id of prev of refid
%Uses the 24-column-data matrix as analyzed from VAST color files

ret=data(refid,16);
