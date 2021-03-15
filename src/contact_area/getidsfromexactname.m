function ids=getidsfromexactname(namesarray,name)

[dummy,ids]=ind2sub(size(namesarray), strmatch(name, namesarray, 'exact'));