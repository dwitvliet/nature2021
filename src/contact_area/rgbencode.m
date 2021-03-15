function dimg=rgbencode(img)
%Encoding of RGB label images from uint32 (max. 24 bits used)
%By Daniel Berger, December 2013

dimg=zeros(size(img,1),size(img,2),3,'uint8');

dimg(:,:,3)=uint8(bitand(img,255));
dimg(:,:,2)=uint8(bitand(bitshift(img,-8),255));
dimg(:,:,1)=uint8(bitand(bitshift(img,-16),255));
