function dimg=rgbdecode(img)
%Decoding of RGB label images saved by VAST from RGB to uint32
%By Daniel Berger, April 2013

dimg=[];

if (size(img,3)==1)
  dimg=img;
end;

if (size(img,3)==3)
  dimg=uint32(img(:,:,1))*65536+uint32(img(:,:,2))*256+uint32(img(:,:,3));
end;