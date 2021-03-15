function [name,data,firstelement]=scanvastcolorfile(filename,replacecharsflag)
%A script to parse VAST color text files.
%By Daniel Berger, August 2013

fid = fopen(filename);
  tline = fgetl(fid);
  y=1;
  while ischar(tline)
    if ((numel(tline)>0)&&(tline(1)~='%'))
      %disp(tline);
      [a,count,errmsg,nextindex]=sscanf(tline, '%d   %d    %d %d %d %d   %d %d %d %d   %d %d %d   %d %d %d %d   %d   %d %d %d %d %d %d   ');
      data(y,:)=a';
      n=tline(nextindex:end); %Rest of line is name
      n=n(2:end-1); %Remove "" from name
      if (replacecharsflag)
        n(n==' ')='_';
        n(n=='?')='_';
        n(n=='*')='_';
        n(n=='\')='_';
        n(n=='/')='_';
        n(n=='|')='_';
        n(n==':')='_';
        n(n=='"')='_';
        n(n=='<')='_';
        n(n=='>')='_';
      end;
      name{y}=n;
      y=y+1;
    end;
    
    tline = fgetl(fid);
  end;
  fclose(fid);
  
  firstelement=data(1,17);
  data=data(2:end,:);
  name=name(2:end);