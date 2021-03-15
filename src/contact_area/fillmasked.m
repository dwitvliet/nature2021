function fseg=fillmasked(seg,mask)
  %This function fills the 0 gaps in seg from nonzero regions, but only where mask is not 0.
  %seg and mask have to be 2d arrays of same size.
  %By Daniel Berger for the CElegans project (Aravi Samuel / Mei Zhen), December 2017
  
  fseg=[];
  
  if ((size(seg,1)~=size(mask,1))||(size(seg,2)~=size(mask,2)))
    %ERROR: seg and mask have different sizes. return
    return;
  end;
  
  seg(mask==0)=0;
  u=unique(seg);
  u(u==0)=[];
  udone=u*0;
  
  if (min(size(u))==0)
    %Error: masked region is empty. return
    return;
  end;
  
  done=0;
  fseg=seg;
  while (done==0)
    %tic;
    %randomize nonzero IDs
    r=randperm(length(u));
    ru=u(r);
    
    done=1;
    for i=1:length(ru)
      if (udone(r(i))==0)
        done=0; %as long as at least one ID is still expanding, we are not done.
        
        %single out each id using randomized order and expand by 1 pixel up,down,left,right
        id=ru(i); %get ID
        iseg=(fseg==id); %mask out ID region
        
        %find bounding box
        xmin=find(max(iseg,[],1)>0,1,'first');
        xmax=find(max(iseg,[],1)>0,1,'last');
        ymin=find(max(iseg,[],2)>0,1,'first');
        ymax=find(max(iseg,[],2)>0,1,'last');
        
        xmin=max([xmin-1 1]);
        xmax=min([xmax+1 size(iseg,2)]);
        ymin=max([ymin-1 1]);
        ymax=min([ymax+1 size(iseg,1)]);
        
        ciseg=iseg(ymin:ymax,xmin:xmax);
        cmask=mask(ymin:ymax,xmin:xmax);
        cfseg=fseg(ymin:ymax,xmin:xmax);
        
        oiseg=ciseg; %store previous region for comparison
        
        %expand ID region
        ciseg(1:end-1,:)=max(ciseg(1:end-1,:),ciseg(2:end,:));
        ciseg(2:end,:)=max(ciseg(2:end,:),ciseg(1:end-1,:));
        ciseg(:,1:end-1)=max(ciseg(:,1:end-1),ciseg(:,2:end));
        ciseg(:,2:end)=max(ciseg(:,2:end),ciseg(:,1:end-1));
        
        ciseg(cmask==0)=0; %constrain to mask
        cfseg(cfseg==0)=uint32(ciseg(cfseg==0))*id; %paste back only into 0 pixels in cfseg
        fseg(ymin:ymax,xmin:xmax)=cfseg; %paste back to full image
        
        %compare previous and new region; if they are the same, this ID is done
        niseg=(cfseg==id);
        if (max(max(abs(oiseg-niseg)))==0)
          udone(r(i))=1;
        end;
      end;
    end;
    %toc
  end;