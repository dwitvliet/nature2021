%celegans_sem2_computetoucharea2.m
%Masked analysis of neighbor touch areas in 2D (XY)
%By Daniel Berger for the C. elegans project (Aravi Samuel, Mei Zhen)
%December 2017, modified January - March 2018
%Modified by Daniel Witvliet, March 2018

dataset = 'Dataset3';

file_extension = '_s%03d.png';

working_dir = 'C:/Path/to/VAST_exports/projects_mip1_exports/';

segfiletemplate = strcat(working_dir, dataset, '_segmentation/', file_extension);
segmentdatafilename = strcat(working_dir, dataset, '_metadata.txt');
excludemaskfiletemplate = strcat(working_dir, dataset, '_artifacts/', file_extension);
cbtagfiletemplate = strcat(working_dir, dataset, '_soma_seeds/', file_extension);
filledfiletemplate = strcat(working_dir, dataset, '_segmentation_expanded/', file_extension);
paddedfilledfiletemplate = strcat(working_dir, dataset, '_segmentation_expanded_padded/', file_extension);

segmentsmaskreduce=10;
segmentsmaskexpand=70;
includedsegmentsdir='Cells';

taggedsegoffset=1; %Tagged region IDs will be shifted up by this amount, which has to be larger than the number of segments in the segmentation

if (strcmp(dataset, 'Dataset1'))
  sections = 0:299;
end
if (strcmp(dataset, 'Dataset2'))
  sections = 0:363;
end
if (strcmp(dataset, 'Dataset3'))
  sections = 0:263;
end
if (strcmp(dataset, 'Dataset4'))
  sections = 0:299;
end
if (strcmp(dataset, 'Dataset5'))
  sections = 0:857;
end
if (strcmp(dataset, 'Dataset6'))
  sections = 0:426;
end
if (strcmp(dataset, 'Dataset8'))
  sections = 0:699;
end
if (strcmp(dataset, 'test'))
  sections = [76 167 203];
end



useiterationlimit=0; %If this flag is 1, filling will be stopped after maxiterations iterations
maxiterations=64;

processfilling=1; %If this flag is 1, it will compute the filling, if not, it will just reload the images from filledfiletemplate

downsample=1; %2 means downsampling by 2 in both directions


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Generate target folder if it doesn't exist
d=dir(fileparts(filledfiletemplate));
if (min(size(d))==0)
  mkdir(fileparts(filledfiletemplate));
end;
d=dir(fileparts(paddedfilledfiletemplate));
if (min(size(d))==0)
  mkdir(fileparts(paddedfilledfiletemplate));
end;

%Load segmentation metadata and find segments to use for analysis
disp(['Loading ' segmentdatafilename ' ...']);
[names,data,firstelement]=scanvastcolorfile(segmentdatafilename,0);


if (taggedsegoffset<size(data,1))
  disp(['Warning: taggedsegoffset too small, according to ' segmentdatafilename '. Corrected.']);
  taggedsegoffset=size(data,1)+1;
end;



fid=getidsfromexactname(names,includedsegmentsdir);
cids=getchildtreeids(data,fid);

% Make dict from segment id to valid id inside dir.
segidx_to_valididx=zeros(size(data,1),1);
n = 0;
for i=1:length(segidx_to_valididx)
  if (ismember(i, cids))
    n = n + 1;
    segidx_to_valididx(i) = n;
  else
    segidx_to_valididx(i) = -1;
  end;
end;


touchareamatrix_all=zeros(length(cids),length(cids));
touchareamatrix_nocb=zeros(length(cids),length(cids));
touchareamatrix_nocbcb=zeros(length(cids),length(cids));

%nocbcb_coords = cell(length(cids),length(cids));

surface_area_matrix=zeros(range(sections)+1,length(cids));
%distance_from_center=zeros(range(sections)+1,length(cids));

matrixidx_to_id = zeros(length(cids), 1);



for i=sections
  %Load segmentation of this section
  segfilename=sprintf(segfiletemplate,i);
  disp(['Loading ' segfilename ' ...']);
  seg=rgbdecode(imread(segfilename));


  includemask=[];

  %Load artifacts mask of this section from file
  maskfilename=sprintf(excludemaskfiletemplate,i);
  disp(['Loading ' maskfilename ' ...']);
  artifactsmask=imread(maskfilename);

  %Determine cell bodies
  shiftcids=[];
  cbtagfilename=sprintf(cbtagfiletemplate,i);
  disp(['Loading ' cbtagfilename ' ...']);
  cbtagimg=imread(cbtagfilename);
  if (sum(cbtagimg(:))>0)

    %Do connected components analysis of tag image
    cc=bwconncomp(cbtagimg);
    cbtagcoords=zeros(length(cc.PixelIdxList),3); %[[x, y, id], ..]
    for c=1:length(cc.PixelIdxList)
      p=cc.PixelIdxList{c}; %linear indices for tagged pixels of connected component c
      [y,x]=ind2sub(size(cbtagimg),p);
      cbtagcoords(c,1)=floor(mean(y));
      cbtagcoords(c,2)=floor(mean(x));
      cbtagcoords(c,3)=seg(cbtagcoords(c,1),cbtagcoords(c,2));

      % Skip tags not marking segments included in proper segments folder
      if (~ismember(cbtagcoords(c,3), cids))
        continue;
      end;
      % Skip tags not marking anything.
      if (cbtagcoords(c,3)==0)
        continue;
      end;

      % For each tag, do connected component analysis of the segment it
      % overlaps with.
      lidx=sub2ind(size(cbtagimg),cbtagcoords(c,1),cbtagcoords(c,2));

      csegimg=(seg==cbtagcoords(c,3));
      scc=bwconncomp(csegimg);
      if (length(scc.PixelIdxList)==1)
        %Only one area - renumber all pixels of this segment
        if (seg(lidx)<taggedsegoffset) %we need this because some connected regions have two tags!
          seg(scc.PixelIdxList{1})=seg(scc.PixelIdxList{1})+taggedsegoffset;
          shiftcids=[shiftcids cbtagcoords(c,3)+taggedsegoffset];
        end;
      else
        %Check which area the tag is in
        n=0;
        for s=1:length(scc.PixelIdxList)
          if (ismember(lidx,scc.PixelIdxList{s}))
            if (seg(lidx)<taggedsegoffset) %we need this because some connected regions have two tags!
              seg(scc.PixelIdxList{s})=seg(scc.PixelIdxList{s})+taggedsegoffset;
              shiftcids=[shiftcids cbtagcoords(c,3)+taggedsegoffset];
              n=n+1;
            end;
          end;
        end;
        if (n~=1)
          disp(sprintf('WARNING: Cell body tag problem in section %d, segment %d: %d regions tagged',i,cbtagcoords(c,3),n));
        end;
      end;
    end;
  end;


  %Compute a mask which allows only the area close to the included segments
  mask=ismember(seg,[cids shiftcids]);
  mask=bwareaopen(mask,segmentsmaskreduce);
  d=bwdist(mask);
  mask=(d<segmentsmaskexpand);
  d=bwdist(~mask);
  mask(d<segmentsmaskexpand)=0;


  mask(artifactsmask>0)=0;

  %Remove segmentation outside of mask
  disp('Masking ...');
  seg(mask==0)=0;

  origsize=size(seg);
  if (downsample~=1)
    seg=seg(1:downsample:end,1:downsample:end);
    mask=mask(1:downsample:end,1:downsample:end);
  end;
  dssize=size(seg);

  %Crop to bounding box
  cxmin=find(max(seg,[],1)>0,1,'first');
  cxmax=find(max(seg,[],1)>0,1,'last');
  cymin=find(max(seg,[],2)>0,1,'first');
  cymax=find(max(seg,[],2)>0,1,'last');

  if (min(size(cxmin))==0)
    disp(sprintf('Section %d is empty. Skipping.',i));
    continue;
  end;

  disp('Filling ...');
  mask=mask(cymin:cymax,cxmin:cxmax);
  figure(1);
  subplot(2,2,1);
  imagesc(mask);
  axis equal;
  title(sprintf('Section %d Mask',i));

  seg=seg(cymin:cymax,cxmin:cxmax);
  figure(1);
  subplot(2,2,2);
  imagesc(seg);
  axis equal;
  title(sprintf('%d Masked Seg',i));

  %fill all gaps between segments within the mask
  filledfilename=sprintf(filledfiletemplate,i);
  if (processfilling==0)
    disp(['Loading ' filledfilename ' ...']);
    fseg=rgbdecode(imread(filledfilename));
  else

    if (useiterationlimit==0)
      fseg=fillmasked(seg,mask);
    else
      fseg=fillmaxdist(seg,mask,maxiterations);
    end;


    pseg=zeros(dssize,'uint32');
    pseg(cymin:cymax,cxmin:cxmax)=fseg;

    if (downsample~=1)
      sfseg=imresize(pseg,origsize,'nearest');
    else
      sfseg=pseg;
    end;

    sfseg(sfseg>taggedsegoffset)=sfseg(sfseg>taggedsegoffset)-taggedsegoffset; %Undo cellbody tagging for padded image stack
    rgbseg=rgbencode(sfseg);
    paddedfilledfilename=sprintf(paddedfilledfiletemplate,i);
    disp(['Writing ' paddedfilledfilename ' ...']);
    imwrite(rgbseg,paddedfilledfilename);

    %Save fseg for reloading if processfilling is 0
    rgbseg=rgbencode(fseg);
    disp(['Writing ' filledfilename ' ...']);
    imwrite(rgbseg,filledfilename);
  end;
  figure(1);
  subplot(2,2,3);

  %remove segments defined as outside folder.
  fseg(ismember(fseg,[cids shiftcids])==0) = 0;

  imagesc(fseg);
  axis equal;
  title(sprintf('%d seg expanded',i));

  %count number of neighbor voxels
  segnr=unique(fseg);
  segnr(segnr==0)=[];
  for j=1:length(segnr)
    id=segnr(j);

    iseg=(fseg==id); %mask out ID region in filled segmentation image

    %find bounding box
    xmin=find(max(iseg,[],1)>0,1,'first');
    xmax=find(max(iseg,[],1)>0,1,'last');
    ymin=find(max(iseg,[],2)>0,1,'first');
    ymax=find(max(iseg,[],2)>0,1,'last');

    %expand bounding box by 1 with border control
    xmin=max([xmin-1 1]);
    xmax=min([xmax+1 size(iseg,2)]);
    ymin=max([ymin-1 1]);
    ymax=min([ymax+1 size(iseg,1)]);

    %Crop to region for faster analysis
    ciseg=iseg(ymin:ymax,xmin:xmax);
    cmask=mask(ymin:ymax,xmin:xmax);
    cfseg=fseg(ymin:ymax,xmin:xmax);

    oiseg=ciseg; %store previous region for comparison

    %expand ID region by 1
    ciseg(1:end-1,:)=max(ciseg(1:end-1,:),ciseg(2:end,:));
    ciseg(2:end,:)=max(ciseg(2:end,:),ciseg(1:end-1,:));
    ciseg(:,1:end-1)=max(ciseg(:,1:end-1),ciseg(:,2:end));
    ciseg(:,2:end)=max(ciseg(:,2:end),ciseg(:,1:end-1));

    ciseg(cmask==0)=0; %constrain to mask
    ciseg(oiseg>0)=0; %leave only border
    m=cfseg.*uint32(ciseg);

    subplot(2,2,4);
    imagesc(m);
    [seg,num]=count_unique(m);
    num(seg==0)=[];
    seg(seg==0)=[];

    for k=1:length(seg)
      if (id<taggedsegoffset)
        id1 = segidx_to_valididx(id);
        id1_iscb = 0;
      else
        id1 = segidx_to_valididx(id-taggedsegoffset);
        id1_iscb = 1;
      end;
      if (seg(k)<taggedsegoffset)
        id2 = segidx_to_valididx(seg(k));
        id2_iscb = 0;
      else
        id2 = segidx_to_valididx(seg(k)-taggedsegoffset);
        id2_iscb = 1;
      end;
      touchareamatrix_all(id1,id2)=touchareamatrix_all(id1,id2)+num(k);
      if (~id1_iscb || ~id2_iscb)
        touchareamatrix_nocbcb(id1,id2)=touchareamatrix_nocbcb(id1,id2)+num(k);
        %[row, col] = find(m == seg(k));
        %coords = [(row+cymin+ymin)*downsample (col+cxmin+xmin)*downsample i*ones(size(row, 1), 1)];
        %nocbcb_coords{id1,id2} = [nocbcb_coords{id1,id2}; coords];
      end;
      if (~id1_iscb && ~id2_iscb)
        touchareamatrix_nocb(id1,id2)=touchareamatrix_nocb(id1,id2)+num(k);
      end;
    end;

    %get area of segment on section
    if (id<taggedsegoffset)
      surface_area = sum(iseg(:));
      idx = segidx_to_valididx(id);
      surface_area_matrix(i-min(sections)+1, idx) = surface_area;
    end;


  end;

  pause(0.1);
end;





valid_names = strrep(strtok(names(segidx_to_valididx~=-1)), '-', '_');

%
% % Join dublicate names.
% for i=1:length(valid_names)
%
%     name = valid_names{i};
%     name_idx = find(strcmp(valid_names, name));
%     if length(name_idx) == 1
%         continue;
%     end;
%     initial_i = name_idx(1);
%     if i == initial_i
%         continue;
%     end;
%
%     disp(strcat('Duplicate cell:', {' '}, name));
%     valid_names{i} = [name num2str(2, '%d')];
%
%     for j=1:length(valid_names)
%
%         touchareamatrix_nocbcb(initial_i, j) = touchareamatrix_nocbcb(initial_i, j) + touchareamatrix_nocbcb(i, j);
%         touchareamatrix_nocbcb(j, initial_i) = touchareamatrix_nocbcb(j, initial_i) + touchareamatrix_nocbcb(j, i);
%
%         touchareamatrix_nocbcb(i, j) = 0;
%         touchareamatrix_nocbcb(j, i) = 0;
%
%         nocbcb_coords{initial_i, j} = [nocbcb_coords{initial_i, j}; nocbcb_coords{i, j}];
%         nocbcb_coords{j, initial_i} = [nocbcb_coords{j, initial_i}; nocbcb_coords{j, i}];
%
%         nocbcb_coords{i, j} = {};
%         nocbcb_coords{j, i} = {};
%
%     end;
%
% end;
%




% save('touchareamatrix_nocbcb.mat','touchareamatrix_nocbcb');


writetable(array2table(touchareamatrix_nocbcb, 'VariableNames', valid_names), strcat(working_dir, dataset, '_adjacency.csv'));
%writetable(array2table(touchareamatrix_nocb, 'VariableNames', valid_names), strcat(working_dir, 'touchareamatrix_nocb_center.csv'));
%writetable(array2table(touchareamatrix_all, 'VariableNames', valid_names), strcat(working_dir, 'touchareamatrix_all_center.csv'));


%disp('saving mat');
%save('myfile.mat', 'nocbcb_coords');

% disp('making json');
% json = jsonencode(nocbcb_coords);


%fid = fopen(strcat(working_dir, dataset, '_adjacency_coords.json'), 'wt');
%split_size = 10000000;
%json_length = length(json);
%for i=1:split_size:json_length
%    fprintf(fid, '%s\n', json(i:min(i+split_size-1, json_length)));
%end;
%fclose(fid);

writetable(array2table(surface_area_matrix, 'VariableNames', valid_names), strcat(working_dir, dataset, '_volume_area.csv'));
%writetable(array2table(distance_from_center, 'VariableNames', valid_names), strcat(working_dir, dataset, '_distance_from_center.csv'));
