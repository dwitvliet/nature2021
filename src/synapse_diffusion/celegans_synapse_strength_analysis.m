%celegans_synapse_strength_analysis.m
%By Daniel Berger, November 08, 2020

%Instructions:
%- Open VAST, enable API
%- Run VastTools.m in Matlab (newest version)
%- Adjust parameters below
%- Run this script.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Switch these flags to 1 for the parts of the script you want to run.
flags.loadlayers = 1; %Run this only once to load up the layers
flags.checkboundingboxes = 1; %Some analysis of bounding boxes and anchor points
flags.computediffusion = 1; %The actual neurotransmitter diffusion simulation

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% These are the parameters for data sources and simulation
param.basefolder='C:\Path\to\VAST_exports\projects\';

param.emfilename='Dataset8.vsv';
param.segmentationfilename='Dataset8_segmentation_withsoma_boundary.vss';
param.synapsefilename='Dataset8_synapses.vss';

param.miplevel=1; %at least in the test data, the resolution at mip1 is 4x4x30 nm (~ 8x lower in z)
%param.volsize_pix=[512,512,64];
param.volboundary_mu=[1,1,1]; %Padding around bounding box for diffusion region
param.synexpandradiusxy_pix=6; %Dilation of the synapse label when it is painted inside the source object to get the emitter surface
param.nrofparticles=1000; %Number of neurotransmitter molecules to simulate
param.diffusionsteplength_nm=4; %For the diffusion simulation, how far does the particle travel each step
param.nrdiffusionvectors=10000; %Number of precomputed diffusion vectors to sample from randomly during the simulation
param.maxiterations=10000; %End simulation after this many diffusion steps. Prevents endless loop in case of holes in segmentation
param.targetfilename='Dataset8_syndiffusion.json';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global vdata;
if (min(size(vdata))==0)
  disp('WARNING: Not connected to VAST. Please first run VAST, enable the API, then VastTools.m and connect.');
  return;
end;

%Disconnect/Reconnect (VAST may have been restarted)
if (vdata.state.isconnected)
  vdata.ui.menu.connectvast.MenuSelectedFcn(); %Disconnect
  vdata.ui.menu.connectvast.MenuSelectedFcn(); %Connect
else
  vdata.ui.menu.connectvast.MenuSelectedFcn(); %Connect
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (flags.loadlayers==1)
  %[emlay,res]=vdata.vast.loadlayer([flags.basefolder 'human_EM_h01_goog8_zmip_bicubic_with_overview.vsvr']);
  [emlay,res]=vdata.vast.loadlayer([param.basefolder param.emfilename]);
  [seglay,res]=vdata.vast.loadlayer([param.basefolder param.segmentationfilename]);
  [synlay,res]=vdata.vast.loadlayer([param.basefolder param.synapsefilename]);

  info=vdata.vast.getinfo();
  mipscalefactors=vdata.vast.getmipmapscalefactors(seglay);

  vdata.vast.setselectedlayernr(synlay);
  syndata=vdata.vast.getallsegmentdatamatrix();
  fids=vdata.vast.getimmediatechildids(0);
  synids=vdata.vast.getimmediatechildids(fids);

  disp(sprintf('%d synapses found.',length(synids)));
end;

if (flags.checkboundingboxes==1)
  bbvol=zeros(length(synids),1);

  for cnr=1:length(synids)
    id=synids(cnr);
    vanchor=syndata(id,11:13);
    vbbox=syndata(id,19:24);
    if ((vanchor(1)>=vbbox(1))&&(vanchor(1)<=vbbox(4))&&(vanchor(2)>=vbbox(2))&&(vanchor(2)<=vbbox(5))&&(vanchor(3)>=vbbox(3))&&(vanchor(3)<=vbbox(6)))
    else
      disp('WARNING: id %d: Anchor point outside of bounding box!',id);
    end;
    minx=vbbox(1); maxx=vbbox(4);
    miny=vbbox(2); maxy=vbbox(5);
    minz=vbbox(3); maxz=vbbox(6);
    bbvol(cnr)=(maxx-minx+1)*(maxy-miny+1)*(maxz-minz+1);
%     if (((maxx-minx)>50*250)||((maxy-miny)>50*250))
%       disp(sprintf('WARNING: Segment %d has a large bounding box: (%d..%d, %d..%d, %d..%d)!',id,minx,maxx,miny,maxy,minz,maxz));
%     end;
  end;
  figure(10);
  plot(bbvol);
  grid on;
  title([param.synapsefilename sprintf(' | %d synapses',length(synids))]);
  xlabel('Synapse Nr.');
  ylabel('Bounding box volume (pixels^3)');
end;


if (flags.computediffusion==1)
  errcount=0;
  msf=[1 1 1];
  if (param.miplevel>0)
    msf=mipscalefactors(param.miplevel,:);
  end;
  xextend=ceil((param.volboundary_mu(1)*1000)/(info.voxelsizex*msf(1)));
  yextend=ceil((param.volboundary_mu(2)*1000)/(info.voxelsizey*msf(2)));
  zextend=ceil((param.volboundary_mu(3)*1000)/(info.voxelsizez*msf(3)));

  mipsize=floor([info.datasizex/msf(1) info.datasizey/msf(2) info.datasizez/msf(3)]);


  %Precompute diffusion vectors
    %param.diffusionsteplength_nm=4; %For the diffusion simulation, how far does the particle travel each step
    %param.nrdiffusionvectors=1000; %Number of precomputed diffusion vectors to sample from randomly during the simulation
  %dvec1 = fibonacci_spiral_sphere(param.nrdiffusionvectors);
  dvec = fibonacci_spiral_sphere2(param.nrdiffusionvectors);
%   figure;
%   plot3(dvec(:,1),dvec(:,2),dvec(:,3),'r.');
%   hold on;
%   plot3(dvec1(:,1),dvec1(:,2),dvec1(:,3),'.');
%   hold off;
%   grid on;

  dsx=param.diffusionsteplength_nm/(info.voxelsizex*msf(1));
  dsy=param.diffusionsteplength_nm/(info.voxelsizey*msf(2));
  dsz=param.diffusionsteplength_nm/(info.voxelsizez*msf(3));
  dvec(:,1)=dvec(:,1)*dsx;
  dvec(:,2)=dvec(:,2)*dsy;
  dvec(:,3)=dvec(:,3)*dsz;
%     figure;
%     plot3(dvec(:,1),dvec(:,2),dvec(:,3),'.');
%     grid on;

  synvolpersection=zeros(info.datasizez,length(synids));


  res = vdata.vast.setapilayersenabled(1);
  for cnr=1:length(synids)
    id=synids(cnr);
    %ids_todo = [571 591 1040 1563 1578 1611 1892 2183 2316 2855 2885 3566];
    %if ~ismember(id, ids_todo)
    %    continue;
    %end;
    disp(sprintf('Synapse %d of %d...',cnr,length(synids)));

    vanchor=syndata(id,11:13);
    vbbox=syndata(id,19:24);

    syn_se=strel('disk',param.synexpandradiusxy_pix);
    seg_se=strel('sphere',1);

    %Define bounding box for analysis at miplevel
    ebbox=[floor(vbbox(1)/msf(1)) floor(vbbox(2)/msf(2)) floor(vbbox(3)/msf(3)) ceil(vbbox(4)/msf(1)) ceil(vbbox(5)/msf(2)) ceil(vbbox(6)/msf(3))];
    ebbox=ebbox+[-xextend -yextend -zextend xextend yextend zextend];
    %Clip to volume
    if (ebbox(1)<0) ebbox(1)=0; end;
    if (ebbox(2)<0) ebbox(2)=0; end;
    if (ebbox(3)<0) ebbox(3)=0; end;
    if (ebbox(4)>=mipsize(1)) ebbox(4)=mipsize(1)-1; end;
    if (ebbox(5)>=mipsize(2)) ebbox(5)=mipsize(2)-1; end;
    if (ebbox(6)>=mipsize(3)) ebbox(6)=mipsize(3)-1; end;

    %Load segmentation in bounding box
    %vdata.vast.setselectedlayernr(seglay);
    res = vdata.vast.setselectedapilayernr(seglay);
    vdata.vast.setsegtranslation([],[]);
    segvol=vdata.vast.getsegimageRLEdecoded(param.miplevel, ebbox(1),ebbox(4),ebbox(2),ebbox(5),ebbox(3),ebbox(6), 0,1);

    %Load synapse annotation in bounding box
    res = vdata.vast.setselectedapilayernr(synlay);
    res = vdata.vast.setsegtranslation(id, 1);
    synvol=vdata.vast.getsegimageRLEdecoded(param.miplevel, ebbox(1),ebbox(4),ebbox(2),ebbox(5),ebbox(3),ebbox(6), 0,1);

    syn{cnr}.vast_id=id;
    syn{cnr}.nrsynapsevoxels=sum(synvol(:));

    %Store how many synapse pixels of this synapse are in each section
    svps=squeeze(sum(sum(synvol,1),2)); %synapse voxels per section in loaded block
    if (ebbox(3)==0)
      if (svps(1)~=0)
        disp(sprintf("  WARNING: Synapse %d has %d pixels in section 0, not stored in synvolpersection!",cnr,svps(1)));
        synvolpersection(ebbox(3)+1:ebbox(6),cnr)=svps(2:end);
      end;
    else
      synvolpersection(ebbox(3):ebbox(6),cnr)=svps; %here we assume that ebbox(3) is never 0
    end;

    %Make list of section-synapse voxel counts for this synapse; first column is slice nr  (counting from 0), second is synapse voxel count
    snr=find(svps);
    syn{cnr}.synapseslicevoxels=[snr+ebbox(3)-1 svps(snr)];

    %Compute surface voxels for neurotransmitter emission
    s=(synvol==1).*(segvol~=0); %remove voxels outside segment
    c=segvol(s==1);
    if (min(size(c))==0)
      if (sum(synvol(:))==0)
        disp(sprintf('ERROR: Synapse with id %d doesn''t exist (not painted)!',id));
        errcount=errcount+1;
      else
        disp(sprintf('ERROR: Synapse with id %d doesn''t exist (painted region not on any object)!',id));
        errcount=errcount+1;
      end;
    else
      [val,num]=count_unique(c);
      preid=val(1);

      if (max(size(val))>1)
        vn=sortrows([val num],-2);
        preid=vn(1,1);
      end;

      %add synapse label to presynaptic segment (D. Witvliet suggested)
      segvol(find(synvol))=preid; %synvol is 1 where the synapse is labeled, 0 everywhere else

      %Expand synapse label in XY (assuming it's not necessary in z because the resolution is much lower)
      esyn=imdilate(synvol,syn_se); %expanded synapse label
      %Find empty voxels next to preid objects which overlap with s
      pseg=(segvol==preid); %Compute boundary voxels around preid
      epseg=imdilate(pseg,seg_se);
      %epseg(pseg==1)=0;
      epseg(segvol~=0)=0; %restrict to extracellular pixels
      epseg(esyn==0)=0; %restrict to synapse

      emittervoxels=find(epseg);
      syn{cnr}.nremittervoxels=length(emittervoxels);
      if (length(emittervoxels)==0)
        disp(sprintf('ERROR: Synapse with id %d has no emitter voxels!',id));
        errcount=errcount+1;
      else

        %Pick random starting locations for neurotransmitter particles from emitter voxels
        sourcevoxnr=randi(length(emittervoxels),[param.nrofparticles,1]);
        sourcevox=emittervoxels(sourcevoxnr);
        [y,x,z]=ind2sub(size(epseg),sourcevox);
        nt=[y x z zeros(length(x),1) ones(length(x),1)]; %XYZ coords, target ID, status
        %figure(11);
        %       plot3(nt(:,1),nt(:,2),nt(:,3),'.');
        %       grid on;


        %      for i=1:1000
        %         vnr=randi(param.nrdiffusionvectors,[param.nrofparticles,1]);
        %         v=dvec(vnr,:);
        %         nt(:,1:3)=nt(:,1:3)+v;
        %
        %         p=round(nt(:,1:3));
        %         vids=segvol(sub2ind(size(segvol),p(:,1),p(:,2),p(:,3))); %extract segment ids where neurotransmitter voxel landed
        count=0;

        while ((count<param.maxiterations)&&(sum(nt(:,5))>0)) %while there are still particles running
          vnr=randi(param.nrdiffusionvectors,[param.nrofparticles,1]);
          v=dvec(vnr,:).*[nt(:,5) nt(:,5) nt(:,5)]; %only move particles which have not settled (5th column is 1)
          nt(:,1:3)=nt(:,1:3)+v;

          p=round(nt(:,1:3));

          minp=min(p); maxp=max(p);
          if ((min(minp)<1)||(max(maxp-size(synvol))>0))
            %Some particles exited the volume. identify and settle at 0
            exited=[find(p(:,1)<1) find(p(:,2)<1) find(p(:,3)<1) find(p(:,1)>size(segvol,1)) find(p(:,2)>size(segvol,2)) find(p(:,3)>size(segvol,3))];
            nt(exited,1:3)=nt(exited,1:3)-v(exited,:); %move back
            nt(exited,5)=0; %flag as settled
            p=round(nt(:,1:3));
          end;

          vids=segvol(sub2ind(size(segvol),p(:,1),p(:,2),p(:,3))); %extract segment ids where neurotransmitter voxel landed

          nt(vids==preid,1:3)=nt(vids==preid,1:3)-v(vids==preid,:); %particles which would have traveled into the presynaptic side are put back
          vids(vids==preid)=0;

          %find and settle particles which ended up in a neighbor segment
          sp=find(vids>0); %settling particles
          nt(sp,4)=vids(sp);
          nt(sp,5)=0;

          count=count+1;
          if (mod(count,1000)==0)
            %Visualize diffusion process
            settled=find(nt(:,5)==0);
            active=find(nt(:,5)==1);
            plot3(nt(settled,1),nt(settled,2),nt(settled,3),'.');
            hold on;
            plot3(nt(active,1),nt(active,2),nt(active,3),'r.');
            hold off;
            pause(0.01);
          end;
        end;

        [val,num]=count_unique(nt(:,4));

        syn{cnr}.preid=preid;
        syn{cnr}.postid=val;
        syn{cnr}.poststrength=num;
        syn{cnr}.iterationsused=count;

        %figure(11);
        plot3(y,x,z,'r.'); %show all (potential) emitter voxels in red
        hold on;
        settled=find(nt(:,5)==0);
        plot3(nt(settled,1),nt(settled,2),nt(settled,3),'.'); %show particle target locations in blue
        hold off;
        grid on;
        title(sprintf('Synapse %d/%d, %d iterations',cnr,length(synids),count));
        pause(0.01);
      end;
    end;
  end;

  vdata.vast.setsegtranslation([],[]);
  res = vdata.vast.setapilayersenabled(0);

  txt=jsonencode(syn);
  fileID = fopen(param.targetfilename,'w');
  fprintf(fileID,txt);
  fclose(fileID);
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% https://dl.acm.org/doi/10.1145/2816795.2818131
% See: https://github.com/bduvenhage/Bits-O-Cpp
% std::vector<Vec3> fibonacci_spiral_sphere(const int num_points) {
%     std::vector<Vec3> vectors;
%     vectors.reserve(num_points);
%     const double gr=(sqrt(5.0) + 1.0) / 2.0;  // golden ratio = 1.6180339887498948482
%     const double ga=(2.0 - gr) * (2.0*M_PI);  // golden angle = 2.39996322972865332
%
%     for (size_t i=1; i <= num_points; ++i) {
%         const double lat = asin(-1.0 + 2.0 * double(i) / (num_points+1));
%         const double lon = ga * i;
%         const double x = cos(lon)*cos(lat);
%         const double y = sin(lon)*cos(lat);
%         const double z = sin(lat);
%         vectors.emplace_back(x, y, z);
%     }
%     return vectors;
% }

function vec = fibonacci_spiral_sphere(num_points)
  vec=zeros(num_points,3);
  gr=(sqrt(5.0) + 1.0) / 2.0;  % golden ratio = 1.6180339887498948482
  ga=(2.0 - gr) * (2.0*pi);    % golden angle = 2.39996322972865332
  for i=1:num_points
    lat = asin(-1.0 + 2.0 * i / (num_points+1));
    lon = ga * i;
    x = cos(lon)*cos(lat);
    y = sin(lon)*cos(lat);
    z = sin(lat);
    vec(i,:)=[x y z];
  end
end

%float3 SF(float i, float n) {
%  //i counts from 0..n-1
%  float u = 2.0*PI*PHI*i;
%  float v = 1.0 - (2.0*i + 1.0)/n;
%  float sinTheta = sqrt(1.0 - v*v);
%  return float3(cos(u)*sinTheta, sin(u)*sinTheta, v);
%}
function vec = fibonacci_spiral_sphere2(num_points)
  vec=zeros(num_points,3);
  gr=(sqrt(5.0) + 1.0) / 2.0;
  for i=0:num_points-1
    u = 2.0*pi*gr*i;
    v = 1.0 - (2.0*i + 1.0)/num_points;
    sinTheta = sqrt(1.0 - v*v);
    vec(i+1,:)=[cos(u)*sinTheta, sin(u)*sinTheta, v];
  end;
end
