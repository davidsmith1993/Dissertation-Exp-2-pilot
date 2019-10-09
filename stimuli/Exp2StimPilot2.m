 function Exp2StimPilot2
%
%generates stimuli for David Exp 2 pilot 1 9/5/2014
%using spatial frequency of large circular grating (at 90 deg) and angle of displaced line
%   need dimensions that are spatially separate so that one dimension can be omitted during inference task
%       this will hurt accuracy (based on relational experiment, but Ss should be able to learn)
%   need them to process both dimensions (i.e., shouldn't be able to ignore a dimension and learn)
%       grating should be at fixed position relative to line varying in angle
%   **try spatially centered grating (maybe some trial/trial jitter) with radial line that can vary from 0-270 degrees
%09/05/18 DBS

close all
%figure stuff
fontname = 'ArialBold'; 
axisLabelSize = 16;
titleSize = 18;
textSize = 14;
linewidth = 2;
markerSize = 5;

%Using this function to generate the uniform stimuli
%function [minMax,X] = GenUnifSample(range,n,cat)

nPerCat = 40;
nCat = 2;
nTrainBlocks = 4;
trialsPerBlock = nPerCat*nCat;


Arangelow = [25 35; 0 100];
Arangemid = [35 45; 0 100];
Arangehigh = [45 55; 0 100];
Brangelow = [60 70; 0 100];
Brangemid = [70 80; 0 100];
Brangehigh = [80 90; 0 100];

%%%%%%%%%%%%Create stim for consistent condition
%Scalein + scaleout*2 needs to equal 1
scaleout = 6/40;
scalein = 28/40;

for b = 1:nTrainBlocks
    %cols are: cat, diameter, orientation, block,  ...
    %stimlus_location (used to matter for 2AFC version of experiment, but just a filler for this version included to minimize changes to experiment program)
    %range(1=low,2=high) ...
    %dimension_to_infer(11=sf, 22=or)
    %A/notA or B/notB (only meaningful for training stimuli)]


    [minMax,Alow] = GenUnifSample(Arangelow,nPerCat*scaleout,1);
    catAlow = [Alow b*ones(nPerCat*scaleout,1) ...
               9999*ones(nPerCat*scaleout,1) ...
               ones(nPerCat*scaleout,1) ...
               randrows([11*ones(nPerCat*scaleout/2,1); 22*ones(nPerCat*scaleout/2,1)]) ...
               randrows([ones(nPerCat*scaleout/2,1) ;2*ones(nPerCat*scaleout/2,1)])];
    [minMax,Amid] = GenUnifSample(Arangemid,nPerCat*scalein,1);
    catAmid = [Amid b*ones(nPerCat*scalein,1) ...
               9999*ones(nPerCat*scalein,1) ...
               ones(nPerCat*scalein,1) ...
               randrows([11*ones(nPerCat*scalein/2,1); 22*ones(nPerCat*scalein/2,1)]) ...
               randrows([ones(nPerCat*scalein/2,1) ;2*ones(nPerCat*scalein/2,1)])];
              
    [minMax,Ahigh] = GenUnifSample(Arangehigh,nPerCat*scaleout,1);               
    catAhigh = [Ahigh b*ones(nPerCat*scaleout,1) ...
               9999*ones(nPerCat*scaleout,1) ...
               ones(nPerCat*scaleout,1) ...
               randrows([11*ones(nPerCat*scaleout/2,1); 22*ones(nPerCat*scaleout/2,1)]) ...
               randrows([ones(nPerCat*scaleout/2,1) ;2*ones(nPerCat*scaleout/2,1)])];
    [minMax,Blow] = GenUnifSample(Brangelow,nPerCat*scaleout,2);
    catBlow = [Blow b*ones(nPerCat*scaleout,1) ...
               9999*ones(nPerCat*scaleout,1) ...
               ones(nPerCat*scaleout,1) ...
               randrows([11*ones(nPerCat*scaleout/2,1); 22*ones(nPerCat*scaleout/2,1)]) ...
               randrows([ones(nPerCat*scaleout/2,1) ;2*ones(nPerCat*scaleout/2,1)])];
    [minMax,Bmid] = GenUnifSample(Brangemid,nPerCat*scalein,2);
    catBmid = [Bmid b*ones(nPerCat*scalein,1) ...
               9999*ones(nPerCat*scalein,1) ...
               ones(nPerCat*scalein,1) ...
               randrows([11*ones(nPerCat*scalein/2,1); 22*ones(nPerCat*scalein/2,1)]) ...
               randrows([ones(nPerCat*scalein/2,1) ;2*ones(nPerCat*scalein/2,1)])];
              
    [minMax,Bhigh] = GenUnifSample(Brangehigh,nPerCat*scaleout,2);               
    catBhigh = [Bhigh b*ones(nPerCat*scaleout,1) ...
               9999*ones(nPerCat*scaleout,1) ...
               ones(nPerCat*scaleout,1) ...
               randrows([11*ones(nPerCat*scaleout/2,1); 22*ones(nPerCat*scaleout/2,1)]) ...
               randrows([ones(nPerCat*scaleout/2,1) ;2*ones(nPerCat*scaleout/2,1)])];


RBcons(b*trialsPerBlock-trialsPerBlock+1:b*trialsPerBlock,:) = [catAlow;catAmid;catAhigh;catBlow;catBmid;catBhigh];
IIcons(b*trialsPerBlock-trialsPerBlock+1:b*trialsPerBlock,:) = [rotate2dDist(catAlow,.785398);rotate2dDist(catAmid,.785398);rotate2dDist(catAhigh,.785398);rotate2dDist(catBlow,.785398);rotate2dDist(catBmid,.785398);rotate2dDist(catBhigh,.785398)];
IIcons2 = [rotate2dDist(RBcons,.785)]
end

%%%%%%%%%%%%Create stim for inconsistent condition
%Scalein + scaleout*2 needs to equal 1
scaleout = 12/40;
scalein = 16/40;

for b = 1:nTrainBlocks
    %cols are: cat, diameter, orientation, block,  ...
    %stimlus_location (used to matter for 2AFC version of experiment, but just a filler for this version included to minimize changes to experiment program)
    %range(1=low,2=high) ...
    %dimension_to_infer(11=sf, 22=or)
    %A/notA or B/notB (only meaningful for training stimuli)]


    [minMax,Alow] = GenUnifSample(Arangelow,nPerCat*scaleout,1);
    catAlow = [Alow b*ones(nPerCat*scaleout,1) ...
               9999*ones(nPerCat*scaleout,1) ...
               ones(nPerCat*scaleout,1) ...
               randrows([11*ones(nPerCat*scaleout/2,1); 22*ones(nPerCat*scaleout/2,1)]) ...
               randrows([ones(nPerCat*scaleout/2,1) ;2*ones(nPerCat*scaleout/2,1)])];
    [minMax,Amid] = GenUnifSample(Arangemid,nPerCat*scalein,1);
    catAmid = [Amid b*ones(nPerCat*scalein,1) ...
               9999*ones(nPerCat*scalein,1) ...
               ones(nPerCat*scalein,1) ...
               randrows([11*ones(nPerCat*scalein/2,1); 22*ones(nPerCat*scalein/2,1)]) ...
               randrows([ones(nPerCat*scalein/2,1) ;2*ones(nPerCat*scalein/2,1)])];
              
    [minMax,Ahigh] = GenUnifSample(Arangehigh,nPerCat*scaleout,1);               
    catAhigh = [Ahigh b*ones(nPerCat*scaleout,1) ...
               9999*ones(nPerCat*scaleout,1) ...
               ones(nPerCat*scaleout,1) ...
               randrows([11*ones(nPerCat*scaleout/2,1); 22*ones(nPerCat*scaleout/2,1)]) ...
               randrows([ones(nPerCat*scaleout/2,1) ;2*ones(nPerCat*scaleout/2,1)])];
    [minMax,Blow] = GenUnifSample(Brangelow,nPerCat*scaleout,2);
    catBlow = [Blow b*ones(nPerCat*scaleout,1) ...
               9999*ones(nPerCat*scaleout,1) ...
               ones(nPerCat*scaleout,1) ...
               randrows([11*ones(nPerCat*scaleout/2,1); 22*ones(nPerCat*scaleout/2,1)]) ...
               randrows([ones(nPerCat*scaleout/2,1) ;2*ones(nPerCat*scaleout/2,1)])];
    [minMax,Bmid] = GenUnifSample(Brangemid,nPerCat*scalein,2);
    catBmid = [Bmid b*ones(nPerCat*scalein,1) ...
               9999*ones(nPerCat*scalein,1) ...
               ones(nPerCat*scalein,1) ...
               randrows([11*ones(nPerCat*scalein/2,1); 22*ones(nPerCat*scalein/2,1)]) ...
               randrows([ones(nPerCat*scalein/2,1) ;2*ones(nPerCat*scalein/2,1)])];
              
    [minMax,Bhigh] = GenUnifSample(Brangehigh,nPerCat*scaleout,2);               
    catBhigh = [Bhigh b*ones(nPerCat*scaleout,1) ...
               9999*ones(nPerCat*scaleout,1) ...
               ones(nPerCat*scaleout,1) ...
               randrows([11*ones(nPerCat*scaleout/2,1); 22*ones(nPerCat*scaleout/2,1)]) ...
               randrows([ones(nPerCat*scaleout/2,1) ;2*ones(nPerCat*scaleout/2,1)])];


RBincons(b*trialsPerBlock-trialsPerBlock+1:b*trialsPerBlock,:) = [catAlow;catAmid;catAhigh;catBlow;catBmid;catBhigh];
IIincons(b*trialsPerBlock-trialsPerBlock+1:b*trialsPerBlock,:) = [rotate2dDist(catAlow,.785398);rotate2dDist(catAmid,.785398);rotate2dDist(catAhigh,.785398);rotate2dDist(catBlow,.785398);rotate2dDist(catBmid,.785398);rotate2dDist(catBhigh,.785398)];
IIincons2 = [rotate2dDist(RBincons,.785)]
end


% %Generate Test stimuli
nTest = 80;
[minMax,TestA] = GenUnifSample(Arangemid,nTest/2,1);    
   TeststimA = [TestA 5*ones(nTest/2,1) ...
               9999*ones(nTest/2,1) ...
               ones(nTest/2,1) ...
               randrows([11*ones(nTest/4,1); 22*ones(nTest/4,1)]) ...
               randrows([ones(nTest/4,1) ;2*ones(nTest/4,1)])]; 
[minMax,TestB] = GenUnifSample(Brangemid,nTest/2,2);    
   TeststimB = [TestB 5*ones(nTest/2,1) ...
               9999*ones(nTest/2,1) ...
               ones(nTest/2,1) ...
               randrows([11*ones(nTest/4,1); 22*ones(nTest/4,1)]) ...
               randrows([ones(nTest/4,1) ;2*ones(nTest/4,1)])]; 
                      
TestRB = [TeststimA; TeststimB];           
%TestII = [rotate2dDist(TeststimA,.785398); rotate2dDist(TeststimB,.785398)]
TestII = [rotate2dDist(TestRB, .785398)];


% 
% %probe stimuli 7/27/2016
% Aprobemean = [650 100];
% Bprobemean = [500 250];
% Cprobemean = [350 100];
% Dprobemean = [500 -50];
% 
% t=60;
% nProbe = 32; %evenly divisible by nClustersprobe and nClustersprobe*2
% nClustersProbe = 4; %number of clusters of test stimuli-4 within training, 4 extrapolation
% 
% Aprobe =  [      ones(1,nProbe/nClustersProbe); ... %A test items = 1, B test items = 2
%                     linspace(Aprobemean(1),Aprobemean(1),nProbe/nClustersProbe);  ... x coords
%                     linspace(Aprobemean(2)-t,Aprobemean(2)+t,nProbe/nClustersProbe); ... %y coords
%                     (nTrainBlocks+1)*ones(1,nProbe/nClustersProbe); ... block number
%                     [randrows([-1*ones(nProbe/(2*nClustersProbe),1);ones(nProbe/(2*nClustersProbe),1)])]'; ... %target location for inference trials -1=left (filler)
%                     ones(1,nProbe/nClustersProbe); ... %in low range = 1; in high range = 2; in mid range = 3;
%                     [randrows([11*ones(nProbe/(2*nClustersProbe),1); 22*ones(nProbe/(2*nClustersProbe),1)])]'; ... %dimension 2 infer 11=SF, 22=OR
%                     9999*ones(1,nProbe/nClustersProbe)]; %filler
% Bprobe =  [      2*ones(1,nProbe/nClustersProbe); ... %A test items = 1, B test items = 2
%                     linspace(Bprobemean(1)+t,Bprobemean(1)-t,nProbe/nClustersProbe);  ... x coords
%                     linspace(Bprobemean(2),Bprobemean(2),nProbe/nClustersProbe); ... %y coords
%                     (nTrainBlocks+1)*ones(1,nProbe/nClustersProbe); ... block number
%                     [randrows([-1*ones(nProbe/(2*nClustersProbe),1);ones(nProbe/(2*nClustersProbe),1)])]'; ... %target location for inference trials -1=left (filler)
%                     ones(1,nProbe/nClustersProbe); ... %in low range = 1; in high range = 2; in mid range = 3;
%                     [randrows([11*ones(nProbe/(2*nClustersProbe),1); 22*ones(nProbe/(2*nClustersProbe),1)])]'; ... %dimension 2 infer 11=SF, 22=OR
%                     9999*ones(1,nProbe/nClustersProbe)]; %filler
% Cprobe =  [      3*ones(1,nProbe/nClustersProbe); ... %A test items = 1, B test items = 2
%                     linspace(Cprobemean(1),Cprobemean(1),nProbe/nClustersProbe);  ... x coords
%                     linspace(Cprobemean(2)-t,Cprobemean(2)+t,nProbe/nClustersProbe); ... %y coords
%                     (nTrainBlocks+1)*ones(1,nProbe/nClustersProbe); ... block number
%                     [randrows([-1*ones(nProbe/(2*nClustersProbe),1);ones(nProbe/(2*nClustersProbe),1)])]'; ... %target location for inference trials -1=left (filler)
%                     ones(1,nProbe/nClustersProbe); ... %in low range = 1; in high range = 2; in mid range = 3;
%                     [randrows([11*ones(nProbe/(2*nClustersProbe),1); 22*ones(nProbe/(2*nClustersProbe),1)])]'; ... %dimension 2 infer 11=SF, 22=OR
%                     9999*ones(1,nProbe/nClustersProbe)]; %filler
% Dprobe =  [      4*ones(1,nProbe/nClustersProbe); ... %A test items = 1, B test items = 2
%                     linspace(Dprobemean(1)+t,Dprobemean(1)-t,nProbe/nClustersProbe);  ... x coords
%                     linspace(Dprobemean(2),Dprobemean(2),nProbe/nClustersProbe); ... %y coords
%                     (nTrainBlocks+1)*ones(1,nProbe/nClustersProbe); ... block number
%                     [randrows([-1*ones(nProbe/(2*nClustersProbe),1);ones(nProbe/(2*nClustersProbe),1)])]'; ... %target location for inference trials -1=left (filler)
%                     ones(1,nProbe/nClustersProbe); ... %in low range = 1; in high range = 2; in mid range = 3;
%                     [randrows([11*ones(nProbe/(2*nClustersProbe),1); 22*ones(nProbe/(2*nClustersProbe),1)])]'; ... %dimension 2 infer 11=SF, 22=OR
%                     9999*ones(1,nProbe/nClustersProbe)]; %filler

%probe stimuli for RB task, block 4.
% RBprobes = [Aprobe(2:3,1:2:end) Bprobe(2:3,1:2:end) Cprobe(2:3,1:2:end) Dprobe(2:3,1:2:end)]';
% nProbes = size(RBprobes,1);
% RBprobes = [10*ones(nProbes,1) RBprobes nTrainBlocks*ones(nProbes,1) ... cat x y block
%             9999*ones(nProbes,1) ... filler
%             10*ones(nProbes,1) ... 10 = probe
%             randrows([11*ones(nProbes/2,1); 22*ones(nProbes/2,1)]) ...%dimension 2 infer 11=SF, 22=OR
%             randrows([ones(nProbes/2,1) ;2*ones(nProbes/2,1)])]; %A/notA or B/notB (only meaningful for training stimuli)]
%         
%scale stimuli(from Maddox et al. 2003)
% diameter_scale = .5;
% or_scale = 180/800;
% RB_low_scaled = [RB_low(:,1)   diameter_scale*RB_low(:,2)  RB_low(:,3)*or_scale];
% test_stim_scaled = [diameter_scale*test_stim(:,2) test_stim(:,3)*or_scale];
% RBprobes_scaled = [diameter_scale*RBprobes(:,2) RBprobes(:,3)*or_scale];
% RB_limits = [min(RB_low_scaled); max(RB_low_scaled)];
% 
% fishers_coeffs_scaled = fisherdiscrim2d([RB_low_scaled],3)
% 
% bestUDx_RBlow = fisherdiscrim2d([RB_low_scaled RB_low_scaled(:,1)],2)
% bestUDy_RBlow = fisherdiscrim2d([RB_low_scaled RB_low_scaled(:,1)],1)
% 
% %calculate mis-classifications
% UDaccuracy = 1 - ((sum(RB_low_scaled(RB_low_scaled(:,1)==1,2)<=bestUDx_RBlow(end)) + sum(RB_low_scaled(RB_low_scaled(:,1)==2,2)>bestUDx_RBlow(end))) / (nTrainBlocks*trialsPerBlock));
% 
% UDaccuracy = ...
%   (sum(RB_low_scaled(RB_low_scaled(:,1)==1,2)<=abs(bestUDx_RBlow(end)) & RB_low_scaled(RB_low_scaled(:,1)==1,3)<=abs(bestUDy_RBlow(end))) + ... sum of correct in lower left quad
%    sum(RB_low_scaled(RB_low_scaled(:,1)==2,2)<=abs(bestUDx_RBlow(end)) & RB_low_scaled(RB_low_scaled(:,1)==2,3)>abs(bestUDy_RBlow(end))) + ...    
%    sum(RB_low_scaled(RB_low_scaled(:,1)==1,2)>abs(bestUDx_RBlow(end)) & RB_low_scaled(RB_low_scaled(:,1)==1,3)>abs(bestUDy_RBlow(end))) + ...
%    sum(RB_low_scaled(RB_low_scaled(:,1)==2,2)>abs(bestUDx_RBlow(end)) & RB_low_scaled(RB_low_scaled(:,1)==2,3)<=abs(bestUDy_RBlow(end)))) ...
%    /(nTrainBlocks*trialsPerBlock)




h1 = figure(30); set(h1,'Position',[146 208 424 837],'name','Extrapolation'); hold on;
    subplot(2,2,1); 
        ax1 = [0 100 0 100];
        plot2dstim(RBcons,ax1,0); axis square; hold on; 
        %plot2dstim(IIcons,ax1,0); axis square; hold on; 
        %plot(TestII(:,2),TestII(:,3),'ro','MarkerSize',markerSize,'MarkerFaceColor','r');
        plot(TestRB(:,2),TestRB(:,3),'ro','MarkerSize',markerSize,'MarkerFaceColor','r');
        %plot(RBprobes(:,2),RBprobes(:,3),'cs','MarkerSize',.75*markerSize,'MarkerFaceColor','c');
        %set(gca,'Xtick',[]); set(gca,'Ytick',[])
        xlabel('Circle Diameter (arbitrary)','FontName',fontname,'FontSize',axisLabelSize); 
        ylabel('Orientation (arbitrary)','FontName',fontname,'FontSize',axisLabelSize);
        title('Unscaled Cons','FontName',fontname,'FontSize',titleSize); 
    subplot(2,2,2); 
        ax1 = [0 100 0 100];
        plot2dstim(RBincons,ax1,0); axis square; hold on; 
        %plot2dstim(IIincons,ax1,0); axis square; hold on;
        %plot(TestII(:,2),TestII(:,3),'ro','MarkerSize',markerSize,'MarkerFaceColor','r');
        plot(TestRB(:,2),TestRB(:,3),'ro','MarkerSize',markerSize,'MarkerFaceColor','r');
        %plot(RBprobes(:,2),RBprobes(:,3),'cs','MarkerSize',.75*markerSize,'MarkerFaceColor','c');
        %set(gca,'Xtick',[]); set(gca,'Ytick',[])
        xlabel('Circle Diameter (arbitrary)','FontName',fontname,'FontSize',axisLabelSize); 
        ylabel('Orientation (arbitrary)','FontName',fontname,'FontSize',axisLabelSize);
        title('Unscaled Incons','FontName',fontname,'FontSize',titleSize);
    subplot(2,2,3); 
        ax1 = [0 100 0 100];
        plot2dstim(IIcons2,ax1,0); axis square; hold on; 
        plot(TestII(:,2),TestII(:,3),'ro','MarkerSize',markerSize,'MarkerFaceColor','r');
        %plot(RBprobes(:,2),RBprobes(:,3),'cs','MarkerSize',.75*markerSize,'MarkerFaceColor','c');
        %set(gca,'Xtick',[]); set(gca,'Ytick',[])
        xlabel('Circle Diameter (arbitrary)','FontName',fontname,'FontSize',axisLabelSize); 
        ylabel('Orientation (arbitrary)','FontName',fontname,'FontSize',axisLabelSize);
        title('Unscaled Cons','FontName',fontname,'FontSize',titleSize); 
    subplot(2,2,4); 
        ax1 = [0 100 0 100]; 
        plot2dstim(IIincons2,ax1,0); axis square; hold on;
        plot(TestII(:,2),TestII(:,3),'ro','MarkerSize',markerSize,'MarkerFaceColor','r');
        %plot(RBprobes(:,2),RBprobes(:,3),'cs','MarkerSize',.75*markerSize,'MarkerFaceColor','c');
        %set(gca,'Xtick',[]); set(gca,'Ytick',[])
        xlabel('Circle Diameter (arbitrary)','FontName',fontname,'FontSize',axisLabelSize); 
        ylabel('Orientation (arbitrary)','FontName',fontname,'FontSize',axisLabelSize);
        title('Unscaled Incons','FontName',fontname,'FontSize',titleSize);        
        
        
%save stimulus files
  RBincons = [RBincons];
   save RBincons.stim RBincons -tabs -ascii
   save IIincons.stim IIincons -tabs -ascii
   save RBcons.stim RBcons -tabs -ascii
   save IIcons.stim IIcons -tabs -ascii
   save TestRB.stim TestRB -tabs -ascii
   save TestII.stim TestII -tabs -ascii
   %save Test_rotate.stim test_stim -tabs -ascii

