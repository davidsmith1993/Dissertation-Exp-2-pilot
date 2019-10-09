 function genStim4extrapolation_interpolation_circles090814
%
%generates stimuli for extrapolation pilot
%using spatial frequency of large circular grating (at 90 deg) and angle of displaced line
%   need dimensions that are spatially separate so that one dimension can be omitted during inference task
%       this will hurt accuracy (based on relational experiment, but Ss should be able to learn)
%   need them to process both dimensions (i.e., shouldn't be able to ignore a dimension and learn)
%       grating should be at fixed position relative to line varying in angle
%   **try spatially centered grating (maybe some trial/trial jitter) with radial line that can vary from 0-270 degrees
%01/16/12   swe     
%2/2/12     swe     modified for circle diameter instead of sf
%                   distractors still generated, but going to try production task
%3/19/12    swe     added interpolation stimuli so that they extrap and interp can be combined in single experiment
%09/08/14   swe     modified to generate RB version
%07/21/16   swe     modified to generate RB version

close all
%figure stuff
fontname = 'ArialBold'; 
axisLabelSize = 16;
titleSize = 18;
textSize = 14;
linewidth = 2;
markerSize = 5;


nPerCat = 40;
nCat = 2;
nTrainBlocks = 4;


% %2/13/12 - shifted diameter to provide more room on low end
% Bmean_low = [300 150];
% Amean_low = [400 50];
% %variance-covariance matrices
% Q = [1/sqrt(2) -1/sqrt(2); 1/sqrt(2) 1/sqrt(2)];
% D = [4000 0; 0 500];
% Kpos = Q*D*Q';
% shift_cat = 400;
% 
% %3/19/12 - increased variance along main diagonal (to decrease accuracy of UD strategies) and shifted means (to preserve spacing between categories)
% Bmean_low = [215 35];
% Amean_low = [285 -35];
% %variance-covariance matrices
% Q = [1/sqrt(2) -1/sqrt(2); 1/sqrt(2) 1/sqrt(2)];
% D = [6000 0; 0 250];
% Kpos = Q*D*Q';
% shift_cat = 500;
%9/8/14 - shifted means to allow for RB structure
%7/19/2016- edited for E3, 4 catagories

Amean_low = [650 250];
Bmean_low = [350 250];
Cmean_low = [350 -50];
Dmean_low = [650 -50];
%variance-covariance matrices
NegQ = [1/sqrt(2) -1/sqrt(2); 1/sqrt(2) 1/sqrt(2)];
PosQ = [-1/sqrt(2) 1/sqrt(2); 1/sqrt(2) 1/sqrt(2)];
D = [7500 0; 0 250];
Kpos = PosQ*D*PosQ';
Kneg = NegQ*D*NegQ';
shift_cat = 500;
%RB stim
Amean_lowRB = [Amean_low(1) Amean_low(2)];
Bmean_lowRB = [Bmean_low(1) Bmean_low(2)];
Cmean_lowRB = [Cmean_low(1) Cmean_low(2)];
Dmean_lowRB = [Dmean_low(1) Dmean_low(2)];


prc = linprobcorr(Kpos,Amean_low',Bmean_low')





RB_low = zeros(nPerCat*nCat*nTrainBlocks,8);

trialsPerBlock = nPerCat*nCat;
for b = 1:nTrainBlocks
    %cols are: cat, diameter, orientation, block,  ...
    %stimlus_location (used to matter for 2AFC version of experiment, but just a filler for this version included to minimize changes to experiment program)
    %range(1=low,2=high) ...
    %dimension_to_infer(11=sf, 22=or)
    %A/notA or B/notB (only meaningful for training stimuli)]


    catA = [genMsamples(Amean_lowRB',Kneg,nPerCat/2,1,1)  b*ones(nPerCat/2,1) ...
               9999*ones(nPerCat/2,1) ...
               ones(nPerCat/2,1) ...
               randrows([11*ones(nPerCat/4,1); 22*ones(nPerCat/4,1)]) ...
               randrows([ones(nPerCat/4,1) ;2*ones(nPerCat/4,1)])];
    catB = [genMsamples(Bmean_lowRB',Kpos,nPerCat/2,2,1)  b*ones(nPerCat/2,1) ...
               9999*ones(nPerCat/2,1) ...
               ones(nPerCat/2,1) ...
               randrows([11*ones(nPerCat/4,1); 22*ones(nPerCat/4,1)]) ...
               randrows([ones(nPerCat/4,1) ;2*ones(nPerCat/4,1)])];
    catC = [genMsamples(Cmean_lowRB',Kneg,nPerCat/2,1,1)  b*ones(nPerCat/2,1) ...
               9999*ones(nPerCat/2,1) ...
               ones(nPerCat/2,1) ...
               randrows([11*ones(nPerCat/4,1); 22*ones(nPerCat/4,1)]) ...
               randrows([ones(nPerCat/4,1) ;2*ones(nPerCat/4,1)])];
    catD = [genMsamples(Dmean_lowRB',Kpos,nPerCat/2,2,1)  b*ones(nPerCat/2,1) ...
               9999*ones(nPerCat/2,1) ...
               ones(nPerCat/2,1) ...
               randrows([11*ones(nPerCat/4,1); 22*ones(nPerCat/4,1)]) ...
               randrows([ones(nPerCat/4,1) ;2*ones(nPerCat/4,1)])];
%     RB_low(b*trialsPerBlock-trialsPerBlock+1:b*trialsPerBlock,:) = [catA;catB;catC;catD]; 
RB_low(b*trialsPerBlock-trialsPerBlock+1:b*trialsPerBlock,:) = [rotate2dDist(catA,1.508);rotate2dDist(catB,1.508);rotate2dDist(catC,1.508);rotate2dDist(catD,1.508)]; 
end
fishers_coeffs_unscaled = fisherdiscrim2d([RB_low],3);


%generate test stimuli
%extrapolation along correlated dimension
nTest = 112; %evenly divisible by nClusters and nClusters*2
nClusters = 8; %number of clusters of test stimuli-4 within training, 4 extrapolation
%x=75
x = 100;
Atest_low =  [      ones(1,nTest/nClusters); ... %A test items = 1, B test items = 2
                    linspace(Amean_low(1)-x,Amean_low(1)+x,nTest/nClusters);  ... x coords
                    linspace(Amean_low(2)-x,Amean_low(2)+x,nTest/nClusters); ... %y coords
                    (nTrainBlocks+1)*ones(1,nTest/nClusters); ... block number
                    [randrows([-1*ones(nTest/(2*nClusters),1);ones(nTest/(2*nClusters),1)])]'; ... %target location for inference trials -1=left (filler)
                    ones(1,nTest/nClusters); ... %in low range = 1; in high range = 2; in mid range = 3;
                    [randrows([11*ones(nTest/(2*nClusters),1); 22*ones(nTest/(2*nClusters),1)])]'; ... %dimension 2 infer 11=SF, 22=OR
                    9999*ones(1,nTest/nClusters)]; %filler
Btest_low =  [      2*ones(1,nTest/nClusters); ... %A test items = 1, B test items = 2
                    linspace(Bmean_low(1)+x,Bmean_low(1)-x,nTest/nClusters);  ... x coords
                    linspace(Bmean_low(2)-x,Bmean_low(2)+x,nTest/nClusters); ... %y coords
                    (nTrainBlocks+1)*ones(1,nTest/nClusters); ... block number
                    [randrows([-1*ones(nTest/(2*nClusters),1);ones(nTest/(2*nClusters),1)])]'; ... %target location for inference trials -1=left (filler)
                    ones(1,nTest/nClusters); ... %in low range = 1; in high range = 2; in mid range = 3;
                    [randrows([11*ones(nTest/(2*nClusters),1); 22*ones(nTest/(2*nClusters),1)])]'; ... %dimension 2 infer 11=SF, 22=OR
                    9999*ones(1,nTest/nClusters)]; %filler
Ctest_low =  [      3*ones(1,nTest/nClusters); ... %A test items = 1, B test items = 2
                    linspace(Cmean_low(1)-x,Cmean_low(1)+x,nTest/nClusters);  ... x coords
                    linspace(Cmean_low(2)-x,Cmean_low(2)+x,nTest/nClusters); ... %y coords
                    (nTrainBlocks+1)*ones(1,nTest/nClusters); ... block number
                    [randrows([-1*ones(nTest/(2*nClusters),1);ones(nTest/(2*nClusters),1)])]'; ... %target location for inference trials -1=left (filler)
                    ones(1,nTest/nClusters); ... %in low range = 1; in high range = 2; in mid range = 3;
                    [randrows([11*ones(nTest/(2*nClusters),1); 22*ones(nTest/(2*nClusters),1)])]'; ... %dimension 2 infer 11=SF, 22=OR
                    9999*ones(1,nTest/nClusters)]; %filler
Dtest_low =  [      4*ones(1,nTest/nClusters); ... %A test items = 1, B test items = 2
                    linspace(Dmean_low(1)+x,Dmean_low(1)-x,nTest/nClusters);  ... x coords
                    linspace(Dmean_low(2)-x,Dmean_low(2)+x,nTest/nClusters); ... %y coords
                    (nTrainBlocks+1)*ones(1,nTest/nClusters); ... block number
                    [randrows([-1*ones(nTest/(2*nClusters),1);ones(nTest/(2*nClusters),1)])]'; ... %target location for inference trials -1=left (filler)
                    ones(1,nTest/nClusters); ... %in low range = 1; in high range = 2; in mid range = 3;
                    [randrows([11*ones(nTest/(2*nClusters),1); 22*ones(nTest/(2*nClusters),1)])]'; ... %dimension 2 infer 11=SF, 22=OR
                    9999*ones(1,nTest/nClusters)]; %filler
                
test_stim = [rotate2dDist(Atest_low',1.508); rotate2dDist(Btest_low',1.508); rotate2dDist(Ctest_low',1.508); rotate2dDist(Dtest_low',1.508);];


%probe stimuli 7/27/2016
Aprobemean = [650 100];
Bprobemean = [500 250];
Cprobemean = [350 100];
Dprobemean = [500 -50];

t=60;
nProbe = 32; %evenly divisible by nClustersprobe and nClustersprobe*2
nClustersProbe = 4; %number of clusters of test stimuli-4 within training, 4 extrapolation

Aprobe =  [      ones(1,nProbe/nClustersProbe); ... %A test items = 1, B test items = 2
                    linspace(Aprobemean(1),Aprobemean(1),nProbe/nClustersProbe);  ... x coords
                    linspace(Aprobemean(2)-t,Aprobemean(2)+t,nProbe/nClustersProbe); ... %y coords
                    (nTrainBlocks+1)*ones(1,nProbe/nClustersProbe); ... block number
                    [randrows([-1*ones(nProbe/(2*nClustersProbe),1);ones(nProbe/(2*nClustersProbe),1)])]'; ... %target location for inference trials -1=left (filler)
                    ones(1,nProbe/nClustersProbe); ... %in low range = 1; in high range = 2; in mid range = 3;
                    [randrows([11*ones(nProbe/(2*nClustersProbe),1); 22*ones(nProbe/(2*nClustersProbe),1)])]'; ... %dimension 2 infer 11=SF, 22=OR
                    9999*ones(1,nProbe/nClustersProbe)]; %filler
Bprobe =  [      2*ones(1,nProbe/nClustersProbe); ... %A test items = 1, B test items = 2
                    linspace(Bprobemean(1)+t,Bprobemean(1)-t,nProbe/nClustersProbe);  ... x coords
                    linspace(Bprobemean(2),Bprobemean(2),nProbe/nClustersProbe); ... %y coords
                    (nTrainBlocks+1)*ones(1,nProbe/nClustersProbe); ... block number
                    [randrows([-1*ones(nProbe/(2*nClustersProbe),1);ones(nProbe/(2*nClustersProbe),1)])]'; ... %target location for inference trials -1=left (filler)
                    ones(1,nProbe/nClustersProbe); ... %in low range = 1; in high range = 2; in mid range = 3;
                    [randrows([11*ones(nProbe/(2*nClustersProbe),1); 22*ones(nProbe/(2*nClustersProbe),1)])]'; ... %dimension 2 infer 11=SF, 22=OR
                    9999*ones(1,nProbe/nClustersProbe)]; %filler
Cprobe =  [      3*ones(1,nProbe/nClustersProbe); ... %A test items = 1, B test items = 2
                    linspace(Cprobemean(1),Cprobemean(1),nProbe/nClustersProbe);  ... x coords
                    linspace(Cprobemean(2)-t,Cprobemean(2)+t,nProbe/nClustersProbe); ... %y coords
                    (nTrainBlocks+1)*ones(1,nProbe/nClustersProbe); ... block number
                    [randrows([-1*ones(nProbe/(2*nClustersProbe),1);ones(nProbe/(2*nClustersProbe),1)])]'; ... %target location for inference trials -1=left (filler)
                    ones(1,nProbe/nClustersProbe); ... %in low range = 1; in high range = 2; in mid range = 3;
                    [randrows([11*ones(nProbe/(2*nClustersProbe),1); 22*ones(nProbe/(2*nClustersProbe),1)])]'; ... %dimension 2 infer 11=SF, 22=OR
                    9999*ones(1,nProbe/nClustersProbe)]; %filler
Dprobe =  [      4*ones(1,nProbe/nClustersProbe); ... %A test items = 1, B test items = 2
                    linspace(Dprobemean(1)+t,Dprobemean(1)-t,nProbe/nClustersProbe);  ... x coords
                    linspace(Dprobemean(2),Dprobemean(2),nProbe/nClustersProbe); ... %y coords
                    (nTrainBlocks+1)*ones(1,nProbe/nClustersProbe); ... block number
                    [randrows([-1*ones(nProbe/(2*nClustersProbe),1);ones(nProbe/(2*nClustersProbe),1)])]'; ... %target location for inference trials -1=left (filler)
                    ones(1,nProbe/nClustersProbe); ... %in low range = 1; in high range = 2; in mid range = 3;
                    [randrows([11*ones(nProbe/(2*nClustersProbe),1); 22*ones(nProbe/(2*nClustersProbe),1)])]'; ... %dimension 2 infer 11=SF, 22=OR
                    9999*ones(1,nProbe/nClustersProbe)]; %filler

%probe stimuli for RB task, block 4.
RBprobes = [Aprobe(2:3,1:2:end) Bprobe(2:3,1:2:end) Cprobe(2:3,1:2:end) Dprobe(2:3,1:2:end)]';
nProbes = size(RBprobes,1);
RBprobes = [10*ones(nProbes,1) RBprobes nTrainBlocks*ones(nProbes,1) ... cat x y block
            9999*ones(nProbes,1) ... filler
            10*ones(nProbes,1) ... 10 = probe
            randrows([11*ones(nProbes/2,1); 22*ones(nProbes/2,1)]) ...%dimension 2 infer 11=SF, 22=OR
            randrows([ones(nProbes/2,1) ;2*ones(nProbes/2,1)])]; %A/notA or B/notB (only meaningful for training stimuli)]
        
%scale stimuli(from Maddox et al. 2003)
diameter_scale = .5;
or_scale = 180/800;
RB_low_scaled = [RB_low(:,1)   diameter_scale*RB_low(:,2)  RB_low(:,3)*or_scale];
test_stim_scaled = [diameter_scale*test_stim(:,2) test_stim(:,3)*or_scale];
RBprobes_scaled = [diameter_scale*RBprobes(:,2) RBprobes(:,3)*or_scale];
RB_limits = [min(RB_low_scaled); max(RB_low_scaled)];

fishers_coeffs_scaled = fisherdiscrim2d([RB_low_scaled],3)

bestUDx_RBlow = fisherdiscrim2d([RB_low_scaled RB_low_scaled(:,1)],2)
bestUDy_RBlow = fisherdiscrim2d([RB_low_scaled RB_low_scaled(:,1)],1)

%calculate mis-classifications
UDaccuracy = 1 - ((sum(RB_low_scaled(RB_low_scaled(:,1)==1,2)<=bestUDx_RBlow(end)) + sum(RB_low_scaled(RB_low_scaled(:,1)==2,2)>bestUDx_RBlow(end))) / (nTrainBlocks*trialsPerBlock));

UDaccuracy = ...
  (sum(RB_low_scaled(RB_low_scaled(:,1)==1,2)<=abs(bestUDx_RBlow(end)) & RB_low_scaled(RB_low_scaled(:,1)==1,3)<=abs(bestUDy_RBlow(end))) + ... sum of correct in lower left quad
   sum(RB_low_scaled(RB_low_scaled(:,1)==2,2)<=abs(bestUDx_RBlow(end)) & RB_low_scaled(RB_low_scaled(:,1)==2,3)>abs(bestUDy_RBlow(end))) + ...    
   sum(RB_low_scaled(RB_low_scaled(:,1)==1,2)>abs(bestUDx_RBlow(end)) & RB_low_scaled(RB_low_scaled(:,1)==1,3)>abs(bestUDy_RBlow(end))) + ...
   sum(RB_low_scaled(RB_low_scaled(:,1)==2,2)>abs(bestUDx_RBlow(end)) & RB_low_scaled(RB_low_scaled(:,1)==2,3)<=abs(bestUDy_RBlow(end)))) ...
   /(nTrainBlocks*trialsPerBlock)




h = figure(30); set(h,'Position',[146 208 424 837],'name','Extrapolation'); hold on;
    subplot(2,1,1); 
        ax1 = [0 1000 -200 800];
        plot2dstim(RB_low,ax1,0); axis square; hold on; 
        plot(test_stim(:,2),test_stim(:,3),'ro','MarkerSize',markerSize,'MarkerFaceColor','r');
        plot(RBprobes(:,2),RBprobes(:,3),'cs','MarkerSize',.75*markerSize,'MarkerFaceColor','c');
        %set(gca,'Xtick',[]); set(gca,'Ytick',[])
        xlabel('Circle Diameter (arbitrary)','FontName',fontname,'FontSize',axisLabelSize); 
        ylabel('Orientation (arbitrary)','FontName',fontname,'FontSize',axisLabelSize);
        title('Unscaled','FontName',fontname,'FontSize',titleSize); 
   subplot(2,1,2); 
        ax2 = [0    diameter_scale*ax1(2)    ax1(3)*or_scale   ax1(4)*or_scale];
        plot2dstim(RB_low_scaled,ax2,0); axis square; hold on; 
        plot(test_stim_scaled(:,1),test_stim_scaled(:,2),'ro','MarkerSize',markerSize,'MarkerFaceColor','r');
        plot(RBprobes_scaled(:,1),RBprobes_scaled(:,2),'cs','MarkerSize',.75*markerSize,'MarkerFaceColor','c');
        hold on; plot2dlinbnd(bestUDx_RBlow,'k-',ax2)
        hold on; plot2dlinbnd(bestUDy_RBlow,'k-',ax2)
        %hold on; plot(abs([bestUDx_RBlow(3) bestUDx_RBlow(3)]),ax2(3:4),'k-') %plots scaled x bound
        %hold on; plot(ax2(1:2),abs([bestUDy_RBlow(3) bestUDy_RBlow(3)]),'k-') %plots scaled y bound
        xlabel('Circle Diameter (pixels)','FontName',fontname,'FontSize',axisLabelSize); 
        ylabel('Orientation (degrees from horizontal)','FontName',fontname,'FontSize',axisLabelSize);
        title('Scaled','FontName',fontname,'FontSize',titleSize); 

% %save stimulus files
  RB_low = [RB_low; RBprobes];
   save RB_rotate.stim RB_low -tabs -ascii
   save Test_rotate.stim test_stim -tabs -ascii

