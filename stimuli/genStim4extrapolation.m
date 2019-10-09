function genStim4extrapolation
%
%generates stimuli for extrapolation pilot
%using spatial frequency of large circular grating (at 90 deg) and angle of displaced line
%   need dimensions that are spatially separate so that one dimension can be omitted during inference task
%       this will hurt accuracy (based on relational experiment, but Ss should be able to learn)
%   need them to process both dimensions (i.e., shouldn't be able to ignore a dimension and learn)
%       grating should be at fixed position relative to line varying in angle
%   **try spatially centered grating (maybe some trial/trial jitter) with radial line that can vary from 0-270 degrees
%01/16/12   swe     
%
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
nTrainBlocks = 6;
%variance-covariance matrices
Q = [1/sqrt(2) -1/sqrt(2); 1/sqrt(2) 1/sqrt(2)];
D = [2500 0; 0 300];
Kpos = Q*D*Q'

%generate training stimuli
Bmean_low = [100 200];
Amean_low = [150 150];
prc = linprobcorr(Kpos,Amean_low',Bmean_low')
shift_cat = 150;
Bmean_high = Bmean_low + shift_cat;
Amean_high = Amean_low + shift_cat;
II_low = zeros(nPerCat*nCat*nTrainBlocks,4); II_high = II_low;
trialsPerBlock = nPerCat*nCat;
for b = 1:nTrainBlocks
    II_low(b*trialsPerBlock-trialsPerBlock+1:b*trialsPerBlock,:) = [genMsamples([Amean_low' Bmean_low'],Kpos,[nPerCat nPerCat],[1 2],nCat)  b*ones(trialsPerBlock,1)];
    II_high(b*trialsPerBlock-trialsPerBlock+1:b*trialsPerBlock,:) = [genMsamples([Amean_high' Bmean_high'],Kpos,[nPerCat nPerCat],[1 2],nCat)  b*ones(trialsPerBlock,1)];
end
MINIMUMVALUES = min(II_low)
sum(II_low(:,2)<1)
%generate test stimuli
%extrapolation along correlated dimension
nTest = trialsPerBlock;
x = 30;
Atest_low =  [      ones(1,nTest/4); ... %A test items = 1, B test items = 2
                    linspace(Amean_low(1)-x,Amean_low(1)+x,nTest/4);  ... x coords
                    linspace(Amean_low(2)-x,Amean_low(2)+x,nTest/4); ... %y coords
                    (nTrainBlocks+1)*ones(1,nTest/4); ... block number
                    [randrows([-1*ones(nTest/8,1);ones(nTest/8,1)])]'; ... %target location for inference trials -1=left
                    ones(1,nTest/4); ... %in low range = 1; in high range = 2;
                    [randrows([11*ones(nTest/8,1); 22*ones(nTest/8,1)])]']; %dimension 2 infer 11=SF, 22=OR
Atest_high = [      ones(1,nTest/4); ... %A test items = 1, B test items = 2
                    linspace(Amean_high(1)-x,Amean_high(1)+x,nTest/4); ... x coords
                    linspace(Amean_high(2)-x,Amean_high(2)+x,nTest/4); ... %y coords
                    (nTrainBlocks+1)*ones(1,nTest/4); ... block number
                    [randrows([-1*ones(nTest/8,1);ones(nTest/8,1)])]'; ... %target location for inference trials -1=left
                    2*ones(1,nTest/4); ... %in low range = 1; in high range = 2;
                    [randrows([11*ones(nTest/8,1); 22*ones(nTest/8,1)])]']; %dimension 2 infer 11=SF, 22=OR
         
Btest_low =  [      2*ones(1,nTest/4); ... %A test items = 1, B test items = 2
                    linspace(Bmean_low(1)-x,Bmean_low(1)+x,nTest/4);  ... x coords
                    linspace(Bmean_low(2)-x,Bmean_low(2)+x,nTest/4); ... %y coords
                    (nTrainBlocks+1)*ones(1,nTest/4); ... block number
                    [randrows([-1*ones(nTest/8,1);ones(nTest/8,1)])]'; ... %target location for inference trials -1=left
                    ones(1,nTest/4); ... %in low range = 1; in high range = 2;
                    [randrows([11*ones(nTest/8,1); 22*ones(nTest/8,1)])]']; %dimension 2 infer 11=SF, 22=OR
Btest_high = [      2*ones(1,nTest/4); ... %A test items = 1, B test items = 2
                    linspace(Bmean_high(1)-x,Bmean_high(1)+x,nTest/4); ... x coords
                    linspace(Bmean_high(2)-x,Bmean_high(2)+x,nTest/4); ... %y coords
                    (nTrainBlocks+1)*ones(1,nTest/4); ... block number
                    [randrows([-1*ones(nTest/8,1);ones(nTest/8,1)])]'; ... %target location for inference trials -1=left
                    2*ones(1,nTest/4); ... %in low range = 1; in high range = 2;
                    [randrows([11*ones(nTest/8,1); 22*ones(nTest/8,1)])]']; %dimension 2 infer 11=SF, 22=OR
         
test_stim = [Atest_low'; Atest_high'; Btest_low'; Btest_high'];
nTest_check = size(test_stim,1)
%scale stimuli(from Maddox et al. 2003)
sf_scale = [.25 50];
or_scale = 180/800;
II_low_scaled = [II_low(:,1) sf_scale(1)+II_low(:,2)/sf_scale(2) II_low(:,3)*or_scale];
II_high_scaled = [II_high(:,1) sf_scale(1)+II_high(:,2)/sf_scale(2) II_high(:,3)*or_scale];
test_stim_scaled = [sf_scale(1)+test_stim(:,2)/sf_scale(2) test_stim(:,3)*or_scale];

h = figure(30); set(h,'Position',[146 208 898 837],'name','Extrapolation'); hold on;
    subplot(2,2,1); 
        ax1 = [0 450 50 500];
        plot2dstim(II_low,ax1,0); axis square; hold on; 
        plot(test_stim(:,2),test_stim(:,3),'ro','MarkerSize',markerSize,'MarkerFaceColor','r');
        title('Lower','FontName',fontname,'FontSize',titleSize); 
        %set(gca,'Xtick',[]); set(gca,'Ytick',[])
    subplot(2,2,3); 
        plot2dstim(II_high,ax1,0); axis square; hold on; 
        plot(test_stim(:,2),test_stim(:,3),'ro','MarkerSize',markerSize,'MarkerFaceColor','r');     
        %set(gca,'Xtick',[]); set(gca,'Ytick',[])
        xlabel('Spatial Frequency (arbitrary)','FontName',fontname,'FontSize',axisLabelSize); 
        ylabel('Orientation (arbitrary)','FontName',fontname,'FontSize',axisLabelSize);
        title('Upper','FontName',fontname,'FontSize',titleSize); 
    subplot(2,2,2); 
        ax2 = [0    sf_scale(1)+ax1(2)/sf_scale(2)    0    ax1(4)*or_scale];
        plot2dstim(II_low_scaled,ax2,0); axis square; hold on; 
        plot(test_stim_scaled(:,1),test_stim_scaled(:,2),'ro','MarkerSize',markerSize,'MarkerFaceColor','r');
    subplot(2,2,4); 
        plot2dstim(II_high_scaled,ax2,0); axis square; hold on; 
        plot(test_stim_scaled(:,1),test_stim_scaled(:,2),'ro','MarkerSize',markerSize,'MarkerFaceColor','r');
        xlabel('Spatial Frequency (cycles/degree)','FontName',fontname,'FontSize',axisLabelSize); 
        ylabel('Orientation (degrees from horizontal)','FontName',fontname,'FontSize',axisLabelSize);
%check overlap of low and high conditions
h=figure(31); set(h,'Position',[146 208 898 837],'name','Extrapolation-check'); hold on;
        plot2dstim(II_low,ax1,0); axis square; hold on; 
        plot2dstim(II_high,ax1,0); axis square; hold on; 
        plot(test_stim(:,2),test_stim(:,3),'ro','MarkerSize',markerSize,'MarkerFaceColor','r');   
        plot([ax1(1) ax1(2)],[90*1/or_scale  90*1/or_scale],'k-')
        plot([ax1(1) ax1(2)],[45*1/or_scale  45*1/or_scale],'k-')
%overlay rotated test distribution for distractors     
fishers_coeffs = fisherdiscrim2d([II_low; II_high],3);
plot2dlinbnd(fishers_coeffs,'-r',ax1)
temp = rotate2dDist([II_low(:,1:3); II_high(:,1:3)],90*pi/180);
fishers_coeffs_rot = fisherdiscrim2d(temp,3);
%plot2dlinbnd(fishers_coeffs_rot,'-b',ax)    
%subject shown SF and asked to select correct OR from 2 possible orientations (1 correct, 1 orientation if there were negative correlation)
distractorsOR = [test_stim(:,2) ((fishers_coeffs_rot(1)*test_stim(:,2) + fishers_coeffs_rot(3))/-fishers_coeffs_rot(2))];
plot(distractorsOR(:,1),distractorsOR(:,2),'b*','MarkerSize',markerSize,'MarkerFaceColor','b');
%subject shown OR and asked to select correct SF from 2 possible SFs (1 correct, 1 SF if there were negative correlation)
distractorsSF = [((fishers_coeffs_rot(2)*test_stim(:,3) + fishers_coeffs_rot(3))/-fishers_coeffs_rot(1)) test_stim(:,3) ];
plot(distractorsSF(:,1),distractorsSF(:,2),'go','MarkerSize',markerSize);

%8 cols of test_stim are: Based_on_what_category(1=A,2=B),   sf,   or,   block,   range(1=low,2=high),   dimension_to_infer(11=sf, 22=or), distractorSF, distractorOR]
%   if dimension_to_infer == 11, the value in the distractorSF column is used, else the value in distractorOR is used
test_stim = [test_stim distractorsSF(:,1) distractorsOR(:,2)];

%add cols to training stimuli to parallel transfer stimuli
%cols are: cat, sf, or, block, range(1=low,2=high) present_target on left (-1) or right (1)  dimension_to_infer(11=sf, 22=or) distractorSF distractorOR]
II_low = [II_low ones(size(II_low,1),1)  9999*ones(size(II_low,1),1) 9999*ones(size(II_low,1),1) 9999*ones(size(II_low,1),1) 9999*ones(size(II_low,1),1)];
II_high = [II_high 2*ones(size(II_high,1),1) 9999*ones(size(II_low,1),1) 9999*ones(size(II_high,1),1) 9999*ones(size(II_high,1),1) 9999*ones(size(II_high,1),1)];
%save stimulus files
save II_train_low.stim II_low -tabs -ascii
save II_train_high.stim II_high -tabs -ascii
save extrapolateTest.stim test_stim -tabs -ascii

checkTestStim = 1;
if checkTestStim
    %double check test stimuli
    h=figure(32); set(h,'Position',[146 208 898 837],'name','Extrapolation-check One at a time'); hold on;
        plot2dstim(II_low,ax1,0); axis square; hold on; 
        plot2dstim(II_high,ax1,0); axis square; hold on; 
        xlabel('Spatial Frequency (arbitrary)','FontName',fontname,'FontSize',axisLabelSize); 
        ylabel('Orientation (arbitrary)','FontName',fontname,'FontSize',axisLabelSize);
    for i=1:nTest
        %infer SF
        %plot(test_stim(i,2), test_stim(i,3),'ro','MarkerSize',markerSize,'MarkerFaceColor','r') %correct
        %plot(test_stim(i,8), test_stim(i,3),'bo','MarkerSize',markerSize,'MarkerFaceColor','b') %incorrect
        %infer OR
        plot(test_stim(i,2), test_stim(i,3),'ro','MarkerSize',markerSize,'MarkerFaceColor','r') %correct
        plot(test_stim(i,2), test_stim(i,9),'go','MarkerSize',markerSize,'MarkerFaceColor','g') %incorrect
        pause
    end
end


