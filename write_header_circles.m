function write_header_circles(data,stim,time,where)
disp('write header circles'); save temp;
%function write_header(data,stim,time,where)
%
%writes header file for pd_wp_vs_ii experiment 
%
%10/16/06   swe     modified for split-brain study, no longer unpack structres b/c inefficient
%5/23/07    swe     modified for new button switch experiment (from splitCatExp_2rule_osxPT3_v2.m and button_switch_nr.m)
%4/22/08    swe     modified button_exp.m (i.e., button switch labels on screen vs. labels on keys) for relational learning experiment (relational_exp.m)
%7/20/2010  swe     modified relationa_bonus.m for focalBG - wm experiment to be run @ Berkeley on PC
%1/17/2012  swe     modified for cat-inference experiment
%2/3/12     swe     modified for new stimuli (circle-line) and production task (still using extrapolateTest.stim, but the distractor values are irrelevant

%disp('here'); save temp

% Create header file - i love redundancy
cd data
date_str = date;
fid = fopen([data.outname '.hdr'],'w');
fprintf(fid,'filename = %s \n',data.outname);
fprintf(fid,'date = %s \n\n',date_str);

fprintf(fid,'\nCONDITION INFORMATION:\n');
fprintf(fid,'Subject Number = %i \n',data.sub);
fprintf(fid,'Gender = %s\n', data.genderText);
fprintf(fid,'Task = %s\n', data.categoryStructureText);
fprintf(fid,'Training Condition = %s\n', data.trainingConditionText);
fprintf(fid,'Transfer Condition = %s\n', data.transferConditionText);

fprintf(fid,'\nTIMING VARIABLES (in seconds):\n');
fprintf(fid,'fixation duration = %6.3f \n',time.fixTime);
fprintf(fid,'stimulus duration = %6.3f \n',time.stimDuration);
fprintf(fid,'Subject must respond within = %6.3f \n',time.responseDeadline);
fprintf(fid,'fb duration = %6.3f \n',time.fbDuration);
fprintf(fid,'ITI = %6.3f \n',time.ITI);
fprintf(fid,'Experiment began = % g % g % g % g % g % g \n',time.experiment_begin);
fprintf(fid,'Experiment duration (minutes) = %6.3f \n',(time.now-time.starttime)/60);

fprintf(fid,'\nPOSITION VARIABLES:\n');
fprintf(fid,'Assumed distance between subject and monitor = %f \n',where.viewing_distance);

fprintf(fid,'\nSTIMULUS VARIABLES:\n');
fprintf(fid,'Stimulus: circles varying in diameter + line varying in angle at fixed length. \n');
fprintf(fid,'Diameter scaling factor = %f \n',1);
fprintf(fid,'(i.e., 1 = unscaled) \n');
fprintf(fid,'Orientation scaling factor = %f \n',data.o_gain);
fprintf(fid,'orientation scaling: tiltInRadians = tiltInDegrees * pi / data.o_gain. \n');
fprintf(fid,'Length unscaled:\n');


fprintf(fid,'\nRESPONSE VARIABLES:\n');
fprintf(fid,'Category-keyboard mapping = %s \n',data.catKeyMapping_text);

fprintf(fid,'\nOTHER VARIABLES:\n');
fprintf(fid,'number of training blocks = %i \n',data.numTrainingBlocks);
fprintf(fid,'number of transfer blocks = %i \n',data.numTransferBlocks);
fprintf(fid,'Trials per block = %i \n',data.trialsPerBlock);
fprintf(fid,'Total number of trials = %i \n',data.numTrials);
fprintf(fid,'Monitor size = %s \n',data.monitor_text);
fprintf(fid,'Monitor resolution = %s \n',data.monitor_res_text);


fprintf(fid,'\nPercent correct by block :\n');
fprintf(fid,'%6.2f\n',data.block_pc(:,1)');
	
fprintf(fid,'\nRT (in ms) by block (includes correct trials with RT < data.rt_cutoff):\n');
fprintf(fid,'%6.2f\n',data.block_rt(:,1)');
fprintf(fid,'RT cutoff used to compute avg RT = %f \n',data.rt_cutoff);

fprintf(fid,'\nRMSE [diameter] by block :\n');
fprintf(fid,'%6.2f\n',data.block_rmse(:,1)');

fprintf(fid,'\nRMSE [orientation] by block :\n');
fprintf(fid,'%6.2f\n',data.block_rmse(:,2)');

fprintf(fid,'\n\nDATA FILE FORMAT:\n');
fprintf(fid,'     columns are columns are [cat, un-scaled_sf, un-scaled_orient, resp, RT, fb, diameter response (inference only), orient response (inference only), currentCondition (1=A/B,2=inf,3=AnA), trainTest (1=train), in low/high range (1=low,2=high),\n'); 
fprintf(fid,'                              inferDimension (11=infer SF,22=infer OR), AnotA question (AnA only), trial, trialWithinBlock, block]\n');
fprintf(fid,'     rows of data file are trials\n');
fclose(fid);
cd ..



