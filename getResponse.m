function [resp,correct,rt,time,data] = getResponse(window,sounds,feedback,color,where,stim,time,data,presExample,trial,grating,trainTest,trialWithinBlock)
%function [resp,correct,rt,time,data] = getResponse(window,sounds,feedback,color,where,stim,time,data,presExample,trial,grating)
%
%converts keyboard response to integer and determines if response was correct. Also computes RT
%
%input
%window - window pointer
%feedback - structure containing feedback images
%color - structure containing color information
%where - structure containing location information
%stim - structure containing stimulus information
%time - structure containing timing information
%data - data structure
%presExample - 0 for experimental trials, >0 for practice trials
%trial - current trial
%
%output
%resp - response
%correct - probabilistic fb
%opt_correct - optimal fb
%rt - response time
%time - updated timing structure
%
% 5/18/05   swe     written for work with set switching experiment
%6/11/05    swe     modified for PD wp v. ii exp
%8/10/05    swe     added flag to denote train/trans phase
%10/16/06   swe     modified for split-brain study, no longer unpack structres b/c inefficient
%12/3/06    swe     modified to work with OpenGL-Psychtoolbox version 3.0.8 on OSX 
%12/5/06    swe     removed manual keyboard buffer clear using KbCheck. This was pointless when also using KbCheck to 
%                       detect 'actual' response and introduced occasional delay (I'm guessing this was due to fast responses that happened to occur during this loop)
%1/23/07    swe     modified for imaging version of select-switch experiment (pilot)
%2/1/07     swe     modified for prob_select - reg focus
%2/22/07    swe     modified to use probabilistically (col 3) and optimally (col 4) correct responses in stimulus file rather than flipping a coin
%5/23/07    swe     modified for new button switch experiment (from splitCatExp_2rule_osxPT3_v2.m and button_switch_nr.m)
%4/22/08    swe     modified button_exp.m (i.e., button switch labels on screen vs. labels on keys) for relational learning experiment (relational_exp.m)
%7/20/2010  swe     modified for focalBG - wm experiment to be run @ Berkeley on PC
%5/14/12    swe     modified for stress cog pilot
%3/10/17    swe     modified to close gratings on the fly
%while(kbCheck);end %Wait for all keys to be released.
keyIsDown=0; 
while (keyIsDown==0); 
    [keyIsDown,response,time.response] = check_key(keyIsDown);
    if (GetSecs-time.stim_onset>time.responseDeadline); 
        response = '='; 
        time.response = GetSecs(); 
        keyIsDown=99;
        %break; 
    end
end %end while keyIsDown
Screen('Close', grating.g);
%GET CATEGORIZATION RESPONSE
correct = 99; resp = 99; rt = 0;
%compute category RT
rt = 1000*(time.response-time.stim_onset); %in ms

%TRANSFORM RESPONSE -- response assignments are counterbalanced
if data.catKeyMapping==1
    if sum(strcmp(response,cellstr(['n']')))==1
        resp = 1; %cat A
    elseif sum(strcmp(response,cellstr(['j']')))==1
        resp = 2; %cat B
    end
else
    if sum(strcmp(response,cellstr(['f']')))==1
        resp = 2; %cat B
    elseif sum(strcmp(response,cellstr(['v,']')))==1
        resp = 1; %cat A
    end
end

%other responses
if resp > 2 %assumes default resp different from 1-4
    EscapeSequence(response,'getResponse.m');
    %if response == '`'
        %quit experiment
        %escape('Q','user broke out of experiment');
    if strcmp(response,'space')
        resp = 10; 
    elseif sum(strcmp(response,cellstr(['njvf']')))==0 && keyIsDown==1 %invalid key press 
        resp = 10;
    end %end if response	
end

%COMPUTE FEEDBACK
if resp == 99 
    correct = 99;
elseif resp == 97
    correct = 97;
elseif resp == 10;
    correct = 2; 
elseif resp<=2
%     if  trainTest==2 || stim.stim(trialWithinBlock,1)==10 %no feedback trials
%         correct = -1;
%     elseif resp==stim.stim(trialWithinBlock,1)
    if resp==stim.stim(trialWithinBlock,1)
        correct = 1;
    else
        correct = 0;
    end
end %end if resp
  
 
if presExample==0 || presExample == 2 %fb for experimental and practice trials
    time = visual_feedback(window,correct,color,where,sounds,feedback,stim,time,trial,trialWithinBlock, trainTest);
elseif presExample==1 %no fb for example stimuli
    WaitSecs(time.fbDuration); time.fb_end = GetSecs;
end

%update ITI
if presExample==0 || presExample == 2 %fb for experimental and practice trials
    time = visual_feedback(window,correct,color,where,sounds,feedback,stim,time,trial,trialWithinBlock, trainTest);
elseif presExample==1 %no fb for example stimuli
    waitsecs(time.fbDuration); time.fb_end = GetSecs;
end

%time.ITI = max(time.trialDuration-(GetSecs-time.stim_onset)+time.ITI_base,time.ITI_base);
