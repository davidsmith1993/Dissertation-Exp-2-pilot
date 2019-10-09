function [time,grating] = present_stim(window,data,color,where,stim,trialWithinBlock,trial,time,presExample,clearmode)
%function [time,grating] = present_stim(window,data,color,where,stim,trialWithinBlock,trial,time,presExample,clearmode)
%
%generates and presents stimulus 
%
%10/15/06   swe     
%12/3/06    swe     modified to work with OpenGL-Psychtoolbox version 3.0.8 on OSX -- using Screen(DrawTexture)
%12/11/06   swe     modified for catch trials (fixation dims)
%1/23/07    swe     modified for imaging version of select-switch experiment (pilot)
%2/1/07     swe     modified for prob selection - regulatory focus experiment
%3/28/08    swe     modified for relational exp 
%4/29/08    swe     added clearmode as input
%11/16/08   swe     present two images in single display
%7/20/2010  swe     modified for focalBG - wm experiment to be run @ Berkeley on PC
%3/10/17    swe     modified to present gratings on the fly for online stressor experiment

try
    %present stimuli w/ ITI
    % Draw fixation point
    %SCREEN('DrawTexture',window,stim.fix,[],[where.xc-stim.fix_size/2, ...
    %                                         where.yc-stim.fix_size/2, ...
    %                                         where.xc+stim.fix_size/2, ...
    %                                         where.yc+stim.fix_size/2],[],0);
                                     %Screen('DrawText',window,stim.fixation,where.xc,where.yc,color.textColor);
    % Show fixation approximately time.initial_delay seconds after time 'starttime', record time of real stimulus onset in vbl: 
    %The clearmode flag allows you to specify what should happen AFTER a Flip:
    %1. A value of clearmode=0 (the default setting) will clear the drawing surface to the background color, so you can draw a completely new stimulus.
    %2. A value of clearmode=1 will restore the drawing surface to exactly the same state as before the flip, so you can incrementally update your stimulus.
    %3. A value of clearmode=2 will save you some millisecond(s) of time by leaving the drawing surface in an undefined state - and leaving the job of reinitializing the surface to you.
    %if trial==1; %synchronize on first trial
    %    [time.fix_flip_begin(trial) time.fix_onset(trial) time.fix_flip_end(trial) Missed Beampos]=...
    %        Screen('Flip', window, time.starttime + time.initial_delay,1); %clearmode = 1;
    %    
    %else
    %    [time.fix_flip_begin(trial) time.fix_onset(trial) time.fix_flip_end(trial) Missed Beampos]=...
    %        Screen('Flip', window, time.ITI_onset(trial-1) + time.ITI(trial-1)-time.ifi,1); %clearmode = 1;
    %end

    %generate grating
    [grating] = genGratings_single(window,stim,trial,data,color,where,trialWithinBlock); 
    %draw stimulus
    Screen('Flip', window);
    Screen('DrawTexture',window,grating.g,[], grating.rect); %orientation and sf
    Screen('DrawingFinished', window);

    
    if ~presExample
        if trialWithinBlock==1; %synchronize on first trial
           [time.stim_flip_begin time.stim_onset time.stim_flip_end Missed Beampos]=...
            Screen('Flip', window, time.starttime + time.initial_delay,clearmode); %clearmode = 0;
        else
           [time.stim_flip_begin time.stim_onset time.stim_flip_end Missed Beampos]=...
            Screen('Flip', window, time.fb_end + time.ITI-time.ifi,clearmode); %clearmode = 0;
         
        end
    else %don't clear screen on practice trials
         [time.stim_flip_begin time.stim_onset time.stim_flip_end Missed Beampos]= ...
                Screen('Flip', window, [],clearmode); %if clearmode = 1, "too slow" appears on the screen. i have no idea why. i cannot find where this text would have been drawn previously
             
    end
    % Real time of onset of second stimulus is returned in 'vbl'
    %The when parameter therefore allows you timing in absolute time units (seconds of system time) or relative to some point in time (starttime + 30 seconds), 
    %or some number of monitor refresh intervals after the last stimulus onset via the formula:
    %when=vbl + (waitframes - 0.5)*ifi 

    %check timing
    %time.observed_ITI(trial) = time.fix_presented-time.end_of_trial;
    %time.time_to_present_fix(trial) = time.fix_presented-time.trial_begin;
    %time.til_stim_cleared(trial) = time.stim_cleared-time.stim_presented;
    %time.between_fix_and_stim(trial) = time.stim_presented-time.fix_presented;
    %if trial>1 ; time.observed_ITI(trial) = time.stim_onset(trial)-time.ITI_onset(trial-1); end


catch
	% ---------- Error Handling ---------- 
	% If there is an error in our code, we will end up here.

	% The try-catch block ensures that Screen will restore the display and return us
	% to the MATLAB prompt even if there is an error in our code.  Without this try-catch
	% block, Screen could still have control of the display when MATLAB throws an error, in
	% which case the user will not see the MATLAB prompt.
    clean_up     
    % Restore preferences
    %Screen('Preference', 'VisualDebugLevel', oldVisualDebugLevel);
    %Screen('Preference', 'SuppressAllWarnings', oldSupressAllWarnings);

	% We throw the error again so the user sees the error description.
	psychrethrow(psychlasterror);
end    
    