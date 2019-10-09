function correctkey(window,color,where,stim,feedback,trial,trialWithinBlock)
%correctkey(window,color,where,stim,trial)
%
%tells subject which category label they should have selected
%
%input:
%window - window pointer
%color - structure containing color information
%sounds - structure containing wrongkey, ching - played if mean(data)<rt_crit, and buzz - played if mean(data)>= rt_crit
%where - structure containing location information
%
%5/22/05        swe     
%12/3/06    swe     modified to work with OpenGL-Psychtoolbox version 3.0.8 on OSX
%7/20/2010  swe     modified for focalBG - wm experiment to be run @ Berkeley on PC


switch stim.stim(trialWithinBlock,1)
case 1,
    corrLabel = 'Category A was correct';
case 2,
    corrLabel = 'Category B was correct';
case 3,
    corrLabel = 'Category C  was correct';
case 4,
    corrLabel = 'Category D was correct';
otherwise,
    corrLabel = [];
end


[normBoundsRect, offsetBoundsRect]= Screen('TextBounds', window, corrLabel);
Screen('DrawText',window,corrLabel, where.xc-normBoundsRect(3)/2, where.yc+2*feedback.fb_size(1)-normBoundsRect(4)/2, color.textColor);

%Screen('Flip', window,[],clearmode);


