%use radarClass to create an object
radarObj = radarClass(1,"resting",2000,tfm_ecg,radar_dist);

DS=radarObj.DownSampleRadar(200);

radarObj.HrFilter();
radarObj.RrFilter();
%% checks

radarObj.FindPeaks();
radarObj.FindRates();
HrB4Corr=radarObj.HrEst;
RrB4Corr=radarObj.RrEst;
additions = radarObj.findMissingBeats();
removals=radarObj.clearFalsePos();
[missed,excess]=radarObj.CorrelatePeaks();
additions = radarObj.findMissingBeats();
radarObj.FindRates();
HrCorr=radarObj.HrEstFinal;
newTime= radarObj.vTimeNew(removals==0);
newTime = 200*newTime(1:length(newTime)-1);
figure;hold on;

plot(HrB4Corr,'DisplayName',"before")
plot(HrCorr,'DisplayName',"after")
plot(radarObj.GtEst,'DisplayName',"GT");
 legend('show');
 hold off;


%{
check_Class: just a playground to check things, added to git
Main:
made 1 instance of filter design, reused on each patient
made the main run on the class radarClass. 
deleted some old code that we commented

radarClass:
added an mse and mae function
this function also makes mse2HrRaw and mse2HrFinal-
each is two nx2 vectors that show GT heart rate, to MSE error. We can plot this as
alon requested.

added fields for HR before and after algorithms, HrEstFinal and HrPeaksFinal
now data is kept for both before and after our compensation.
CalcError gives error for both before and after compensation.

this push is in the mainClass branch

Left %TODO with reminders on line  147 in main

%}