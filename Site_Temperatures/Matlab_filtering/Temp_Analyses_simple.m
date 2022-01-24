% Script for bandpass filtering Palau Temperature Time Series datasets
% The imported data had reformatted using a combination of excel and
% R in order to get everything into one long vector form with serial
% datetime values from the logger output files.
% Leap years are accounted for. Missing values within the time series have been retained.
% Objective:
% Filter the time series to keep mainly daily signal (i.e. remove frequencies
% lower than 36 hours and higher than 5 hours, meaning you retain patterns that cycle between 5 and 36 hours)
% This will remove seasonal and sporadic variability so that we can compare daily range across sites with data
% from various times.

%% Import the Data
DropOff=importdata('DropOff_Matlab.txt');
Ngerd=importdata('Ngerdiluches_Matlab.txt');
Helen=importdata('Helen_Matlab.txt');
Kayangel=importdata('Kayan_Matlab.txt');
Mecherchar=importdata('Mecherchar_Matlab.txt');
Ngerchelong=importdata('Ngerchelong_Matlab.txt');
Risong=importdata('Risong_Matlab.txt');
Taoch=importdata('Taoch_Matlab.txt');
Ngermid=importdata('Ngermid_Matlab.txt');

%% Filter Time Series
fs=1/(30*60);      % sampling frequency in Hertz % for half hour sample
fc=1/(30*60*60);   % cutoff frequency for ~1 day cut off frequency
fc2=1/(5*60*60);   % cutoff frequency for > 5 hour freq

% Filter frequencies
Wn=fc/(fs/2);
Wn2=fc2/(fs/2);
[b1,a1]=butter(4, [Wn,Wn2], 'bandpass'); % high in order to keep the higher freq stuff , % low in order to keep low freq stuff, % bandpass middle

filt_DropOff=filtfilt(b1,a1, DropOff(:,2));
filt_Ngerd=filtfilt(b1,a1, Ngerd(:,2));
filt_Kayangel=filtfilt(b1,a1, Kayangel(:,2));
filt_Ngerchelong=filtfilt(b1,a1, Ngerchelong(:,2));
filt_Helen=filtfilt(b1,a1, Helen(:,2));
filt_Mecherchar=filtfilt(b1,a1, Mecherchar(:,2));

%% For data with 15 minute sampling intervals
cfs=1/(15*60);      % sampling frequency in Hertz % for half hour sample

% Filter frequencies
cWn=fc/(cfs/2);
cWn2=fc2/(cfs/2);
[b2,a2]=butter(4, [cWn,cWn2], 'bandpass'); % high in order to keep the higher freq stuff , % low in order to keep low freq stuff, % bandpass middle

filt_Risong=filtfilt(b2,a2, Risong(:,2));
filt_Taoch=filtfilt(b2,a2, Taoch(:,2));
filt_Ngermid=filtfilt(b2,a2, Ngermid(:,2));

%% Save the vectors to for stat analyses in R
save DO_filtered_temp.txt filt_DropOff -ascii
save Ngerdiluches_filtered_temp.txt filt_Ngerd -ascii
save Kayangel_filtered_temp.txt filt_Kayangel -ascii
save Ngerchelong_filtered_temp.txt filt_Ngerchelong -ascii
save Helen_filtered_temp.txt filt_Helen -ascii
save Mecherchar_filtered_temp.txt filt_Mecherchar -ascii
save Risong_filtered_temp.txt filt_Risong -ascii
save Toach_filtered_temp.txt filt_Taoch -ascii
save Ngermid_filtered_temp.txt filt_Ngermid -ascii
