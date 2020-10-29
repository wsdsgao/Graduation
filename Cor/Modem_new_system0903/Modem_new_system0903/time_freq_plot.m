clear;clc
%load('C:\Users\zj\Documents\MATLAB\12.xls');
load ('signal_trans.mat');
X=resample(signal_trans,320,1024);%Êä³öµÚ2ÐÐ
%E=E';
%n=size(X);
%s=E(1:n(2));
T=0:1/320:length(X)/320;
%T = 0:0.001:2;
%X = chirp(T,0,1,150);
[S,F,T,P] = spectrogram(X,256,255,256,320);
%surf(T,F,10*log10(P),'edgecolor','none'); axis tight;
mesh(T,F,10*log10(P)); axis tight;
colorbar;  
view(0,90);
xlabel('Time (Seconds)'); ylabel('MHz');