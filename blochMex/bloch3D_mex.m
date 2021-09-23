function [mag,data] = bloch3D_mex(b0,dist,pulse,b1maps,g,ts,mag,target)
    if ~exist('target','var')
        target = [];
    end
% ==== Show simulation of M(time,position,freq)
%
%	-> Try changing mode from 3 to 2...
%

%clear;
%close all;

% try mode=1 for steady-state
mode=0;
b1 = pulse.';
% b1maps from T/V to G/V
s = b1maps.'*1e4;
% grads from T/m to G/cm
g = g.'*1e4/1e2;
% pos from m to cm
x = dist.' * 100;

% only 1 freq
f = 0;
nf = size(f,2);

t = ts;
% t1=1.0;
% t2=0.2;
t1=0;
t2=0;

[mx,my,mz] = bloch_pTx_mex(b1*1,g,t,t1,t2,f,x,b0,s,mode);

mag = [mx.';my.';mz.'];
end
