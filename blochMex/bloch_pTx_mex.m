%
%	Bloch simulator.
%	[mx,my,mz] = bloch_pTx_mex(b1,gr,tp,t1,t2,df,dp,b0,s,mode,mx,my,mz)
%
%	Bloch simulation of rotations due to B1, gradient and
%	off-resonance, including relaxation effects.  At each time
%	point, the rotation matrix and decay matrix are calculated.
%	Simulation can simulate the steady-state if the sequence
%	is applied repeatedly, or the magnetization starting at m0.
%
%	INPUT:
%		b1 = (nTx x M) RF pulse in V.  Can be complex. 
%		gr = (M x 1,2,or 3) 1,2 or 3-dimensional gradient in G/cm.
%		tp = (M x 1) time duration of each b1 and gr point, in seconds,
%				or 1x1 time step if constant for all points
%				or monotonically INCREASING endtime of each
%				interval..
%		t1 = T1 relaxation time in seconds. A negative value disables
%		        simulation of T1/T2 decay.
%		t2 = T2 relaxation time in seconds. A negative value disables
%		        simulation of T1/T2 decay.
%		df = (N x 1) Array of off-resonance frequencies (Hz)
%		dp = (P x 1,2,or 3) Array of spatial positions (cm).  
%			    Width should match width of gr.
%       b0 = (P x 1) Array of B0 off-resonance per voxel (Hz)
%       s  = (nTx x P) Array of Tx coil sensitivities (G/V)
%               if empty/invalid size we're assuming unity
%		mode= Bitmask mode:
%			    Bit 0:  0-Simulate from start or M0, 1-Steady State
%               Bit 1:  1-Record m at time points.  0-just end time.
%
%	(optional)
%		mx,my,mz (PxN) arrays of starting magnetization, where N
%	            is the number of frequencies and P is the number
%               of spatial positions.
%
%	OUTPUT:
%		mx,my,mz = PxN arrays of the resulting magnetization
%				components at each position and frequency.
%
%	B. Hargreaves.	Nov 2003.
%   Updated and extended for pTx:
%   D. Bosch.       Feb 2021.
%

function [mx,my,mz] = bloch_pTx_mex(b1,gr,tp,t1,t2,df,dp,b0,s,mode,mx,my,mz)
    % If the user reaches this it means the mex file was not found.
    error('Check if the mex file was compiled!');
end




