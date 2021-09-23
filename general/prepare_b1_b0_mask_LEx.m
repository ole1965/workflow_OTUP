function [mapB1,mapB0,mask,spatial] = prepare_b1_b0_mask_LEx(loadpath)
%function to load and postprocess the following data for a subject:
%    - design mask for the entire head, including the subcutaneous tissue (inportant for LEx)
%    - B1 map
%    - B0 map
%    - spatial location of each voxel

%load data
load(loadpath);

%Rotate everything whats necessay to have it in the same orientation as on the scanner computer
deltaB0_Hz = rot90(deltaB0_Hz);
if exist('tissueMask','var')
    mask = tissueMask;
elseif exist('tissuemask','var')
    mask = tissuemask;
end
mask = rot90(mask);
nT_per_V = rot90(nT_per_V);
meta.matSizeOld = meta.matSize;
meta.matSize = [meta.matSize(2);meta.matSize(1);meta.matSize(3)];
meta.FOV = [meta.FOV(2);meta.FOV(1);meta.FOV(3)]; 
meta.voxSz = [meta.voxSz(2);meta.voxSz(1);meta.voxSz(3)];
meta.isocenter = [meta.matSizeOld(2)-meta.isocenter(2)+1;meta.isocenter(1);meta.isocenter(3)];%be carefull, orientation of x axis changes because of rotation

%clean and reshape the data
deltaB0_Hz(isnan(deltaB0_Hz)) = 0;
nT_per_V(isnan(nT_per_V)) = 0;

mapB0 = -deltaB0_Hz(:);
mask = mask(:);

nCh = size(nT_per_V,4);
mapB1 = nT_per_V*1e-9; %convert from nT to T
mapB1 = reshape(mapB1,[meta.matSize(1)*meta.matSize(2)*meta.matSize(3) nCh]); %reshape to a 2D matrix with one column per transmit channel

%create spatial location matrix
x = zeros(1,meta.matSize(1));
y = zeros(1,meta.matSize(2));
z = zeros(1,meta.matSize(3));

for i = meta.isocenter(1)-1:-1:1
    x(i) = x(i+1) - meta.voxSz(1);
end
for i = meta.isocenter(1)+1:1:meta.matSize(1)
    x(i) = x(i-1) + meta.voxSz(1);
end
for i = meta.isocenter(2)-1:-1:1
    y(i) = y(i+1) - meta.voxSz(2);
end
for i = meta.isocenter(2)+1:1:meta.matSize(2)
    y(i) = y(i-1) + meta.voxSz(2);
end
for i = meta.isocenter(3)-1:-1:1
    z(i) = z(i+1) - meta.voxSz(3);
end
for i = meta.isocenter(3)+1:1:meta.matSize(3)
    z(i) = z(i-1) + meta.voxSz(3);
end
[axisX3D,axisY3D,axisZ3D]=meshgrid(y,x,z);
dpDimZ3D=reshape(axisZ3D,[],1);  dpDimXY3D=reshape(axisX3D+1i*axisY3D,[],1); 
spatial = [dpDimXY3D,dpDimZ3D];  

end
