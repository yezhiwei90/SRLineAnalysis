%
% This software is released under the GPL v3. It is provided AS-IS and no
% warranty is given.
%
% Author: Ye Zhiwei, 2024

function [ bin_x, bin_y ] = Binning_WideFieldData(widefield_x, normalized_widefield_y, magnification)
%BINNING_WIDEFIELDDATA Summary of this function goes here
%   Detailed explanation goes here
y = normalized_widefield_y;
x = widefield_x;
mag = magnification;

bin_x = binningdata(x, mag);
bin_y = binningdata(y, mag);

end

function bin_data = binningdata(data, magnification)
y = data;
mag = magnification;
j=0;
temp=0;
while j*mag<length(y)-mag
    num = fix((j+1)*mag)-fix(j*mag);
    tempmax=sum(y((1+fix(j*mag)):(num+fix(j*mag)) ));
    j=j+1;
    temp(j)=tempmax./num;
end
tempmax=0;
ii=round(j*mag)+1;
k=0;
while ii<=length(y)
    tempmax=tempmax+y(ii);ii=ii+1;k=k+1;
end
if k >= mag/2
    temp(j+1)=tempmax./k;
else 
end

bin_data=temp;
end
