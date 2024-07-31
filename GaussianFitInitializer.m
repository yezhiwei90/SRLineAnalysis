%
% This software is released under the GPL v3. It is provided AS-IS and no
% warranty is given.
%
% Author: Ye Zhiwei, 2024

function [ fit_type,fit_option ] = GaussianFitInitializer( length, gaussian_num )
%GAUSSIANFITINITIALIZER Summary of this function goes here
%   Detailed explanation goes here
num = gaussian_num;
if num < 8
fittypename = ['gauss' num2str(num)];
else
    errordlg('Input gaussian number exceeds maximun limit!');
    return;
end
fto = fitoptions(fittypename);
ftt = fittype(fittypename);
lowerlimits = [0 0 0];
highlimits = [1 length length*2];
initialvalue = [];
for i = 1 : num
    initialvalue = [initialvalue,1, length/num*(i-1), length/num/2];
end
lowerlimits = repmat(lowerlimits, 1, num);
highlimits = repmat(highlimits, 1, num);
fto.Lower = lowerlimits;
fto.Upper = highlimits;
fto.StartPoint = initialvalue;
fit_type = ftt;
fit_option = fto;
end

