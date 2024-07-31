%
% This software is released under the GPL v3. It is provided AS-IS and no
% warranty is given.
%
% Author: Ye Zhiwei, 2024

function [ Uni_data ] = DataUnitUnify( data, unittype, pixelsize )
%DATAUNITNORM Summary of this function goes here
%   unittype: 1, um; 2, nm; 3, pixel
switch unittype
    case 1
        Uni_data = data .* 1000;
    case 2
        Uni_data = data;
    case 3
        Uni_data = data .* pixelsize;
    case 4
        Uni_data = data .* 100 .* pixelsize;
end
end

