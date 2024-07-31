%
% This software is released under the GPL v3. It is provided AS-IS and no
% warranty is given.
%
% Author: Ye Zhiwei, 2024

function [ fit_optionS ] = StartPointInitializer( xs, nys, gauss_num, fit_optionS )
%STARTPOINTINITIALIZER Summary of this function goes here
%   Detailed explanation goes here
limit = 0.3; rectime = 0;
[pos, PeakLoc, flag] = FindPeaks(nys, gauss_num, limit, rectime);
if flag
    Peaknum = numel(pos) + 1;
    pos = [0;pos;numel(PeakLoc)];
    StartPeakPos = zeros(1, Peaknum);
    for i = 1: Peaknum
    StartPeakPos(i) = round( mean(pos(i)+1:pos(i+1)) );
    CrPeakLoc = PeakLoc(StartPeakPos(i));
    fit_optionS.StartPoint(3*i - 2:3*i-1) = [nys(CrPeakLoc), xs(CrPeakLoc )];
    end
end
end

function [pos,PeakLoc, flag] = FindPeaks(nys, Targetnum, limit, rectime)
PeakLoc = find(nys > limit);
pos = find(abs(diff(PeakLoc)) > 1);
Peaknum = numel(pos) + 1; flag = false;
if Peaknum > Targetnum
    limit = min([1, limit + 0.1]);
    if rectime < 20
        [pos,PeakLoc, flag] = FindPeaks(nys, Targetnum, limit, rectime);
    end
else
    if Peaknum < Targetnum
        limit = min([0.1, limit - 0.1]);
        if rectime < 20
            [pos,PeakLoc, flag] = FindPeaks(nys, Targetnum, limit, rectime);
        end
    else
        flag = true;
    end
end
end