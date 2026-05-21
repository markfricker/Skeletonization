function bw = wsNmsSkeletonize(im, orf, p)
%WSNMSSKELETONIZE  Watershed ridge skeleton combined with Kovesi NMS.
%
%   bw = wsNmsSkeletonize(im, orf, p)
%
% Generates a closed-loop skeleton via 8-connected watershed, then uses
% Kovesi's featureorient + smoothorient + nonmaxsup pipeline to detect
% additional open-ended tubule ridges, and combines both into a final
% binary mask.
%
% Designed for ER network images where tubules form both closed loops and
% free-ending branches.  The watershed handles loops; NMS captures branches.
%
% IMPORTANT: the caller (erSkeletonRun) should apply h-minima pre-processing
% and the cell-boundary intensity barrier to im before calling this function.
%
% INPUTS
%   im  – 2-D single enhanced image, normalised [0, 1].
%   orf – ignored; orientation is computed internally from im via featureorient.
%         Accepted for dispatcher interface consistency.
%   p   – parameter struct (unused; accepted for dispatcher interface
%         consistency).
%
% OUTPUT
%   bw – logical binary mask combining watershed ridge lines and NMS
%        detections.  skeletonPostProcess will apply bwskel + area filter.
%
% DEPENDENCIES (must be on MATLAB path)
%   Kovesi phase congruency toolkit: featureorient, smoothorient, nonmaxsup

im  = double(im);
WS  = watershed(im, 8);
sk1 = WS == 0;

% Merge WS ridge lines back into image so NMS also fires along them.
im2 = max(im, double(sk1));

% Kovesi orientation estimate (normal-to-ridge, degrees) + smooth + NMS.
or  = featureorient(im2, 0, 1, 3, 0);
or  = smoothorient(or, 1.5);
nms = nonmaxsup(im2, or, 3);

% Union of WS and NMS; keep largest component; fill diagonal gaps.
tmp = bwareafilt(sk1 | logical(nms), 1);
bw  = bwmorph(tmp, 'majority') | tmp;
end
