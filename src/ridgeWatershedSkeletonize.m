function bw = ridgeWatershedSkeletonize(im, p)
%RIDGEWATERSHEDSKELETONIZE  8-connected watershed ridge lines as skeleton.
%
%   bw = ridgeWatershedSkeletonize(im, p)
%
% The zero-labelled pixels of the 8-connected watershed on im form a
% single-pixel-wide ridge skeleton.  Because every intensity ridge becomes
% a watershed divide, this method naturally captures closed loops as well
% as free-ending branches.
%
% IMPORTANT: the caller (erSkeletonRun) is responsible for applying any
% h-minima pre-processing and imposing a cell-boundary intensity barrier
% before passing im to this function.  Without an intensity barrier at the
% cell edge the watershed basins will span the entire image and the ridge
% lines will not be confined to the cell interior.
%
% INPUTS
%   im – 2-D single image normalised [0, 1].
%   p  – parameter struct (unused; accepted for dispatcher interface
%        consistency).
%
% OUTPUT
%   bw – logical, same size as im.  True where watershed label == 0.

bw = watershed(double(im), 8) == 0;
end
