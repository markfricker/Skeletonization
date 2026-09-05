function bw = wsOnlySkeletonize(im, p)
%WSONLYSKELETONIZE  Pure 8-connected watershed ridge skeleton, no NMS.
%
%   bw = wsOnlySkeletonize(im, p)
%
% Generates a closed-loop skeleton via 8-connected watershed only -- no
% NMS branch-detection union (see wsNmsSkeletonize.m for that combined
% variant, which runs this exact same watershed step internally as its
% own first stage before unioning in NMS).
%
% Exact port of AnalyzER_v2's 'watershed' method (sk = watershed(im,8)==0,
% AnalyzER_v2.m's fnc_process_skeleton). h-minima suppression is not a
% separate algorithm branch here -- per this dispatcher family's own
% convention (see funcSkeletonDispatcher.m's header comment), it is
% applied by the CALLER before this function is invoked. So this method
% with hminUse enabled reproduces v2's 'WS + hmin' method exactly, and
% with hminUse disabled reproduces v2's plain 'watershed' method exactly.
%
% NOTE: not to be confused with watershedSkeletonize.m (a different,
% pre-existing function -- marker-controlled watershed on a distance
% transform, backing the dispatcher's 'distWatershed' method). This is a
% deliberately distinctly-named sibling, not a replacement.
%
% IMPORTANT: the caller (erSkeletonRun) should apply h-minima
% pre-processing and the cell-boundary intensity barrier to im before
% calling this function -- same convention as every other method in this
% dispatcher family.
%
% INPUTS
%   im - 2-D single enhanced image, normalised [0, 1].
%   p  - parameter struct (unused; accepted for dispatcher interface
%        consistency, matching ridgeWatershed/wsNmsSkeletonize's own
%        signatures).
%
% OUTPUT
%   bw - logical binary mask from the watershed ridge lines.
%        skeletonPostProcess (called by funcSkeletonDispatcher) applies
%        bwskel + area filter on top of this.

im = double(im);
WS = watershed(im, 8);
bw = WS == 0;
end
