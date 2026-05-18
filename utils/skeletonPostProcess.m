function skel = skeletonPostProcess(bw, p)
%SKELETONPOSTPROCESS  Shared post-processing for all skeleton methods.
%
%   skel = skeletonPostProcess(bw, p)
%
% Applies area filtering then bwskel to convert a binary foreground mask
% into a single-pixel-wide skeleton.  Called by skeletonSegment and
% funcSkeletonDispatcher after every segmentation method so that
% skeletonization behaviour is identical regardless of which method
% produced the mask.
%
% Spur pruning is intentionally NOT performed here.  Pruning before the
% ER mask is applied would remove branches of interest that happen to
% extend into background regions.  Instead, pruning happens after masking
% inside erSkeletonTruncate, where only mask-boundary stubs remain.
%
% INPUTS
%   bw  – 2-D logical mask.
%   p   – struct with optional fields:
%           p.minArea – integer >= 1.  Isolated fragments smaller than
%                       this area are removed before skeletonisation
%                       (bwareaopen).  1 = no area filter (default 10).
%
% OUTPUT
%   skel – logical, single-pixel-wide skeleton, same size as bw.
%
% NOTE
%   bwskel is idempotent: applying it to an already-thin binary image
%   (e.g. the output of NMS or hessianCenterline) is safe and returns
%   the image unchanged.

minArea = getf(p, 'minArea', 10);

% --- Area filter ---------------------------------------------------------
if minArea > 1
    bw = bwareaopen(logical(bw), minArea);
end

% --- Skeletonise (no pruning — see header) -------------------------------
skel = bwskel(logical(bw));
end

% ---- local helper -------------------------------------------------------
function v = getf(s, field, default)
    if isfield(s, field) && ~isempty(s.(field))
        v = s.(field);
    else
        v = default;
    end
end
