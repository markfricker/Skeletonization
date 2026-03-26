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
% INPUTS
%   bw  – 2-D logical mask.
%   p   – struct with optional fields:
%           p.pruneLength  – integer >= 0.  Spurs shorter than this many
%                            pixels are removed by bwskel(...,
%                            'MinBranchLength', pruneLength).
%                            0 = no pruning (default 10).
%           p.minArea      – integer >= 1.  Isolated fragments smaller
%                            than this area are removed before
%                            skeletonisation (bwareaopen).
%                            1 = no area filter (default 10).
%
% OUTPUT
%   skel – logical, single-pixel-wide skeleton, same size as bw.
%
% NOTE
%   bwskel is idempotent: applying it to an already-thin binary image
%   (e.g. the output of NMS or hessianCenterline) is safe and returns
%   the image unchanged except for spur pruning.

pruneLength = getf(p, 'pruneLength', 10);
minArea     = getf(p, 'minArea',     10);

% --- Area filter ---------------------------------------------------------
if minArea > 1
    bw = bwareaopen(logical(bw), minArea);
end

% --- Skeletonise ---------------------------------------------------------
if pruneLength > 0
    skel = bwskel(logical(bw), 'MinBranchLength', pruneLength);
else
    skel = bwskel(logical(bw));
end
end

% ---- local helper -------------------------------------------------------
function v = getf(s, field, default)
    if isfield(s, field) && ~isempty(s.(field))
        v = s.(field);
    else
        v = default;
    end
end
