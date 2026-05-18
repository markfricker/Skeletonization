function bwOut = skeletonJunctionRepair(bwIn, imRef, p)
%SKELETONJUNCTIONREPAIR  Reconnect skeleton endpoints at tubule junctions.
%
%   bwOut = skeletonJunctionRepair(bwIn, imRef, p)
%
% Two-stage repair:
%
%   Stage A — connect to existing edges
%     For each degree-1 endpoint, search for any skeleton pixel from a
%     different connected component within repairRadius.  Connects the
%     endpoint to the best candidate: the path with highest mean intensity
%     that also clears the repairMinIntensity floor on every pixel.
%     After each connection the two component labels are merged in compL so
%     that subsequent searches treat the newly joined pixels as one
%     component, preventing cycles.
%
%   Stage B — junction node for remaining endpoints
%     Endpoints not resolved by Stage A are clustered by proximity.  For
%     each cluster of >= 2 endpoints, a single junction node is chosen
%     (the intensity maximum within the search region) and all endpoints
%     connect to that one point rather than to each other, avoiding the
%     triangle artefact that pairwise reconnection produces.
%     Each path is also subject to the repairMinIntensity gate and a
%     same-component guard (skips the connection if endpoint and junction
%     node already belong to the same component).
%
% INPUTS
%   bwIn   – 2-D logical skeleton.
%   imRef  – 2-D single intensity image (same size as bwIn), used to assess
%            path quality and to locate the junction node.  Typically the
%            h-min + boundary-barrier preprocessed image from erSkeletonRun.
%   p      – struct with optional fields:
%              .repairRadius       – endpoint search radius (pixels).
%                                   ~2–3× tubule FWHM.  Default: 10.
%              .repairMinIntensity – minimum pixel intensity permitted on
%                                   any proposed connection path [0,1].
%                                   Paths with any pixel below this are
%                                   rejected.  Default: 0.1.
%              .repairMinCluster   – minimum endpoints per cluster for
%                                   Stage B.  Default: 2.
%
% OUTPUT
%   bwOut – logical skeleton with junction gaps filled.

gapRadius        = double(getf(p, 'repairRadius',       10));
minPathIntensity = double(getf(p, 'repairMinIntensity', 0.1));
minCluster       = double(getf(p, 'repairMinCluster',   2));

bwOut = bwIn;
if ~any(bwIn(:)), return; end

[nY, nX] = size(bwIn);

% Label connected components of the input skeleton
[compL, ~] = bwlabel(bwIn, 8);

% Degree-1 endpoints
ep = bwmorph(bwIn, 'endpoints');
[epR, epC] = find(ep);
nEp = numel(epR);
if nEp == 0, return; end

connected = false(nEp, 1);   % tracks which endpoints Stage A resolved
anyAdded  = false;

% =========================================================================
% Stage A: connect each endpoint to the nearest skeleton pixel from a
%          different connected component (joins free ends to intact edges)
% =========================================================================
for i = 1:nEp
    r = epR(i);  c = epC(i);
    myLabel = compL(r, c);

    % Search window (clamped to image bounds)
    r0 = max(1,  r - gapRadius);  r1 = min(nY, r + gapRadius);
    c0 = max(1,  c - gapRadius);  c1 = min(nX, c + gapRadius);

    winBw = bwIn(r0:r1, c0:c1);
    winL  = compL(r0:r1, c0:c1);

    % Skeleton pixels belonging to any other component
    otherMask = winBw & (winL ~= myLabel) & (winL > 0);
    if ~any(otherMask(:)), continue; end

    [oRl, oCl] = find(otherMask);
    oRg = oRl + r0 - 1;
    oCg = oCl + c0 - 1;

    % Rank candidates by proximity; evaluate up to 10 nearest
    dist2    = (oRg - r).^2 + (oCg - c).^2;
    [~, ord] = sort(dist2);
    nCand    = min(10, numel(ord));

    bestScore = -1;  bestR = NaN;  bestC = NaN;
    for j = 1:nCand
        ti  = ord(j);
        pxi = bresenhamPixels(r, c, oRg(ti), oCg(ti), nY, nX);
        pathInt = double(imRef(pxi));
        if min(pathInt) >= minPathIntensity
            score = mean(pathInt);
            if score > bestScore
                bestScore = score;
                bestR = oRg(ti);
                bestC = oCg(ti);
            end
        end
    end

    if ~isnan(bestR)
        bwOut(bresenhamPixels(r, c, bestR, bestC, nY, nX)) = true;
        connected(i) = true;
        anyAdded = true;
        % Merge component labels so later iterations don't reconnect pixels
        % that are now in the same component (prevents cycle creation).
        mergedLabel = compL(bestR, bestC);
        if mergedLabel ~= myLabel
            compL(compL == mergedLabel) = myLabel;
        end
    end
end

% =========================================================================
% Stage B: cluster unresolved endpoints and connect through a single
%          intensity-guided junction node (avoids triangle artefact)
% =========================================================================
remIdx = find(~connected);
if numel(remIdx) >= minCluster

    % Recompute connected components from the updated skeleton so that
    % Stage B checks reflect any connections already made in Stage A.
    [compL, ~] = bwlabel(bwOut, 8);

    % Build endpoint image for remaining endpoints only
    epImg = false(nY, nX);
    for i = 1:numel(remIdx)
        epImg(epR(remIdx(i)), epC(remIdx(i))) = true;
    end

    se = strel('disk', max(1, ceil(gapRadius / 2)));
    [Lep, nClusters] = bwlabel(imdilate(epImg, se), 8);

    for k = 1:nClusters
        clMask    = (Lep == k) & epImg;
        [cR, cC]  = find(clMask);
        if numel(cR) < minCluster, continue; end

        % Junction node: intensity maximum in a circle around the centroid
        centR = mean(cR);  centC = mean(cC);
        r0 = max(1,  floor(centR - gapRadius));
        r1 = min(nY, ceil (centR + gapRadius));
        c0 = max(1,  floor(centC - gapRadius));
        c1 = min(nX, ceil (centC + gapRadius));

        roi = double(imRef(r0:r1, c0:c1));
        [rg, cg] = ndgrid(r0:r1, c0:c1);
        roi((rg - centR).^2 + (cg - centC).^2 > gapRadius^2) = -Inf;

        [~, idx]    = max(roi(:));
        [jrL, jcL] = ind2sub(size(roi), idx);
        jR = r0 + jrL - 1;
        jC = c0 + jcL - 1;

        % Connect each endpoint to the junction node, subject to intensity
        % gate and same-component guard (prevents cycles).
        jLabel = compL(jR, jC);
        for i = 1:numel(cR)
            if cR(i) == jR && cC(i) == jC, continue; end
            % Skip if already in the same connected component as the node.
            if compL(cR(i), cC(i)) == jLabel && jLabel > 0, continue; end
            pxi     = bresenhamPixels(cR(i), cC(i), jR, jC, nY, nX);
            pathInt = double(imRef(pxi));
            if min(pathInt) >= minPathIntensity
                bwOut(pxi) = true;
                anyAdded   = true;
                % Update label so subsequent endpoints in this cluster
                % don't create a cycle through the junction node.
                epLabel = compL(cR(i), cC(i));
                if jLabel == 0
                    jLabel = epLabel;
                elseif epLabel ~= jLabel && epLabel > 0
                    compL(compL == epLabel) = jLabel;
                end
            end
        end
    end
end

% Re-skeletonise to thin any pixels thickened by converging connections
if anyAdded
    bwOut = bwskel(bwOut);
end
end


% =========================================================================
function pixIdx = bresenhamPixels(r1, c1, r2, c2, nY, nX)
%BRESENHAMPIXELS  Linear indices of pixels on the integer Bresenham line.

dr  = abs(r2 - r1);  dc  = abs(c2 - c1);
sr  = sign(r2 - r1); sc  = sign(c2 - c1);
r   = r1;            c   = c1;
err = dr - dc;

pixIdx = zeros(dr + dc + 1, 1);
n      = 0;

while true
    if r >= 1 && r <= nY && c >= 1 && c <= nX
        n = n + 1;
        pixIdx(n) = r + (c - 1) * nY;
    end
    if r == r2 && c == c2, break; end
    e2 = 2 * err;
    if e2 > -dc,  err = err - dc;  r = r + sr;  end
    if e2 <  dr,  err = err + dr;  c = c + sc;  end
end

pixIdx = pixIdx(1:n);
end


% =========================================================================
function v = getf(s, field, default)
    if isfield(s, field) && ~isempty(s.(field))
        v = s.(field);
    else
        v = default;
    end
end
