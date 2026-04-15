function bwSkel = erSkeletonBuild(imIn, cellBoundary, p)
%ERSKELETONBUILD  Build ER network skeleton using ER-specific methods.
%
%   bwSkel = erSkeletonBuild(imIn, cellBoundary, p)
%
% Constructs a single-pixel-wide skeleton from an enhanced ER image using
% one of four methods.  Applies optional h-minima pre-processing and uses
% the cell boundary as a barrier / clipping mask.
%
% METHODS
%   'watershed'   – standard 8-connected watershed; zero-pixels form skeleton
%   'hysteresis'  – dual-threshold with connectivity propagation (Canny-like)
%   'WS + NMS'    – watershed combined with non-maxima suppression
%   'hybrid'      – hysteresis + watershed loop-handling; reconnects free
%                   branches to the loop skeleton via Bresenham lines
%
% REQUIREMENTS
%   'WS + NMS' requires featureorient, smoothorient, nonmaxsup (Kovesi toolkit).
%   'hybrid'   requires bresenham (line-drawing utility).
%   Both must be on the MATLAB path.
%
% INPUTS
%   imIn         – [nY nX nC nZ nT] single/double enhanced image, values in
%                  [0,1].
%   cellBoundary – [nY nX ...] logical cell-interior mask (1 = inside cell),
%                  or [] to skip boundary handling.  May have fewer C/Z/T
%                  dimensions than imIn (broadcast via min indexing).
%   p            – parameter struct (typically app.parameters.er.network.skeleton):
%                    p.method            – method string (see above)
%                    p.hminUse           – logical; apply h-minima pre-processing
%                    p.hmin              – h-minima threshold
%                    p.thresholdMax      – upper hysteresis threshold
%                    p.hysteresis.threshold – lower threshold for 'hysteresis'
%                    p.hybrid.threshold     – lower threshold for 'hybrid'
%
% OUTPUT
%   bwSkel – [nY nX nC nZ nT] logical skeleton

[nY, nX, nC, nZ, nT] = size(imIn);

hasBoundary = ~isempty(cellBoundary);
if hasBoundary
    [~, ~, bC, bZ, bT] = size(cellBoundary);
end

bwSkel = false([nY, nX, nC, nZ, nT]);

for iT = 1:nT
    for iZ = 1:nZ
        for iC = 1:nC
            im = double(imIn(:,:,iC,iZ,iT));

            % outer region = pixels NOT inside the cell
            if hasBoundary
                outerRegion = ~cellBoundary(:,:, min(iC,bC), min(iZ,bZ), min(iT,bT));
            else
                outerRegion = false(nY, nX);
            end

            % impose perimeter as intensity barrier for watershed-based methods
            switch p.method
                case {'watershed','WS + NMS','hybrid'}
                    im(bwperim(outerRegion, 8)) = 1;
            end

            % optional h-minima to suppress spurious local minima
            if p.hminUse
                imexmin = imextendedmin(im, p.hmin);
                im      = imimposemin(im, imexmin);
                im(isinf(im)) = 0;
            end

            % ---- method-specific skeleton extraction --------------------
            switch p.method

                case 'watershed'
                    WS = watershed(im, 8);
                    bw = WS == 0;

                case 'hysteresis'
                    loT = p.hysteresis.threshold;
                    hiT = p.thresholdMax;
                    aboveLo = im > loT;
                    [hiR, hiC] = find(im > hiT);
                    bw = bwselect(aboveLo, hiC, hiR, 8);

                case 'WS + NMS'
                    WS  = watershed(im, 8);
                    sk1 = WS == 0;
                    im  = max(im, single(sk1));
                    or  = featureorient(im, 0, 1, 3, 0);
                    or  = smoothorient(or, 1.5);
                    nms = nonmaxsup(im, or, 3);
                    tmp = bwareafilt(sk1 | logical(nms), 1);
                    tmp = bwmorph(tmp, 'majority') | tmp;
                    bw  = bwskel(tmp);

                case 'hybrid'
                    bw = localHybrid(im, outerRegion, p);

                otherwise
                    error('erSkeletonBuild:unknownMethod', ...
                        ['Unknown method "%s". ' ...
                         'Use ''watershed'', ''hysteresis'', ''WS + NMS'' or ''hybrid''.'], ...
                        p.method);
            end

            % ---- shared post-processing ---------------------------------
            bw = bwskel(logical(bw));
            bw = bwmorph(bw, 'fill');
            bw = bwskel(bw);
            bw(outerRegion) = false;

            bwSkel(:,:,iC,iZ,iT) = bw;
        end
    end
end
end


% =========================================================================
function bw = localHybrid(im, outerRegion, p)
%LOCALHYBRID  Hybrid hysteresis+watershed skeleton for loop-containing ER.

loT = p.hybrid.threshold;
hiT = p.thresholdMax;

% --- hysteresis skeleton -------------------------------------------------
aboveLo = im > loT;
[hiR, hiC] = find(im >= hiT);
bw = bwselect(aboveLo, hiC, hiR, 8);

% --- watershed skeleton --------------------------------------------------
WS            = watershed(im, 8);
skLoopInitial = WS == 0;
skLoopInitial = bwmorph(skLoopInitial, 'thin', Inf);
skLoop        = skLoopInitial;

% identify isolated rings (components not part of the largest component)
skRing = skLoop & ~bwareafilt(skLoop, 1);

if any(skRing(:))
    skLoopFill   = imfill(skRing, 'holes');
    skLoopFill   = skLoopFill & ~skRing;          % interior of loops
    skLoopPoints = bwulterode(skLoopFill);        % medial point in each loop
    bwLoop       = imfill(bw | skRing, find(skLoopPoints));
    skHyst       = bwmorph(bwLoop, 'thin', Inf);
else
    skHyst = bwmorph(bw, 'thin', Inf);
end

% --- extract tree branches from thinned skeleton -------------------------
W1     = watershed(skHyst, 4);
skLoop2 = W1 == 0;
skTree = xor(skHyst, skLoop2);
skTree(skLoop) = false;

if any(skRing(:))
    skTree = skTree & ~skLoopFill;
    skTree = skTree | skRing;
end

% --- reconnect free branches to loop skeleton via Bresenham lines --------
epsk = bwmorph(skHyst, 'endpoints');
[r, c] = find(epsk);
skTree = bwselect(skTree, c, r);

epskTree = bwmorph(skTree, 'endpoints');
epskTree = xor(epsk, epskTree);

connected = bwareafilt(epskTree | skLoop, 1);
epskTree  = epskTree & ~connected;

[~, skW_idx]  = bwdist(skLoopInitial);
[y1, x1]      = ind2sub(size(skRing), skW_idx(epskTree));
[y2, x2]      = find(epskTree);

lineFn = @(xa,ya,xb,yb) bresenham(xa,ya,xb,yb);
[xs, ys] = arrayfun(lineFn, x1, y1, x2, y2, 'UniformOutput', false);
pixIdx  = cellfun(@(xi,yi) sub2ind(size(skLoop), yi, xi), xs, ys, 'UniformOutput', false);
if ~isempty(pixIdx)
    skTree(cat(1, pixIdx{:})) = true;
end

% --- merge and clean -----------------------------------------------------
if any(skRing(:))
    skTree(skLoopFill) = false;
    skLoop(skRing)     = false;
end

skFinal = skTree | skLoop;
skFinal = bwmorph(skFinal, 'thin',  Inf);
skFinal = bwmorph(skFinal, 'diag');
skFinal = bwmorph(skFinal, 'fill');
skFinal = bwmorph(skFinal, 'thin', Inf);

bw = skFinal;
end
