function [bwSkelOut, cisternaeSkeleton, cisternaeNode, perimeterNode, perimeter] = ...
        erSkeletonCisternae(bwSkel, erCisternae, imBackground)
%ERSKELETONCISTERNAE  Construct cisternae nodes and perimeter from ER skeleton.
%
%   [bwSkelOut, cisternaeSkeleton, cisternaeNode, perimeterNode, perimeter] = ...
%       erSkeletonCisternae(bwSkel, erCisternae, imBackground)
%
% Cisternae are large ER regions whose internal structure is not resolved by
% skeletonisation.  This function:
%   1. Extracts the skeleton within each cisterna (cisternaeSkeleton).
%   2. Identifies perimeter nodes — skeleton endpoints at the cisterna edge.
%   3. Places a single central node inside each cisterna (intensity-weighted
%      centroid when imBackground is supplied, geometric centroid otherwise).
%   4. Removes the interior of each cisterna from the skeleton, retaining
%      only the perimeter attachment points.
%
% If erCisternae is empty or all-zero, all outputs except bwSkelOut are
% returned as empty ([]) and the skeleton is passed through unchanged.
%
% INPUTS
%   bwSkel      – [nY nX nC nZ nT] logical skeleton from erSkeletonTruncate.
%   erCisternae – [nY nX nC nZ nT] logical cisternae mask.  Must be the same
%                 size as bwSkel.
%   imBackground – [nY nX nC nZ nT] intensity image used to weight the
%                  cisterna centroid.  Pass [] to use unweighted centroid.
%
% OUTPUTS
%   bwSkelOut        – [nY nX nC nZ nT] logical skeleton with cisternae
%                      removed; perimeter attachment points retained.
%   cisternaeSkeleton – [nY nX nC nZ nT] logical; skeleton pixels that fell
%                       inside the cisternae before removal.
%   cisternaeNode    – [nY nX nC nZ nT] uint32; each cisterna labelled at its
%                      centroid pixel (label = region index).  0 = background.
%   perimeterNode    – [nY nX nC nZ nT] uint32; perimeter attachment pixels
%                      labelled with the index of their parent cisterna.
%   perimeter        – [nY nX nC nZ nT] logical; bwperim of each cisterna,
%                      excluding skeleton pixels.
%
% NOTE
%   The function modifies the skeleton by removing cisterna interiors.  All
%   computations are performed per (iC, iZ, iT) 2-D slice.

% --- early exit if no cisternae -----------------------------------------
if isempty(erCisternae) || ~any(erCisternae, 'all')
    bwSkelOut         = bwSkel;
    cisternaeSkeleton = [];
    cisternaeNode     = [];
    perimeterNode     = [];
    perimeter         = [];
    return
end

[nY, nX, nC, nZ, nT] = size(bwSkel);

bwSkelOut         = bwSkel;
cisternaeSkeleton = false(nY, nX, nC, nZ, nT);
cisternaeNode     = zeros(nY, nX, nC, nZ, nT, 'uint32');
perimeterNode     = zeros(nY, nX, nC, nZ, nT, 'uint32');
perimeter         = false(nY, nX, nC, nZ, nT);

useWeighted = ~isempty(imBackground);

for iT = 1:nT
    for iZ = 1:nZ
        for iC = 1:nC
            sk       = bwSkel(:,:,iC,iZ,iT);
            cisternae = erCisternae(:,:,iC,iZ,iT);

            % ----- extract skeleton within cisternae ---------------------
            skCist = sk & cisternae;

            % perimeter nodes are new endpoints created inside the cisterna
            Pnodes = bwmorph(skCist, 'endpoints');
            % keep only those within one pixel of the cisterna boundary
            Pnodes = Pnodes & ~imerode(cisternae, ones(3));

            % remove cisterna interior from skeleton
            sk = sk & ~cisternae;

            % strip orphaned cisterna pixels that lost skeleton contact
            cisternae = cisternae & ~sk;

            % guarantee overlap at attachment points
            sk        = sk        | Pnodes;
            cisternae = cisternae | Pnodes;

            % also attach any existing skeleton endpoints on the cisterna edge
            edgeEP = bwmorph(sk, 'endpoints') & ~Pnodes & imdilate(cisternae, ones(3));
            Pnodes    = Pnodes    | edgeEP;
            sk        = sk        | edgeEP;
            cisternae = cisternae | edgeEP;

            % split k=2 (degree-2) perimeter nodes by removing them from the
            % skeleton and re-finding endpoints
            otherEP = bwmorph(sk, 'endpoints') & ~Pnodes;
            K2      = Pnodes & ~bwmorph(sk, 'endpoints');
            sk(K2)  = false;

            newEP     = bwmorph(sk, 'endpoints') & ~otherEP;
            Pnodes    = Pnodes    | newEP;
            cisternae = cisternae | newEP;

            % tidy isolated single-pixel artefacts
            sk        = bwmorph(sk,        'clean');
            cisternae = bwmorph(cisternae, 'clean');

            % ----- cisternae central node (intensity-weighted centroid) --
            L = bwlabel(cisternae);
            if useWeighted
                bg    = double(imBackground(:,:,iC,iZ,iT));
                stats = regionprops(L, bg, 'WeightedCentroid');
                centField = 'WeightedCentroid';
            else
                stats = regionprops(L, 'Centroid');
                centField = 'Centroid';
            end

            nodeIm = zeros(nY, nX, 'uint32');
            if ~isempty(stats)
                pos = round(cat(1, stats.(centField)));    % [nRegions × 2] x,y
                idx = sub2ind([nY nX], pos(:,2), pos(:,1));
                for k = 1:numel(idx)
                    nodeIm(idx(k)) = uint32(k);
                end
            end

            % ----- perimeter node labels (parent cisterna index) ---------
            perimIm = uint32(L) .* uint32(Pnodes & sk);

            % ----- store outputs -----------------------------------------
            bwSkelOut(:,:,iC,iZ,iT)         = sk;
            cisternaeSkeleton(:,:,iC,iZ,iT) = skCist;
            cisternaeNode(:,:,iC,iZ,iT)     = nodeIm;
            perimeterNode(:,:,iC,iZ,iT)     = perimIm;
            perimeter(:,:,iC,iZ,iT)         = bwperim(erCisternae(:,:,iC,iZ,iT)) ...
                                               & ~sk;
        end
    end
end
end
