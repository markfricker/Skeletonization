function bwOut = erSkeletonTruncate(bwIn, erMask, erCisternae, cellBoundary, p)
%ERSKELETONTRUNCATE  Filter ER skeleton by mask, size and cisternae connectivity.
%
%   bwOut = erSkeletonTruncate(bwIn, erMask, erCisternae, cellBoundary, p)
%
% Post-processes a raw ER skeleton by:
%   1. Optionally masking with the ER mask and cell boundary.
%   2. Removing fragments smaller than p.minAreaPixels.
%   3. Optionally retaining only the connected component(s) that overlap
%      the cisternae (greatest-connected-component, GCC).
%
% INPUTS
%   bwIn         – [nY nX nC nZ nT] logical skeleton from erSkeletonBuild.
%   erMask       – [nY nX ...] logical ER mask, or [] to skip.
%                  May have fewer C/Z/T dimensions (broadcast via min indexing).
%   erCisternae  – [nY nX ...] logical cisternae mask, or [] to skip GCC.
%   cellBoundary – [nY nX ...] logical cell-interior mask (1 = inside cell),
%                  or [] to skip.
%   p            – parameter struct:
%                    p.maskUse       – logical; apply erMask and cellBoundary
%                    p.gccUse        – logical; retain only component(s)
%                                      connected to cisternae
%                    p.minAreaPixels – integer >= 1; minimum skeleton fragment
%                                      length in pixels (use
%                                      fwhmTarget_px * 2 as a guide)
%
% OUTPUT
%   bwOut – [nY nX nC nZ nT] filtered logical skeleton, same size as bwIn.

[nY, nX, nC, nZ, nT] = size(bwIn);

hasMask       = ~isempty(erMask);
hasCisternae  = ~isempty(erCisternae);
hasBoundary   = ~isempty(cellBoundary);

if hasMask,      [~,~,mC,mZ,mT] = size(erMask);      end
if hasCisternae, [~,~,fC,fZ,fT] = size(erCisternae); end
if hasBoundary,  [~,~,bC,bZ,bT] = size(cellBoundary); end

minArea = max(1, round(p.minAreaPixels));

bwOut = bwIn;

for iT = 1:nT
    for iZ = 1:nZ
        for iC = 1:nC
            sk = bwIn(:,:,iC,iZ,iT);

            % --- apply ER mask and cell boundary -------------------------
            if p.maskUse
                if hasMask
                    sk = sk & erMask(:,:, min(iC,mC), min(iZ,mZ), min(iT,mT));
                end
                if hasBoundary
                    sk = sk & cellBoundary(:,:, min(iC,bC), min(iZ,bZ), min(iT,bT));
                end
            end

            % also always clip to cell boundary regardless of maskUse
            if hasBoundary
                sk(~cellBoundary(:,:, min(iC,bC), min(iZ,bZ), min(iT,bT))) = false;
            end

            % --- remove small fragments ----------------------------------
            sk = bwareafilt(sk, [minArea, Inf]);

            % --- keep only component(s) connected to cisternae -----------
            if p.gccUse && hasCisternae
                cist = erCisternae(:,:, min(iC,fC), min(iZ,fZ), min(iT,fT));
                if any(cist(:))
                    tmp = sk;
                    tmp(cist) = true;
                    [r, c] = find(cist);
                    tmp = bwselect(tmp, c, r);
                    sk = tmp & sk;
                end
            end

            bwOut(:,:,iC,iZ,iT) = sk;
        end
    end
end
end
