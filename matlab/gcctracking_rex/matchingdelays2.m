function [uref,urefk] = matchingdelays2(u,settings)
%MATCHINGDELAYS - Limits the delay candidates
%
%    This function...
%
%    uref = MATCHINGDELAYS(u,settings)
%
%    Input:
%    u  - delay data output from getdelays
%    settings - struct that must have...
%
%    Output:
%    uref - delays after cleanup

channels = settings.channels;
[m,n] = size(u{channels(1),channels(1)});
uref = cell(settings.mm,settings.mm);
urefk = cell(settings.mm,settings.mm);
for refChannel = settings.channels;
    for ch = channels(channels~=refChannel)
        [refChannel ch]
        %keyboard;
        candidates = NaN(1+m*m*(numel(channels)-2),n);
        candidates(1:m,:) = u{refChannel,ch};
        %Using u{refChannel,ch} = urefChannel{refChannel,k}-u{ch,k}:
        pp = 1; %loop counter
        for k = setdiff(channels,[refChannel,ch])
            candidates(2+(pp-1)*m*m:pp*m*m+1,:) = dfcombs(u{refChannel,k},u{ch,k});
            pp = pp+1;
        end
        %Finding the candidates that have at least minNbrOfInliers inliers:
        [bin,binNbr] = histc(candidates,-settings.sw:settings.binSize:settings.sw);
        ind = bin > settings.minNbrOfInliers;
        urch = NaN(max(sum(ind)),size(bin,2));
        urchk = NaN(max(sum(ind)),size(bin,2));
        for cc = find(sum(ind) > 0)
            pp = 1; %loop counter
            for k = find(ind(:,cc))'
%                 if ((refChannel == 7) & (ch == 8) & (cc == 1862)),
%                     keyboard;
%                 end
                tmp = u{refChannel,ch}(:,cc);
                tmp = tmp(find(isfinite(tmp)));
                tmp2 = intersect(tmp,candidates(binNbr(:,cc) == k,cc));
                if length(tmp2)==1,
                    urch(pp,cc)=tmp2;
                    urchk(pp,cc)=bin(k,cc);
                else
                    urch(pp,cc) = median(candidates(binNbr(:,cc) == k,cc));
                    urchk(pp,cc)=-bin(k,cc);
                end
                pp = pp+1;
            end
        end
        uref{refChannel,ch} = urch;
        urefk{refChannel,ch} = urchk;
    end
end
end

%Help function:
function combs = dfcombs(A,B)
%...
%  A, B, same size

[m,n] = size(A);
combs = NaN(m*m,n);
for k = 1:m
    combs(1+(k-1)*m:k*m,:) = A-B;
    B = circshift(B,-1);
end
end