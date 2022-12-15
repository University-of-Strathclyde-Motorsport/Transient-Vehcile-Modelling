function [hs] = perp_dist(Zrcf,Zrcr,Scg2D,l)
 for i = 1:length(Zrcf)
    for j = 1:length(Zrcr)
        grad(i,j) = (Zrcr(j)-Zrcf(i))/l;
        
        hs(i,j) = abs(Scg2D(1)*grad(i,j) - Scg2D(2) + Zrcf(i))/sqrt((grad(i,j)^2 + 1));
    end
 end
end

