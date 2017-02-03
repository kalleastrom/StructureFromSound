%
if 0,
end

matches = result.matches;
rns = scores;
cmatches = result.matches;
cres = result.res;
settings = result.settings;
[dmatches,dres] = find_more_inliers_uij(matches,rns,cmatches,cres,settings);

result.res.



myplot(scores,result.matches,result.res,result.settings);
