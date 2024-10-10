function [tvals, pvals] = jvh_unpaired_ttest_pooled_sd(Amean_vec,Bmean_vec,Asd_vec,Bsd_vec,NA,NB)
% iterate over vectors of means and sd of group A and group B to compute the t-value and p-value for each element.
% Sample size NA of group A and NB of group B needs to be given.
% This is for 2 unpaired groups, using a two-tailed t-test, equal variance between groups
arguments
    Amean_vec {mustBeVector}
    Bmean_vec {mustBeVector}
    Asd_vec   {mustBeVector}
    Bsd_vec   {mustBeVector}
    NA        {mustBePositive, mustBeInteger}
    NB        {mustBePositive, mustBeInteger}
end

assert(length(Amean_vec)==length(Asd_vec))
assert(length(Bmean_vec)==length(Bsd_vec))
assert(length(Amean_vec)==length(Bmean_vec))

pvals = nan(length(Amean_vec),1);
tvals = nan(length(Amean_vec),1);

for i = 1:length(Amean_vec)
    Amean =Amean_vec(i); 
    Asd= Asd_vec(i); 
    Bmean= Bmean_vec(i); 
    Bsd= Bsd_vec(i);

    v = NA + NB -2; %degrees of freedom
    Sp_square = ((NA-1)*Asd^2 + (NB-1)*Bsd^2)/v; %pooled sd
    
    tval = (Amean-Bmean) / (sqrt(Sp_square)*sqrt(1/NA + 1/NB));       % Calculate T-Statistic
    if isnan(tval)
        tval = 0;
    end
    tdist2T = @(t,v) (1-betainc(v/(v+t^2),v/2,0.5));    % 2-tailed t-distribution
    pvals(i) = 1-[tdist2T(tval,v)];
    tvals(i) = tval;
end

