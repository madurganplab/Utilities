## polynomial of arbitrary length compatible with curve_fit

function pfit(x,p)
    pval=0;
    for i in eachindex(p)
        pval=pval.+p[i]*x.^(i-1)
    end
    return pval
end