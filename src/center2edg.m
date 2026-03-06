function binedges = center2edg(bincenters)

% defaults
ndigit = 10;

binsz = bincenters(2)-bincenters(1);
logbinsz = log10(bincenters(2))-log10(bincenters(1));

if all(round(diff(log10(bincenters)),ndigit)==round(logbinsz,ndigit)) % log-spaced
    logbincenters = log10(bincenters);
    logbinedges = mean([logbincenters(1:end-1);logbincenters(2:end)],1);
    logbinedges = ...
        [logbinedges(1)-logbinsz,logbinedges,logbinedges(end)+logbinsz];
    binedges = 10.^logbinedges;
    
elseif all(round(diff(bincenters),ndigit)==round(binsz,ndigit)) % lineraly spaced
    binedges = mean([bincenters(1:end-1);bincenters(2:end)],1);
    binedges = [binedges(1)-binsz,binedges,binedges(end)+binsz];

else
    nbins = length(bincenters);
    logbincenters = log10(bincenters);
    logbinedges = zeros(1,nbins+1);
    logbinedges(1) = logbincenters(1) - ...
        (logbincenters(2)-logbincenters(1))/2;
    logbinedges(end) = logbincenters(nbins) + ...
        (logbincenters(nbins)-logbincenters(nbins-1))/2;
    for bin = 2:nbins
        logbinedges(bin) = logbincenters(bin-1) + ...
            (logbincenters(bin)-logbincenters(bin-1))/2;
    end
    binedges = 10.^logbinedges;
end