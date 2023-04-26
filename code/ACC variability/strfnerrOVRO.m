function [sf, stdsf, errsfwalter, tlag, nsf] = strfnerrOVRO(days, S, Serr, dtlag, nlag)

% compute structure function on data defined by:
% time array: day
% data array: flux
% error array: erflux

% data is not assumed to be equi-spaced
% bin the sf into nlag bins each of width dtlag
 
% tlag(1) = 0  when there are no intervals < dtlag/2  and sf(1)=0
% tlag(1) = dtlag/2  sf(1) has average sf over all non-zero
% intervals < dtlag/2

% initially I do nothing with erflux

time = days - days(1);  % assume first measurement is at time 0

nd = length(days);     % find out number of array columns, equivalent to number of plot points in lightcurve
sf = zeros(1,nlag+1); % generate an array of zeros to store, each element represents a time-lag bin of integer multiples 
% of smallest time-lag % nlag = % number of timelag bins
nsf = sf;   % 為nsf設一全0矩陣
stdsf = sf; % 為stdsf設一全0矩陣

caliberr = 0.00;    % let calibration error = 0

flux = S./mean(S);  % calculate flux
erflux = sqrt((Serr.^2)+((caliberr.*mean(S)).^2))./mean(S); % 計算flux error

%testttt = erflux(1)

for i=1:nd         % loop through for each point on the lightcurve
    ti = time(i);               % time for this point, relative to first point
    datai = flux(i);            % normalized flux density for this point1    
    dataerri = erflux(i);
    
    dti = time(i:end) - ti;     %time difference between every other point and this point
    
    indti = round(dti./dtlag);  % bin time differences to nearest integer, as multiples of the smallest time-lag
    indtip = find(indti>=0 & indti <= nlag);  % find only positive time differences
    
    fluxj = flux(i:end);        % array to store other points that will form pair with this point
    erfluxj = erflux(i:end);
    
    
    dataj = fluxj(indtip);      % use only normalized flux densities with positive time differences - negative ones taken care of earlier
    dataerrj = erfluxj(indtip);
    
    ind = 1 + indti(indtip);      % integer time-lags, as multiples of smallest time-lag scale; add 1 to leave a zero lag in sf(1)
    diffsq = (datai-dataj).^2;
  % errwalter = (dataerri+dataerrj).*2.*(datai-dataj);
  % errwalter = sqrt(dataerri.^2 + dataerrj.^2).*2.*(datai-dataj);
  % errwalter = (dataerri - dataerrj).^2 + 2.*(abs(datai-dataerri)-abs(dataj-dataerrj)).*abs(dataerri - dataerrj);
    errwalter = (dataerri + dataerrj).^2 + 2.*abs((datai+dataerri)-(dataj-dataerrj)).*(dataerri + dataerrj);
    
    if (i == 1)
        indall = ind;
        diffsqall = diffsq;
        errwalterall = errwalter;
    else
        indall = horzcat(indall, ind);
        diffsqall = horzcat(diffsqall, diffsq);
        errwalterall = horzcat(errwalterall, errwalter);
    end
end
  
for k = 1: nlag
    findlagk = find(indall == k);

    diffsqk = diffsqall(findlagk);
    errwalterk = errwalterall(findlagk);
    
    sumdiffsq(k) = sum(diffsqk);
    sumerrwalter(k) = sum(errwalterk);

    nsf(k) = numel(diffsqk);
    stddiffsq(k) = std(diffsqk);  
end

inp = find(nsf>0);
sf(inp) = sumdiffsq(inp)./nsf(inp);   %summation divided by the total number of pairs in each bin
errsfwalter(inp) = sumerrwalter(inp)./nsf(inp);
stdsf(inp) = stddiffsq(inp);
tlag = [0:nlag].*dtlag;        % multiply back integer time-lags with smallest time-lag

return
