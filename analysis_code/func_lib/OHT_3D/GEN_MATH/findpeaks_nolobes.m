function [pks_locs,pks_vals] = findpeaks_nolobes(sig,x,relpower_req,neighbormult_req)

%findpeaks_nolobes: This function performs similar analysis to findpeaks,
%but uses two criteria to select only those peaks that are central to a
%signal with sidelobes. The function uses two strategies
%to eliminate spurious sidelobes that can occur in noisy data:
%   1) Only peaks with amplitude (height above minimum value) greater than
%   some relative power value are selected
%   2) Only peaks that are a multiplicative factor larger than their
%   neighboring peaks are selected.
%
%[pks_locs,pks_vals] =
%findpeaks_nolobes(sig,x,relpower_req,neighbormult_req)
%where:
%   -pks_locs & pks_vals are the location and values of the selected peaks,
%   respectively (vector of all peaks that meet the given criteria)
%   -sig is the input signal (the signal on which the peaks are being found
%   -x is the x axis values for the signal
%   -relpower_req is a scalar <1 representing the relative power (relative
%   to the overall maximum amplitude of the signal) that a peak must have
%   to be selected
%   -neighbormult_req is a scalar >1 representing the relative ratio
%   (relative to both neighbors) that a peak must have in order to be
%   selected.

min_sig = min(sig);
max_amp = max(sig) - min_sig;

peak_thresh = (min_sig + relpower_req*max_amp);
[pks_all, locs_all] = findpeaks(sig,x,'MinPeakHeight',peak_thresh);

num_findpks = numel(pks_all);
%Calculate relative size, relative to neighbor peaks
pks_rel = zeros(num_findpks,1);
%Calculate the average peak height ratio, relative to its neighbors
pks_rel(1) = (pks_all(1)-min_sig)/(pks_all(2)-min_sig);
for k = 2:1:(num_findpks-1)
    pks_rel(k) = min([((pks_all(k)-min_sig)/(pks_all(k+1)-min_sig)), ...
        ((pks_all(k)-min_sig)/(pks_all(k-1)-min_sig))]);
end
pks_rel(num_findpks) = (pks_all(num_findpks)-min_sig)/(pks_all(num_findpks-1)-min_sig);
%Filter to peaks that are at least neighbormult_req * their
%neighbor in height
pks_filt = pks_rel > neighbormult_req;
pks_locs = locs_all(pks_filt);
pks_vals = pks_all(pks_filt);
%Sort the peaks that meet all requiremnts in descending order
[~,I] = sort(pks_vals,'descend');
pks_vals = pks_vals(I);
pks_locs = pks_locs(I);