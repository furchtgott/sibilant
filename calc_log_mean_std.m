function [loggenemeans, loggenestds] = calc_log_mean_std(data, iunique)

loggenemeans = zeros(size(data,1), max(iunique));
loggenestds = zeros(size(data,1), max(iunique));

for i=1:max(iunique)
    loggenemeans(:,i) = mean(data(:,ismember(iunique,i)),2);
    loggenestds(:,i) = std(data(:,ismember(iunique,i)),[],2);
end

end
