% Example array
data = [1 2 3 4 5];

% Calculate auto-correlation
auto_corr = corrcoef(data, data);

% The correlation coefficient between the same array will be at index (1,1) of the correlation matrix
auto_corr_coefficient = auto_corr(1, 2);

disp(['Auto-correlation coefficient: ', num2str(auto_corr_coefficient)]);