clear
load('datelist.mat');

years = unique(Year);
len = length(years);
c = cell(len,2)

for i = 1:len
    c{i,1} = years(i);
    correct_rows = Year == years(i);
    c{i,2} = [Day(correct_rows), Week(correct_rows)];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

months = unique(Month);
len_months = length(months);
S = struct('month',{},'day',{},'weekday',{},'year',{});

for j = 1:len_months
    correct_rows_month = Month == months(j);
    S(j).month = months(j);
    S(j).day = Day(correct_rows_month);
    S(j).weekday = Weekday(correct_rows_month);
    S(j).year = Year(correct_rows_month);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%niet zelf bedacht
day_week_year = [Day, Week, Year];
dup_index = find(all(diff(day_week_year) == 0, 2));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%niet zelf bedacht
fid = fopen('duplicates.csv', 'wt');
fprintf(fid, 'Year,Month,Day,Weekday,Week\n');
for i=1:length(dup_index)
    k = dup_index(i);
    
    fprintf(fid, '%d,%s,%d,%s,%d\n', ...
        Year(k), Month(k), Day(k), Weekday(k), Week(k));
end
fclose(fid);

    