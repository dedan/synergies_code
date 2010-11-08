function plot_location(loc, field)

disp(field)

for i = 1:length(loc)
    subplot(length(loc),1,i)
    
    % plot relation in range
    dxy = loc(i).xy_dist(loc(i).range);
    dac = loc(i).(field)(loc(i).range);
    plot(dxy, dac , '.');
    
    % compute correlation coefficient
    [a b] = corrcoef([dxy' dac']);
    title(['corrcoef is: ' num2str(a(1,2)) ' with p: ' num2str(b(1,2))]);
    disp(['corrcoef is: ' num2str(a(1,2)) ' with p: ' num2str(b(1,2))]);
    
    % plot the means
    u = unique(dxy);
    m = NaN(size(u));
    for j = 1:length(u)
        m(j) = mean(dac(dxy == u(j)));
    end
    hold on
    plot(u,m,'*r');
    hold off
    corrcoef(u,m)
end