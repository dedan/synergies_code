function plot_rose(f, syns, pds, types)

colors = ['r' 'g' 'b'];

for j = 1:size(syns, 1)

    subplot(3,2,j);
    h_fake = rose(ones(1,100));
    hold on;

    for t = 1:3 	% different muscle types
        rose_agg = [];
	    for k = find(types == t)
	        rose_agg = [rose_agg ones(1, floor(syns(j,k) * 100)) * pds(k)]; %#ok<AGROW>
	    end
        h = rose(rose_agg, 30);
        x = get(h, 'XData');
        y = get(h, 'YData');
        p = patch(x, y, colors(t));
	end
    set(h_fake, 'Visible', 'Off');
    title(['# ' num2str(j)]);
end

subplot(3,2,5)
for t = 1:3
    k = find(types == t);
    h = rose(pds(1, k), 180);
    hold on
    x = get(h, 'XData');
    y = get(h, 'YData');
    p = patch(x, y, colors(t));
end

title('pd distribution');
