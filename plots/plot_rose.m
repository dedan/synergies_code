function plot_rose(h, syn, pds)


if size(syns, 1) ~= 3
    warning('dedan:dim', 'this function is for 3d synergies');
    return;
end


figure(h);

for j = 1:3
    
    subplot(2,2,j);
    rose_agg = [];
    
    for k = 1:size(syns,2)
        rose_agg = [rose_agg ones(1, floor(syn(j,k) * 100)) * pds(k)]; %#ok<AGROW>
    end
    
    h_fake = rose(ones(1,100));
    hold on;
    rose(rose_agg,30);
    set(h_fake, 'Visible', 'Off');
    title(['# ' num2str(j)]);
end

subplot(2,2,4)
rose(pds(1,:), 360);
title('pd distribution');
