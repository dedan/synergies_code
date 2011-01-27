function plot_rose(h, syns, pds)


if size(syns, 1) ~= 3
    warning('dedan:dim', 'this function is for 3d synergies');
    return;
end


figure(h);

for j = 1:3
    
    subplot(2,2,j);
    [~, id] = sort(pds);
    polar(pds(id), syns(j,:))
    title(['# ' num2str(j)]);
end

subplot(2,2,4)
rose(pds(1,:), 360);
title('pd distribution');
