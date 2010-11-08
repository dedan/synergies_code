
%% check distributions of dot products between positiv vectores

x = normr(randn(1000,2));
y = normr(randn(1000,2));

for i = 1:1000
    bla(i) = abs(x(i,:)*y(i,:)');
end

hist(bla,100)


% das hier stand in nonevoked_syn_analysis drin im consistency over sessions teil 
% und ich habe keine ahnung mehr was es bedeuten soll
% wollte es wohl irgendwie zur quantifierzung nehmen, das versuche ich jetz dafür
% mal mit der k-means cluster verteilung

shuf = struct;
for j = 1:100
    for k = 1:size(fin_syns.(all_names{i}),2)
        ran = randperm(size(all(i).flat,1));
        shuf(j).syns(:,k) = all(i).flat(ran(1:3),k);
    end
end

n = 1;
for j = 1:length(shuf)
    for k = j+1:length(shuf)
        scores(n,:) = matchNscore(shuf(j).syns', shuf(k).syns'); %#ok<AGROW>
        n = n+1;
    end
end



n = 1;
for j = 1:length(all_syns)
    for k = j+1:length(all_syns)
        sc2(n,:)      = matchNscore(all_syns(j).(all_names{i})', all_syns(k).(all_names{i})'); %#ok<AGROW>
        n = n+1;
    end
end

