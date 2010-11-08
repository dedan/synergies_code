
%old stuff


%% residuals over handpositions

ratios = NaN(4,length(sig_sessions));
for i = 1:length(sig_sessions)
    tmp     = sig_sessions(i).res;
    h = figure('Visible','off');
    val     = [mean(tmp.r_dat); mean(tmp.r_max); mean(tmp.rr_dat); mean(tmp.rr_max)]';
    eval    = [std(tmp.r_dat);  std(tmp.r_max);  std(tmp.rr_dat);  std(tmp.rr_max)]';
    
    ratios(1,i)                     = val(1,1) / (val(1,1) - val(1,2));
    sig_sessions(i).test.ratio_pro  = ratios(1,i);
    ratios(2,i)                     = val(1,3) / (val(1,3) - val(1,4));
    sig_sessions(i).test.rratio_pro = ratios(2,i);

    ratios(3,i)                     = val(2,1) / (val(2,1) - val(2,2));
    sig_sessions(i).test.ratio_pro  = ratios(3,i);
    ratios(4,i)                     = val(2,3) / (val(2,3) - val(2,4));
    sig_sessions(i).test.rratio_pro = ratios(4,i);
    
    barweb( val,eval,.5,[], [],'Hand Position','Explained variance', 'default', [], {'dat','max','Rdat', 'Rmax'});
    saveas(h, [config.outpath 'sep' filesep  'sep_' sig_sessions(i).info.name  '.' config.image_format]);
    close(h);
end

h = figure('Visible','off');
subplot 211
plot(ratios(1:2,:)');
legend('normal','random');
title('pronation');
subplot 212
plot(ratios(3:4,:)');
legend('normal','random');
title('supination');
saveas(h, [config.outpath 'sep_all_ratio.' config.image_format]);
close(h);

clear tmp val eval ans ratios



