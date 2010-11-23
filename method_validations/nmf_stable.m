



if config.nmf_stab_plots == 1
            h = figure('Visible','off');
            imagesc(nmf_res.flat);
            title(['standard deviation of group size: ' num2str(nmf_res.std)]);
            saveas(h, [config.outpath  'nmf_expl_stab' int2str(i) '_' int2str(j) '.' config.image_format]);
            close(h);
        end

    % show stability
    disp(['standard deviation of group size: ' num2str(nmf_res.std)]);