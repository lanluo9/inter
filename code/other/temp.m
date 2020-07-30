for i = 1:nStim
%     subplot(length(Els),length(Azs),i); 
    figure
    imagesc(data_dfof_avg_all(:,:,i)); 
    colormap(gray)
    title(num2str(Stims(i,:)))
    
    set(gcf, 'Position', get(0, 'Screensize'));
    cd C:\Users\lan\Documents\repos\inter\code
    print(['map-' num2str(Stims(i,:))],'-dpdf','-fillpage')
    close

%     saveas(gcf, ['ori_tuning_fit_', num2str(icell)], 'jpg')
%             img_avg_resp(i) = mean(mean(mean(data_dfof_avg_all(:,:,i),3),2),1);
    %clim([0 max(data_dfof_avg_all(:))./2])
end