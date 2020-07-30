for i = 1:nStim
%     subplot(length(Els),length(Azs),i); 
    subplot(4,3,i); 
%     figure
    imagesc(data_dfof_avg_all(:,:,i)); 
    colormap(gray)
    title(num2str(Stims(i,:)))
    
end

set(gcf, 'Position', get(0, 'Screensize'));
cd C:\Users\lan\Documents\repos\inter\code
print(['map'],'-dpdf','-fillpage')
% close
