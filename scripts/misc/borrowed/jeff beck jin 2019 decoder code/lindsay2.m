
%%%%%%%%%%%%THE ABOVE SHOULD BE CROSSVALIDATED TO MAKE COMPARISON FAIR>>>>>
clear all
close all
% load 170510_i574_runs-002-003_testResp_adapt.mat
% load 170510_i574_runs-002-003_testResp_forJeff.mat
% 
% load 170323_i689_runs-002-003_testResp_adapt.mat
% load 170323_i696_runs-002-003_testResp_forJeff.mat  
% 
% load 170327_i684_runs-002-003_testResp_adapt.mat  
% load 170327_i684_runs-002-003_testResp_forJeff.mat  
% 
% load 170503_i712_runs-002-003_testResp_adapt.mat 



% load 170323_i689_runs-002-003_testResp_adapt
% load 170323_i696_runs-002-003_testResp_adapt  
% load 170324_i674_runs-002-003_testResp_adapt  
% load 170327_i684_runs-002-003_testResp_adapt  
 load 170503_i711_runs-002-003_testResp_adapt  
% load 170503_i712_runs-002-003_testResp_adapt  
% load 170510_i574_runs-002-003_testResp_adapt

% Preprocess
% note that 3 is control case


% prepare data matrices
% (1) normalize based upon z-scores from the control condition 
% (2) make Y labels 1 to 8 with 8 being vertical
for j=1:3
    data{j}.X=[];
    data{j}.Y=[];

    for n=1:8   % uses one neurons in ind
        data{j}.X = [data{j}.X;ppResp{j,n}';];
        data{j}.Y = [data{j}.Y;n*ones(size(ppResp{j,n},2),1);];
    end
    idx=~any(isnan(data{j}.X),2);
    data{j}.X=data{j}.X(idx,:);
    data{j}.Y=data{j}.Y(idx,1);
    NS(j)=size(data{j}.Y,1);
    data{j}.Xraw=[data{j}.X,ones(size(data{j}.X,1),1)];
end

mu=mean([data{1}.X;data{2}.X;data{3}.X]);
stdX=std([data{1}.X;data{2}.X;data{3}.X]);

for j=1:3
    data{j}.X=bsxfun(@plus,data{j}.X,-mu);
    data{j}.X=bsxfun(@times,data{j}.X,1./stdX);
    data{j}.X=[data{j}.X,ones(size(data{j}.X,1),1)];
end
issparse=10;
maxiters=2000;
NR = size(data{3}.X,2);
m3=VonMesisRegression(NR,issparse);
m3.fit(data{3}.Y/8*2*pi,data{3}.X,2000);
ind = find(abs(m3.a.mean)./sqrt(diag(m3.a.var))>2 | abs(m3.b.mean)./sqrt(diag(m3.b.var))>2);
NR=Inf;
NR=min(length(ind),NR);
issparse=0;

dataAll.X=[];
dataAll.Y=[];
dataAll.cond=[];
for j=1:3
    data{j}.X=data{j}.X(:,ind(1:NR));
    dataAll.X=[dataAll.X;data{j}.X];
    dataAll.Y=[dataAll.Y;data{j}.Y];
    dataAll.cond=[dataAll.cond;j*ones(size(data{j}.Y));];
end

%do cross validations for self-model comparisons
issparse=NR/2;
maxiters=500;
ntheta=500;
theta = [1:ntheta]/ntheta*2*pi-pi;
thetaS = [1:8]/8*2*pi-pi;
dtheta=theta(2)-theta(1);
for j=1:3
    for n=1:size(data{j}.X)
        m1 = VonMesisRegression(NR,issparse);
        idx = [1:n-1,n+1:size(data{j}.X,1)];
        m1.fit(data{j}.Y(idx,1)*2*pi/8,data{j}.X(idx,:),maxiters);
        [postmu{j,j}(n,1),postkappa{j,j}(n,1)]=m1.getPredictions(data{j}.X(n,:));
        
        pr{j,j}(n,:)=m1.getPdf(theta,data{j}.X(n,:));
        DV{j,j}(n,1) = 1-sum(pr{j,j}(n,abs(theta)<pi/8),2)/sum(pr{j,j}(n,:),2);
        DV8{j,j}(n,1) = 1-m1.getPdf(0,data{j}.X(n,:))/sum(m1.getPdf(thetaS,data{j}.X(n,:)));
        thetaDV{j,j}(n,1) = abs(postmu{j,j}(n,1));
%        B1 = mnrfit(data{j}.X(idx,:),(data{j}.Y(idx,1)~=8)+1);
%        temp=mnrval(B1,data{j}.X(n,:));
%        logitDV{j,j}(n,1) = temp(2);
%        temp=mnrval(B1,data{j}.X(n,:));
%        logitDV{j,j}(n,1) = temp(2);
        B1 = glmfit(data{j}.X(idx,:),(data{j}.Y(idx,1)~=8),'binomial','link','logit','Constant','off');
        logitDV{j,j}(n,1) = glmval(B1,data{j}.X(n,:),'logit','Constant','off');
    end
    model{j} = VonMesisRegression(NR,issparse);
    model{j}.fit(data{j}.Y*2*pi/8,data{j}.X,maxiters)
    inSample.DV{j} = 1-sum(model{j}.getPdf(theta,data{j}.X))*dtheta;
    [inSample.postmu{j},inSample.postkappa{j}]=model{j}.getPredictions(data{j}.X);
    sumDV{j} = mean(data{j}.Xraw(:,ind(1:NR)),2);
    [B{j},dev,stats] = glmfit(data{j}.X,(data{j}.Y~=8),'binomial','link','logit','Constant','off');
    4-j
end

[Ball,dev,stats] = glmfit(dataAll.X,(dataAll.Y~=8),'binomial','link','logit','Constant','off');
mAll = VonMesisRegression(NR,issparse);
mAll.fit(dataAll.Y*2*pi/8,dataAll.X,maxiters)


for j=1:3
for k=1:3
    if(j~=k)
        pr{k,j}=model{k}.getPdf(theta,data{j}.X);
        DV{k,j} = 1-sum(pr{k,j}(:,abs(theta)<2*pi/16),2)./sum(pr{k,j},2);        
        DV8{k,j} = 1-model{k}.getPdf(0,data{j}.X)./sum(model{k}.getPdf(thetaS,data{j}.X),2);
        [postmu{k,j},postkappa{k,j}]=model{k}.getPredictions(data{j}.X);
        thetaDV{k,j} = abs(postmu{k,j})/pi;
        logitDV{k,j} = glmval(B{k},data{j}.X,'logit','Constant','off');
    end
end
end

for j=1:3
    prAll{j}=mAll.getPdf(theta,data{j}.X);
    DVAll{j} = 1-sum(prAll{j}(:,abs(theta)<2*pi/16),2)./sum(prAll{j},2);        
    DV8All{j} = 1-mAll.getPdf(0,data{j}.X)./sum(mAll.getPdf(thetaS,data{j}.X),2);
    [postmuAll{j},postkappaAll{j}]=mAll.getPredictions(data{j}.X);
    thetaDVAll{j} = abs(postmuAll{j})/pi;
    logitDVAll{j} = glmval(Ball,data{j}.X,'logit','Constant','off');
    sumDVall{j}=sumDV{j};
end

% l=0;
% for k=1:3
% for j=1:3
%     l=l+1;
%     figure(l)
%     subplot(5,1,1), hist(DV{k,j}(data{j}.Y==8))
%     ax=axis; ax(1)=0;ax(2)=1;axis(ax);     
%     title(['OptDV Train = ',num2str(k),' Test = ',num2str(j)])
%     ylabel('\theta=0')
%     subplot(5,1,2), hist(DV{k,j}(data{j}.Y==1 | data{j}.Y==7))
%     ax=axis; ax(1)=0;ax(2)=1;axis(ax);    
%     ylabel('\theta=+- 22.5')
%     subplot(5,1,3), hist(DV{k,j}(data{j}.Y==2 | data{j}.Y==6))
%     ax=axis; ax(1)=0;ax(2)=1;axis(ax);    
%     ylabel('\theta=+- 45')
%     subplot(5,1,4), hist(DV{k,j}(data{j}.Y==3 | data{j}.Y==5))
%     ax=axis; ax(1)=0;ax(2)=1;axis(ax);    
%     ylabel('\theta=+- 67.5')
%     subplot(5,1,5), hist(DV{k,j}(data{j}.Y==4))
%     ax=axis; ax(1)=0;ax(2)=1;axis(ax);    
%     ylabel('\theta=90')
% end
% end
    
NDC=100;
l=0;
figure
not8 = [1,2,6,7];
% not8 = [1,7];
for k=1:3
for j=1:3
    l=l+1;
    subplot(3,3,l) 
    for n=1:NDC
        dc(n)=(n-1)/(NDC-1);
        fpDV{k,j}(n) = mean(DV{k,j}(data{j}.Y==8)>dc(n));
%        cdDV{k,j}(n) = mean(DV{k,j}(data{j}.Y~=8)>dc(n));
%        cdDV{k,j}(n) = mean(DV{k,j}(data{j}.Y==7 | data{j}.Y==1)>dc(n));
        cdDV{k,j}(n) = mean(DV{k,j}(logical(sum(data{j}.Y==not8,2))) > dc(n));
        ndist=sum(data{j}.Y==8);
        ntarget=sum(sum(data{j}.Y==not8));
        perfDV{k,j}(n) = ((1-fpDV{k,j}(n))*ndist + cdDV{k,j}(n)*ntarget)/(ndist+ntarget);        
    end
    plot(fpDV{k,j},cdDV{k,j},'-o')
    hold on
    
    for n=1:NDC
        fpDV8{k,j}(n) = mean(DV8{k,j}(data{j}.Y==8)>dc(n));
%        cdDV8{k,j}(n) = mean(DV8{k,j}(data{j}.Y~=8)>dc(n));
        cdDV8{k,j}(n) = mean(DV8{k,j}(logical(sum(data{j}.Y==not8,2)))>dc(n));
        ndist=sum(data{j}.Y==8);
        ntarget=sum(sum(data{j}.Y==not8));
        perfDV8{k,j}(n) = ((1-fpDV8{k,j}(n))*ndist + cdDV8{k,j}(n)*ntarget)/(ndist+ntarget);        
    end
    plot(fpDV8{k,j},cdDV8{k,j},'-x')
    
    for n=1:NDC
        fpthetaDV{k,j}(n) = mean(thetaDV{k,j}(data{j}.Y==8)>dc(n));
%        cdthetaDV{k,j}(n) = mean(thetaDV{k,j}(data{j}.Y~=8)>dc(n)*2*pi/4);
        cdthetaDV{k,j}(n) = mean(thetaDV{k,j}(logical(sum(data{j}.Y==not8,2)))>dc(n));
        ndist=sum(data{j}.Y==8);
        ntarget=sum(sum(data{j}.Y==not8));
        perfthetaDV{k,j}(n) = ((1-fpthetaDV{k,j}(n))*ndist + cdthetaDV{k,j}(n)*ntarget)/(ndist+ntarget);        
    end
    plot(fpthetaDV{k,j},cdthetaDV{k,j},'-+')

    for n=1:NDC
        fpsumDV{k,j}(n) = mean(sumDV{j}(data{j}.Y==8)>dc(n)*(max(sumDV{j}-min(sumDV{j})))+min(sumDV{j}));
%        cdsumDV{k,j}(n) = mean(sumDV{j}(data{j}.Y~=8)>dc(n)*(max(sumDV{j}-min(sumDV{j})))+min(sumDV{j}));
        cdsumDV{k,j}(n) = mean(sumDV{j}(logical(sum(data{j}.Y==not8,2)))>dc(n)*(max(sumDV{j}-min(sumDV{j})))+min(sumDV{j}));
        ndist=sum(data{j}.Y==8);
        ntarget=sum(sum(data{j}.Y==not8));
        perfsumDV{k,j}(n) = ((1-fpsumDV{k,j}(n))*ndist + cdsumDV{k,j}(n)*ntarget)/(ndist+ntarget);        
    end
    plot(fpsumDV{k,j},cdsumDV{k,j},'-s')
    
    for n=1:NDC
        fplogit{k,j}(n) = mean(logitDV{k,j}(data{j}.Y==8)>dc(n));
%        cdlogit{k,j}(n) = mean(logitDV{k,j}(data{j}.Y~=8)>dc(n));
        cdlogit{k,j}(n) = mean(logitDV{k,j}(logical(sum(data{j}.Y==not8,2)))>dc(n));        
        ndist=sum(data{j}.Y==8);
        ntarget=sum(sum(data{j}.Y==not8));
        perflogit{k,j}(n) = ((1-fplogit{k,j}(n))*ndist + cdlogit{k,j}(n)*ntarget)/(ndist+ntarget);        
    end
    plot(fplogit{k,j},cdlogit{k,j},'-*')
    plot(0:1,0:1,'k-')
    xlabel('FP')
    ylabel('CD')
    title(['Train = ',num2str(k),' Test = ',num2str(j)])
    legend('Optimal','8 Category opt','Theta est','Sum','logit','refline')
    hold off

end
end



NDC=100;
l=0;
figure
for j=1:3
    l=l+1;
    subplot(3,1,l) 
    for n=1:NDC
        dc(n)=(n-1)/(NDC-1);
        fpDVAll{j}(n) = mean(DVAll{j}(data{j}.Y==8)>dc(n));
%        cdDVAllj}(n) = mean(DVAll{j}(data{j}.Y~=8)>dc(n));
%        cdDVAll{j}(n) = mean(DVAll{j}(data{j}.Y==7 | data{j}.Y==1)>dc(n));
        cdDVAll{j}(n) = mean(DVAll{j}(logical(sum(data{j}.Y==not8,2))) > dc(n));
        ndist=sum(data{j}.Y==8);
        ntarget=sum(sum(data{j}.Y==not8));
        perfDVAll{j}(n) = ((1-fpDVAll{j}(n))*ndist + cdDVAll{j}(n)*ntarget)/(ndist+ntarget);        
    end
    plot(fpDVAll{j},cdDVAll{j},'-o')
    hold on
    
    for n=1:NDC
        fpDV8All{j}(n) = mean(DV8All{j}(data{j}.Y==8)>dc(n));
%        cdDV8All{j}(n) = mean(DV8All{j}(data{j}.Y~=8)>dc(n));
        cdDV8All{j}(n) = mean(DV8All{j}(logical(sum(data{j}.Y==not8,2)))>dc(n));
        ndist=sum(data{j}.Y==8);
        ntarget=sum(sum(data{j}.Y==not8));
        perfDV8All{j}(n) = ((1-fpDV8All{j}(n))*ndist + cdDV8All{j}(n)*ntarget)/(ndist+ntarget);        
    end
    plot(fpDV8All{j},cdDV8All{j},'-x')
    
    for n=1:NDC
        fpthetaDVAll{j}(n) = mean(thetaDVAll{j}(data{j}.Y==8)>dc(n));
%        cdthetaDVAll{j}(n) = mean(thetaDVAll{j}(data{j}.Y~=8)>dc(n));
        cdthetaDVAll{j}(n) = mean(thetaDVAll{j}(logical(sum(data{j}.Y==not8,2)))>dc(n));
        ndist=sum(data{j}.Y==8);
        ntarget=sum(sum(data{j}.Y==not8));
        perfthetaDVAll{j}(n) = ((1-fpthetaDVAll{j}(n))*ndist + cdthetaDVAll{j}(n)*ntarget)/(ndist+ntarget);        
    end
    plot(fpthetaDVAll{j},cdthetaDVAll{j},'-+')

    for n=1:NDC
        fpsumDVAll{j}(n) = mean(sumDV{j}(data{j}.Y==8)>dc(n)*(max(sumDV{j}-min(sumDV{j})))+min(sumDV{j}));
%        cdsumDVAll{j}(n) = mean(sumDV{j}(data{j}.Y~=8)>dc(n)*(max(sumDV{j}-min(sumDV{j})))+min(sumDV{j}));
        cdsumDVAll{j}(n) = mean(sumDV{j}(logical(sum(data{j}.Y==not8,2)))>dc(n)*(max(sumDV{j}-min(sumDV{j})))+min(sumDV{j}));
        ndist=sum(data{j}.Y==8);
        ntarget=sum(sum(data{j}.Y==not8));
        perfsumDVAll{j}(n) = ((1-fpsumDVAll{j}(n))*ndist + cdsumDVAll{j}(n)*ntarget)/(ndist+ntarget);        
    end
    plot(fpsumDVAll{j},cdsumDVAll{j},'-s')
    
    for n=1:NDC
        fplogitAll{j}(n) = mean(logitDVAll{j}(data{j}.Y==8)>dc(n));
%        cdlogitAll{j}(n) = mean(logitDVAll{j}(data{j}.Y~=8)>dc(n));
        cdlogitAll{j}(n) = mean(logitDVAll{j}(logical(sum(data{j}.Y==not8,2)))>dc(n));        
        ndist=sum(data{j}.Y==8);
        ntarget=sum(sum(data{j}.Y==not8));        
        perflogitAll{j}(n) = ((1-fplogitAll{j}(n))*ndist + cdlogitAll{j}(n)*ntarget)/(ndist+ntarget);
    end
    plot(fplogitAll{j},cdlogitAll{j},'-*')
    plot(0:1,0:1,'k-')
    xlabel('FP')
    ylabel('CD')
    title(['Train = All',' Test = ',num2str(j)])
    legend('Optimal','8 Category opt','Theta est','Sum','logit','refline')
    hold off

end




hand=figure;
l=hand.Number-1;
for k=1:3
    figure(l+2*k-1)
    clf
    figure(l+2*k)
    clf
for j=1:3
    idx = find(fpDV{k,j} ~= cdDV{k,j});
    minidx=min(idx)-1;
    maxidx=max(idx)+1;
    idx=[max(minidx,1):min(maxidx,length(dc))];
    figure(l+2*k-1)
    subplot(3,2,1), scatter(dc(idx),fpDV{k,j}(idx)), hold on
    title(['Optimal with train = ',num2str(k)])
    %legend('FP 250','FP 500','FP 750')
    figure(l+2*k)
    subplot(3,2,1), scatter(dc(idx),cdDV{k,j}(idx)), hold on
    title(['Optimal with train = ',num2str(k)])
    %legend('CD 250','CD 500','CD 750')

    idx = find(fpDV8{k,j} ~= cdDV8{k,j});
    minidx=min(idx)-1;
    maxidx=max(idx)+1;
    idx=[max(minidx,1):min(maxidx,length(dc))];
    figure(l+2*k-1)
    subplot(3,2,2), scatter(dc(idx),fpDV8{k,j}(idx)), hold on
    %legend('FP 250','FP 500','FP 750')
    title('Optimal8')
    figure(l+2*k)
    subplot(3,2,2), scatter(dc(idx),cdDV8{k,j}(idx)), hold on
    title('Optimal8')
    %legend('CD 250','CD 500','CD 750')

    idx = find(fpthetaDV{k,j} ~= cdthetaDV{k,j});
    minidx=min(idx)-1;
    maxidx=max(idx)+1;
    idx=[max(minidx,1):min(maxidx,length(dc))];
    figure(l+2*k-1)
    subplot(3,2,3), scatter(dc(idx)*360/2/2,fpthetaDV{k,j}(idx)), hold on
    title('theta estimated')
    %legend('FP 250','FP 500','FP 750')
    figure(l+2*k)
    subplot(3,2,3), scatter(dc(idx)*360/2/2,cdthetaDV{k,j}(idx)), hold on
    title('theta estimated')
    %legend('CD 250','CD 500','CD 750')

    idx = find(fpsumDV{k,j} ~= cdsumDV{k,j});
    minidx=min(idx)-1;
    maxidx=max(idx)+1;
    idx=[max(minidx,1):min(maxidx,length(dc))];
    figure(l+2*k-1)    
    subplot(3,2,4), scatter(dc(idx)*(max(sumDV{j}-min(sumDV{j})))+min(sumDV{j}),fpsumDV{k,j}(idx)), hold on
    title('sumDV')
    %legend('FP 250','FP 500','FP 750')
    figure(l+2*k)
    subplot(3,2,4), scatter(dc(idx)*(max(sumDV{j}-min(sumDV{j})))+min(sumDV{j}),cdsumDV{k,j}(idx)), hold on
    title('sumDV')
    %legend('CD 250','CD 500','CD 750')    
    
    idx = find(fplogit{k,j} ~= cdlogit{k,j});
    minidx=min(idx)-1;
    maxidx=max(idx)+1;
    idx=[max(minidx,1):min(maxidx,length(dc))];
    figure(l+2*k-1)    
    subplot(3,2,5), scatter(dc(idx),fplogit{k,j}(idx)), hold on
    if(j==3)
        legend('FP 250','FP 750','FP Control')    
    end
    title('weighted sum')
    xlabel('decision criterion')
    figure(l+2*k)
    subplot(3,2,5), scatter(dc(idx),cdlogit{k,j}(idx)), hold on
    title('weighted sum')
    xlabel('decision criterion')
    if(j==3)
        legend('CD 250','CD 750','CD Control')    
    end
end
end





hand=figure;
l=hand.Number-1;
k=1;
figure(l+2*k-1)
clf
figure(l+2*k)
clf
for j=1:3
    idx = find(fpDVAll{j} ~= cdDVAll{j});
    minidx=min(idx)-1;
    maxidx=max(idx)+1;
    idx=[max(minidx,1):min(maxidx,length(dc))];
    figure(l+2*k-1)
    subplot(3,2,1), scatter(dc(idx),fpDVAll{j}(idx)), hold on
    title(['Optimal with train = All'])
    %legend('FP 250','FP 500','FP 750')
    figure(l+2*k)
    subplot(3,2,1), scatter(dc(idx),cdDVAll{j}(idx)), hold on
    title(['Optimal with train = All'])
    %legend('CD 250','CD 500','CD 750')

    idx = find(fpDV8All{j} ~= cdDV8All{j});
    minidx=min(idx)-1;
    maxidx=max(idx)+1;
    idx=[max(minidx,1):min(maxidx,length(dc))];
    figure(l+2*k-1)
    subplot(3,2,2), scatter(dc(idx),fpDV8All{j}(idx)), hold on
    %legend('FP 250','FP 500','FP 750')
    title('Optimal8')
    figure(l+2*k)
    subplot(3,2,2), scatter(dc(idx),cdDV8All{j}(idx)), hold on
    title('Optimal8')
    %legend('CD 250','CD 500','CD 750')

    idx = find(fpthetaDVAll{j} ~= cdthetaDVAll{j});
    minidx=min(idx)-1;
    maxidx=max(idx)+1;
    idx=[max(minidx,1):min(maxidx,length(dc))];
    figure(l+2*k-1)
    subplot(3,2,3), scatter(dc(idx)*360/2/2,fpthetaDVAll{j}(idx)), hold on
    title('theta estimated')
    %legend('FP 250','FP 500','FP 750')
    figure(l+2*k)
    subplot(3,2,3), scatter(dc(idx)*360/2/2,cdthetaDVAll{j}(idx)), hold on
    title('theta estimated')
    %legend('CD 250','CD 500','CD 750')

    idx = find(fpsumDVAll{j} ~= cdsumDVAll{j});
    minidx=min(idx)-1;
    maxidx=max(idx)+1;
    idx=[max(minidx,1):min(maxidx,length(dc))];
    figure(l+2*k-1)    
    subplot(3,2,4), scatter(dc(idx)*(max(sumDV{j}-min(sumDV{j})))+min(sumDV{j}),fpsumDVAll{j}(idx)), hold on
    title('sumDV')
%    legend('FP 250','FP 500','FP 750')
    figure(l+2*k)
    subplot(3,2,4), scatter(dc(idx)*(max(sumDV{j}-min(sumDV{j})))+min(sumDV{j}),cdsumDVAll{j}(idx)), hold on
    title('sumDV')
%    legend('CD 250','CD 500','CD 750')    
    
    idx = find(fplogitAll{j} ~= cdlogitAll{j});
    minidx=min(idx)-1;
    maxidx=max(idx)+1;
    idx=[max(minidx,1):min(maxidx,length(dc))];
    figure(l+2*k-1)    
    subplot(3,2,5), scatter(dc(idx),fplogitAll{j}(idx)), hold on
    if(j==3)
        legend('FP 250','FP 750','FP Control')
    end
    title('weighted sum')
    xlabel('decision criterion')
    figure(l+2*k)
    subplot(3,2,5), scatter(dc(idx),cdlogitAll{j}(idx)), hold on
    title('weighted sum')
    if(j==3)
        legend('CD 250','CD 750','CD Control')    
    end
    xlabel('decision criterion')
end







hand=figure;
l=hand.Number-1;
for k=1:3
    figure(l+2*k-1)
    clf
    figure(l+2*k)
    clf
for j=1:3
    idx = find(fpDV{k,j} ~= cdDV{k,j});
    minidx=min(idx)-1;
    maxidx=max(idx)+1;
    idx=[max(minidx,1):min(maxidx,length(dc))];
    figure(l+2*k-1)
    dc=perfDV{k,j};
    subplot(3,2,1), scatter(dc(idx),fpDV{k,j}(idx)), hold on
    title(['Optimal with train = ',num2str(k)])
    %legend('FP 250','FP 500','FP 750')
    figure(l+2*k)
    subplot(3,2,1), scatter(dc(idx),cdDV{k,j}(idx)), hold on
    title(['Optimal with train = ',num2str(k)])
    %legend('CD 250','CD 500','CD 750')

    idx = find(fpDV8{k,j} ~= cdDV8{k,j});
    minidx=min(idx)-1;
    maxidx=max(idx)+1;
    idx=[max(minidx,1):min(maxidx,length(dc))];
    figure(l+2*k-1)
    dc=perfDV8{k,j};
    subplot(3,2,2), scatter(dc(idx),fpDV8{k,j}(idx)), hold on
    %legend('FP 250','FP 500','FP 750')
    title('Optimal8')
    figure(l+2*k)
    subplot(3,2,2), scatter(dc(idx),cdDV8{k,j}(idx)), hold on
    title('Optimal8')
    %legend('CD 250','CD 500','CD 750')

    idx = find(fpthetaDV{k,j} ~= cdthetaDV{k,j});
    minidx=min(idx)-1;
    maxidx=max(idx)+1;
    idx=[max(minidx,1):min(maxidx,length(dc))];
    figure(l+2*k-1)
    dc=perfthetaDV{k,j};
    subplot(3,2,3), scatter(dc(idx),fpthetaDV{k,j}(idx)), hold on
    title('theta estimated')
    %legend('FP 250','FP 500','FP 750')
    figure(l+2*k)
    subplot(3,2,3), scatter(dc(idx),cdthetaDV{k,j}(idx)), hold on
    title('theta estimated')
    %legend('CD 250','CD 500','CD 750')

    idx = find(fpsumDV{k,j} ~= cdsumDV{k,j});
    minidx=min(idx)-1;
    maxidx=max(idx)+1;
    idx=[max(minidx,1):min(maxidx,length(dc))];
    figure(l+2*k-1)    
    dc=perfsumDV{k,j};
    subplot(3,2,4), scatter(dc(idx)*(max(sumDV{j}-min(sumDV{j})))+min(sumDV{j}),fpsumDV{k,j}(idx)), hold on
    title('sumDV')
    %legend('FP 250','FP 500','FP 750')
    figure(l+2*k)
    subplot(3,2,4), scatter(dc(idx)*(max(sumDV{j}-min(sumDV{j})))+min(sumDV{j}),cdsumDV{k,j}(idx)), hold on
    title('sumDV')
    %legend('CD 250','CD 500','CD 750')    
    
    idx = find(fplogit{k,j} ~= cdlogit{k,j});
    minidx=min(idx)-1;
    maxidx=max(idx)+1;
    idx=[max(minidx,1):min(maxidx,length(dc))];
    figure(l+2*k-1)    
    dc=perflogit{k,j};
    subplot(3,2,5), scatter(dc(idx),fplogit{k,j}(idx)), hold on
    if(j==3)
        legend('FP 250','FP 750','FP Control')    
    end
    title('weighted sum')
    xlabel('Performance')
    figure(l+2*k)
    subplot(3,2,5), scatter(dc(idx),cdlogit{k,j}(idx)), hold on
    title('weighted sum')
    if(j==3)
        legend('CD 250','CD 750','CD Control')    
    end
    xlabel('Performance')
end
end





hand=figure;
l=hand.Number-1;
k=1;
figure(l+2*k-1)
clf
figure(l+2*k)
clf
for j=1:3
    idx = find(fpDVAll{j} ~= cdDVAll{j});
    minidx=min(idx)-1;
    maxidx=max(idx)+1;
    idx=[max(minidx,1):min(maxidx,length(dc))];
    figure(l+2*k-1)
    dc=perfDVAll{j};
    subplot(3,2,1), scatter(dc(idx),fpDVAll{j}(idx)), hold on
    title(['Optimal with train = All'])
    %legend('FP 250','FP 500','FP 750')
    figure(l+2*k)
    subplot(3,2,1), scatter(dc(idx),cdDVAll{j}(idx)), hold on
    title(['Optimal with train = All'])
    %legend('CD 250','CD 500','CD 750')

    idx = find(fpDV8All{j} ~= cdDV8All{j});
    minidx=min(idx)-1;
    maxidx=max(idx)+1;
    idx=[max(minidx,1):min(maxidx,length(dc))];
    figure(l+2*k-1)
    dc=perfDV8All{j};
    subplot(3,2,2), scatter(dc(idx),fpDV8All{j}(idx)), hold on
    %legend('FP 250','FP 500','FP 750')
    title('Optimal8')
    figure(l+2*k)
    subplot(3,2,2), scatter(dc(idx),cdDV8All{j}(idx)), hold on
    title('Optimal8')
    %legend('CD 250','CD 500','CD 750')

    idx = find(fpthetaDVAll{j} ~= cdthetaDVAll{j});
    minidx=min(idx)-1;
    maxidx=max(idx)+1;
    idx=[max(minidx,1):min(maxidx,length(dc))];
    figure(l+2*k-1)
    dc=perfthetaDVAll{j};
    subplot(3,2,3), scatter(dc(idx),fpthetaDVAll{j}(idx)), hold on
    title('theta estimated')
    %legend('FP 250','FP 500','FP 750')
    figure(l+2*k)
    subplot(3,2,3), scatter(dc(idx),cdthetaDVAll{j}(idx)), hold on
    title('theta estimated')
    %legend('CD 250','CD 500','CD 750')

    idx = find(fpsumDVAll{j} ~= cdsumDVAll{j});
    minidx=min(idx)-1;
    maxidx=max(idx)+1;
    idx=[max(minidx,1):min(maxidx,length(dc))];
    figure(l+2*k-1)    
    dc=perfsumDVAll{j};
    subplot(3,2,4), scatter(dc(idx),fpsumDVAll{j}(idx)), hold on
    title('sumDV')
%    legend('FP 250','FP 500','FP 750')
    figure(l+2*k)
    subplot(3,2,4), scatter(dc(idx),cdsumDVAll{j}(idx)), hold on
    title('sumDV')
%    legend('CD 250','CD 500','CD 750')    
    
    idx = find(fplogitAll{j} ~= cdlogitAll{j});
    minidx=min(idx)-1;
    maxidx=max(idx)+1;
    idx=[max(minidx,1):min(maxidx,length(dc))];
    figure(l+2*k-1)
    dc=perflogitAll{j};
    subplot(3,2,5), scatter(dc(idx),fplogitAll{j}(idx)), hold on
    if(j==3)
        legend('FP 250','FP 750','FP Control')
    end
    xlabel('performance')
    title('weighted sum')
    figure(l+2*k)
    subplot(3,2,5), scatter(dc(idx),cdlogitAll{j}(idx)), hold on
    title('weighted sum')
    if(j==3)
        legend('CD 250','CD 750','CD Control')    
    end
    xlabel('performance')
end


