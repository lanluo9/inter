
clear all
close all
load behavior_11mice
% 
% 
% 1) ID stores the mouse ID, it has total of 4 mice;
% 2) offs are store from 250, 500, 750;
% 3) Orien (n_mice x orientation ), start from 0 to 90 degree
% 4) Data(n_mice x offs x orientations x (yes/no))
% For instance: 
% First mouse FA number for 250 ms ISI  = Data(1,1,1,1)
%             CR number for 250 ms ISI = Data(1,1,1,2)
% First mouse Hit number for 90 degree 250 ms ISI  = Data(1,1,6,1)
%             Miss number for 90 degree 250 ms ISI = Data(1,1,6,2)
%
nummice=11;

for i=1:nummice
for j=1:3
    NTemp(i,j,:)=squeeze((Data(i,j,:,1)+Data(i,j,:,2)));
    temp(i,j,:)=squeeze(Data(i,j,:,1)./(Data(i,j,:,1)+Data(i,j,:,2)));    
    pp=spline(Orien(i,:),squeeze(temp(i,j,:)));
    prtemp(i,j,:)=ppval(pp,[0:22.5:90]);    
end
end

premp=squeeze(mean(prtemp));
prstd=squeeze(std(prtemp))/sqrt(nummice);
clear prtemp
premp=premp([1,3],:);



figure
cc='bry';
for i=1:nummice
for j=1:3
    subplot(4,3,i), plot(Orien(i,:),squeeze(Data(i,j,:,1)./(Data(i,j,:,1)+Data(i,j,:,2)))), hold on
end
    hold off
    title(['mouse #',num2str(i)])
    axis([0 90 0 1])
    legend('250','500','750')
end

figure
cc='bry';
for i=1:nummice
for j=1:3
    plot(Orien(i,:),squeeze(Data(i,j,:,1)./(Data(i,j,:,1)+Data(i,j,:,2))),strcat(cc(j),'-')), hold on
end
end
    hold off
%    title(['mouse #',num2str(i),' at ',num2str(offs(j))])
    axis([0 90 0 1])
    legend('250','500','750')
figure
plot([0:22.5:90],premp(1,:)), hold on
plot([0:22.5:90],premp(2,:)), hold off
axis([0 90 0 1])
legend('250','750')


figure
errorbar([0:22.5:90],premp(1,:),prstd(1,:)/sqrt(nummice)), hold on
errorbar([0:22.5:90],premp(2,:),prstd(2,:)/sqrt(nummice)), hold off
axis([0 90 0 1])
title('Pr detect')
xlabel('Orientation Change')
legend('250','750')


title('interpolated average')

filename{1} = '170323_i689_runs-002-003';%_testResp_adapt.mat';
filename{2} = '170323_i696_runs-002-003';%_testResp_adapt.mat';  
filename{3} = '170324_i674_runs-002-003';%_testResp_adapt.mat';  
filename{4} = '170327_i684_runs-002-003';%_testResp_adapt.mat';  
filename{5} = '170503_i711_runs-002-003';%_testResp_adapt.mat';  
filename{6} = '170503_i712_runs-002-003';%_testResp_adapt.mat';  
filename{7} = '170510_i574_runs-002-003';%_testResp_adapt.mat';
filename{8} = '170808_i720_runs-002-003';%_testResp_adapt.mat';
filename{9} = '170810_i738_runs-002-003';%_testResp_adapt.mat';
filename{10} = '170811_i739_runs-002-003';%_testResp_adapt.mat';
filename{11} = '170816_i745_runs-002-003';%_testResp_adapt.mat';
filename{12} = '170826_i746_runs-002-003';%_testResp_adapt.mat';

for n=1:length(filename)
for j=1:7
    data{n,j}.X=[];
    data{n,j}.Y=[];
end
end
NRmax=15;
for n=1:length(filename)

    load([filename{n},'_newFits.mat']);
%    load([filename{n},'_fits'])
    [m,loc]=max(ori_fit);
    
    
    
    prefs{n}=loc/180*2*pi;
%    if(exist('theta_90','var'))
        idxn{n}=find(theta_90<22.5);% & -cos(prefs{n})<0.5);
       ['Dataset ',num2str(n),' has ',num2str(length(idxn{n})),' good units using theta_90<22.5!']
        if(length(idxn{n})>NRmax)
%            theta_90(-cos(prefs{n})>0.5)=Inf;
            [m,idxn{n}]=sort(theta_90,'ascend');            
%            [m,idxn{n}]=sort(theta_90,'descend');            
            idxn{n}=idxn{n}(1:NRmax);
%            idxn{n}=idxn{n}(length(idxn{n})-NRmax+1:end);
%            idxn{n}=idxn{n}(randperm(length(idxn{n}),NRmax));
            ['Using NRmax = ',num2str(NRmax),' units']
        end
        
%    else
%        idxn=find(max(ori_fit)-min(ori_fit)>0.05*max(max(ori_fit)));
%       ['Dataset ',num2str(n),' has ',num2str(length(idxn)),' good units using Amplitude!']
%    end            
    prefs{n}=prefs{n}(idxn{n});
    f{n}=ori_fit(:,idxn{n});
    kappa{n}=abs(fft(log(f{n})));
    kappa{n}=kappa{n}(end,:);
    
    SNR{n}=(sqrt(kappa{n}).*max(f{n}))';
    
    for j=1:3

        for k=1:8   % uses one neurons in ind
            data{n,j}.X = [data{n,j}.X;ppResp{j,k}';];
%            data{n,j}.X = [data{n,j}.X;shuffle(ppResp{j,k}');];
            data{n,j}.Y = [data{n,j}.Y;k*ones(size(ppResp{j,k},2),1);];
        end
        
        idx=~any(isnan(data{n,j}.X),2);        
        data{n,j}.X=data{n,j}.X(idx,:);
        data{n,j}.Y=data{n,j}.Y(idx,1);       
        data{n,j}.X=data{n,j}.X(:,idxn{n});
        for k=1:8
            mu{n,j}(k,:)=mean(data{n,j}.X(data{n,j}.Y==k,:));
        end        
    end
    Am{n}=sqrt(sum((mu{n,2}-mu{n,1}).^2,1));
end



temp=data;
clear data
for n=1:length(filename)
%     data{n}.X = [temp{n,1}.X;temp{n,2}.X;];
%     data{n}.Y = [temp{n,1}.Y;temp{n,2}.Y;];
%     data{n}.cond = [ones(size(temp{n,1}.Y));2*ones(size(temp{n,2}.Y));];
    data{n}.X = [temp{n,1}.X;temp{n,2}.X;temp{n,3}.X;];
    data{n}.Y = [temp{n,1}.Y;temp{n,2}.Y;temp{n,3}.Y;];
    data{n}.cond = [ones(size(temp{n,1}.Y));2*ones(size(temp{n,2}.Y));3*ones(size(temp{n,3}.Y));];
%    data{n}.X = bsxfun(@plus,data{n}.X,-mean(data{n}.X));
    
    V{n}=eye(size(data{n}.X,2));
    D{n}=(std(data{n}.X)').^2;
    data{n}.X = bsxfun(@times,data{n}.X,1./std(data{n}.X));
    data{n}.Am = Am{n}; %is 1xnumber of cells
end

DATA=temp;
clear temp
    BLall = [];
    BRall = [];
    Ball = [];
    BAall = [];
    wall = [];
    prefsall = [];
    SNRall = [];
    kappaall = [];    
    Amall=[];
    datasetall=[];
%    usedatasets=[2:4,6:12];
usedatasets=[2:4,6:12];
for n=usedatasets
    

    %initialize 
    X=data{n}.X;
    Y=data{n}.Y;        
    cond=data{n}.cond;
    
    lidx = (Y>=4);
    ridx = (Y<=4 | Y==8);
        
    idx=(cond==2|cond==1|cond==3);%&(Y==7|Y==8|Y==1|Y==2|Y==6);
    idx=(cond==1|cond==2);
%    idx=cond==3;
    [BL{n},dev,stats] = glmfit(X(lidx&idx,:),Y(lidx&idx)~=8,'binomial','link','logit');
    [BR{n},dev,stats] = glmfit(X(ridx&idx,:),Y(ridx&idx)~=8,'binomial','link','logit');
    [B{n},dev,stats] = glmfit(X(idx,:),Y(idx)~=8,'binomial','link','logit');

    PR{n}=glmval(B{n},X(idx,:),'logit');
    
    idx=(cond==1);%&(Y==7|Y==8|Y==1|Y==2|Y==6);
    [BA{n},dev,stats] = glmfit(X(idx,:),Y(idx)~=8,'binomial','link','logit');
    
    

    w=B{n}(2:end);%zeros(size(X,2),1);
    b=B{n}(1);
    alpha=4;
    
    clear idx prc dprcDV
    ds = 22.5/180*pi;
    for j=1:2
    for k=1:5
        switch k
            case 1
                idx{j,k}=(Y==8) & (cond==j);
            case 2
                idx{j,k}=(Y==7|Y==1) & (cond==j);
%                idx{j,k}=(Y==7) & (cond==j);
%                 idx{j,k}=(Y==1) & (cond==j);
            case 3
                idx{j,k}=(Y==6|Y==2) & (cond==j);
%                idx{j,k}=(Y==6) & (cond==j);
%                idx{j,k}=(Y==2) & (cond==j);
            case 4
                idx{j,k}=(Y==5|Y==3) & (cond==j);
%                idx{j,k}=(Y==5) & (cond==j);
%                idx{j,k}=(Y==3) & (cond==j);
            case 5
                idx{j,k}=(Y==4) & (cond==j);
        end
    end
    end
    dv=0:0.001:1;
    for j=1:2
    for k=1:5
        DVLR=glmval(B{n},X(idx{j,k},:),'logit');
        pLR(n,j,k)=mean(DVLR>0.875);
        if(k==1)
            FP=mean(DVLR>dv);
            CD=mean(DVLR>dv);
        else
            CD=mean(DVLR>dv);
        end
        AUROC(n,j,k)=-trapz(FP,CD);
        pDetect(n,j,k)=mean(DVLR>mean(Y(~(cond==3))~=8));
%        %pLR(j,k)=mean(log(glmval(B{n},X(idx{j,k},:)),'logit'));
    end
    end
    DVLRsave{n}=X*B{n}(2:end,1)+B{n}(1);
    
        %START
    DF=-Inf;
    F=Inf;
    tol=1e-12;
    eps=0.2;
    iters=0;
    maxiters=50000;
    while ((DF < -tol & iters < maxiters) | iters<1000)
        iters=iters+1;
        Fold = F;
        DV=(X*w+b);
%        signDV=sign(DV);
%        DV=abs(DV);
        pr = 1./(1+exp(-DV));
        dprdDV = 1./(1+exp(-DV))./(1+exp(DV));
%        pr = 1./(1+exp(-DV.^2));
%        dprdDV = 1./(1+exp(-DV.^2))./(1+exp(DV.^2)).*DV*2;
 
        F=0;
        dKLdalpha=0;
        dKLdb=0;        
        dKLdw=zeros(size(X,2),1);
        
        for j=1:2
        for k=1:5
            
            prc = mean(pr(idx{j,k}));
            F = F - alpha*premp(j,k)*log(prc); 
            F = F - (1-premp(j,k))*log(1-prc^alpha);
            
            dKLdalpha = dKLdalpha - premp(j,k)*log(prc);
            dKLdalpha = dKLdalpha + (1-premp(j,k))/(1-prc^alpha)*log(prc)*prc^alpha;
            
            dKLdprc = -alpha*premp(j,k)/prc;
            dKLdprc = dKLdprc + alpha*(1-premp(j,k))/(1-prc^alpha)*prc^(alpha-1);
             
            dprcdb = mean(dprdDV(idx{j,k}));
            dprcdw = X(idx{j,k},:)'*dprdDV(idx{j,k})/sum(idx{j,k});
            
            dKLdb = dKLdb + dKLdprc*dprcdb;            
            dKLdw = dKLdw + dKLdprc*dprcdw;
        end
        end
        
        w = w - eps*dKLdw;
%        w = w - eps*sum(dKLdw)*w;
        b = b - eps*dKLdb;
        alpha = alpha - eps*dKLdalpha;
        DF=F-Fold;
    end
        
%    B{n}=B{n})/norm(B{n});
    B{n}(2:end)=diag(1./sqrt(D{n}))*B{n}(2:end)
    wsave{n}=diag(1./sqrt(D{n}))*w;    
    alpahsave(n)=alpha;
    bsave(n)=b;
    corsave(n)=corr(w,B{n}(2:end));
    behDVsave{n}=DV;
    figure(3)
    scatter(B{n}(2:end),wsave{n}), hold on
    
    for j=1:2
    for k=1:5
        prbar(n,j,k) = mean(pr(idx{j,k})).^alpha;
        behprsave(n,j,k) = mean(pr(idx{j,k}));
    end
    end
    figure
    plot([0:22.5:90],squeeze(prbar(n,1,:)),'r-o',[0:22.5:90],squeeze(prbar(n,2,:)),'b-o'), hold on
    plot([0:22.5:90],premp(1,:),'r-*',[0:22.5:90],premp(2,:),'b-*'), hold off
    title(['Dataset ',num2str(n),' with ',num2str(length(idxn{n})),' units and corr = ',num2str(corsave(n))])
    legend('model 250','model 750','behavior 250','behavior 750')

    
    
%    figure(20)
%    scatter((cos(prefs{n})*(cos(ds)-1)+sin(prefs{n})*(sin(ds))).*kappa{n} , wsave{n}')
%    plot([-pi:0.01:pi],((1-cos(ds))*cos([-pi:0.01:pi])-sin(ds)*sin([-pi:0.01:pi])*max(wsave{n})))
%    scatter(mod(prefs{n}-pi,2*pi)-pi,wsave{n}), hold on
    drawnow
    
    if(sum(usedatasets==n))
    BLall = [BLall;BL{n}(2:end);];
    BRall = [BRall;BR{n}(2:end);];
    Ball = [Ball;B{n}(2:end);];
    BAall = [BAall;BA{n}(2:end);];
    wall = [wall;wsave{n};];
    datasetall = [datasetall;size(wall)*n;];
    prefsall = [prefsall;prefs{n}';];
    SNRall = [SNRall;SNR{n};];
    kappaall = [kappaall;kappa{n}';];
    Amall = [Amall;data{n}.Am'];
    end
end


figure
%behprsave=behprsave(usedatasets,:,:);
behprsave=prbar(usedatasets,:,:);
LRprsave=pDetect(usedatasets,:,:);

errorbar([0:22.5:90],squeeze(mean(behprsave(:,1,:))),squeeze(std(behprsave(:,1,:)))/sqrt(length(usedatasets)),'r'), hold on
errorbar([0:22.5:90],squeeze(mean(behprsave(:,2,:))),squeeze(std(behprsave(:,2,:)))/sqrt(length(usedatasets)),'b')
errorbar([0:22.5:90],squeeze(mean(LRprsave(:,1,:))),squeeze(std(LRprsave(:,1,:)))/sqrt(length(usedatasets)),'r--'), hold on
errorbar([0:22.5:90],squeeze(mean(LRprsave(:,2,:))),squeeze(std(LRprsave(:,2,:)))/sqrt(length(usedatasets)),'b--')


prtot = squeeze(mean(prbar(usedatasets,:,:)));
prstd = squeeze(std(prbar(usedatasets,:,:)));
    figure
%    plot([0:22.5:90],squeeze(prtot(1,:)),'r-o',[0:22.5:90],squeeze(prtot(2,:)),'b-o'), hold on
    errorbar([0:22.5:90],prtot(1,:),prstd(1,:),'r'), hold on
    errorbar([0:22.5:90],prtot(2,:),prstd(2,:),'b')
    plot([0:22.5:90],premp(1,:),'r-*',[0:22.5:90],premp(2,:),'b-*'), hold off
    title('All Datasets')
    legend('model 250','model 750','behavior 250','behavior 750')
    drawnow

behavPC(1) = ((1-prtot(1,1))*7 + sum(prtot(1,2:end))/4)/8;    
behavPC(2) = ((1-prtot(2,1))*7 + sum(prtot(2,2:end))/8)/8    
    
    pDtot=squeeze(mean(pDetect(usedatasets,:,:)));
    pDstd=squeeze(std(pDetect(usedatasets,:,:)));
    figure
%    plot([0:22.5:90],squeeze(prtot(1,:)),'r-o',[0:22.5:90],squeeze(prtot(2,:)),'b-o'), hold on
    errorbar([0:22.5:90],pDtot(1,:),pDstd(1,:)/sqrt(length(usedatasets)),'r'), hold on
    errorbar([0:22.5:90],pDtot(2,:),pDstd(2,:)/sqrt(length(usedatasets)),'b')
    plot([0:22.5:90],premp(1,:),'r*',[0:22.5:90],premp(2,:),'b*'), hold off
    title('All Datasets')
    legend('LR 250','LR 750','behavior 250','behavior 750')
    axis([0 90 0 1])
    drawnow
    

    LRPC(1) = ((1-pDtot(1,1))*7 + sum(pDtot(1,2:end))/4)/8;
    LRPC(2) = ((1-pDtot(2,1))*7 + sum(pDtot(2,2:end))/4)/8
    figure
    subplot(2,1,1), hold on
    subplot(2,1,2), hold on
    for n=usedatasets
        subplot(2,1,1), plot([0:22.5:90],squeeze(pDetect(n,1,:)))
        subplot(2,1,2), plot([0:22.5:90],squeeze(pDetect(n,2,:)))
    end
    subplot(2,1,1), plot([0:22.5:90],premp(1,:),'k*'), hold off
    xlabel('Orientation Difference')
    ylabel('Detection Probability')
    title('Adapted (250ms)')
    subplot(2,1,2), plot([0:22.5:90],premp(2,:),'k*'), hold off
    xlabel('Orientation Difference')
    ylabel('Detection Probability')
    title('Un-adapted (750ms)')
    
    figure
    hold on
    for n=usedatasets
        plot([0:22.5:90],squeeze(pDetect(n,1,:)+pDetect(n,2,:))/2)
    end
    plot([0:22.5:90],(premp(1,:)+premp(2,:))/2,'k*'), hold off
    xlabel('Orientation Difference')
    ylabel('Detection Probability')
    title('Average Performance')
    
    
    perf(:,1) = sum(squeeze(pDetect(:,1,2:5)),2);
    perf(:,2) = sum(squeeze(pDetect(:,2,2:5)),2);
    perf(:,1) = perf(:,1) + 1-squeeze(pDetect(:,1,1));
    perf(:,2) = perf(:,2) + 1-squeeze(pDetect(:,2,1));
    perf=perf(usedatasets,:)/5;
    
    
    
    tempD=pDetect(usedatasets,:,:);
    [n1,n2,n3]=size(tempD);
    tempD = reshape(tempD,[n1*n2,n3]);
    
    figure
    errorbar([0:22.5:90],mean(tempD),std(tempD)/sqrt(length(usedatasets)),'k')
    hold on
    plot([0:22.5:90],premp(1,:),'r*',[0:22.5:90],premp(2,:),'b*'), hold off
    title('All Datasets')
    legend('LR','behavior 250','behavior 750')
    axis([0 90 0 1])
    drawnow
       
    dv=[0:0.00001:1];
    for n=usedatasets
        Y=data{n}.Y~=8;
%            cond=logical(sum(data{n}.cond==[1,2],2));
%            behDVsave{n}=behDVsave{n}-mean(behDVsave{n});
%            behDVsave{n}=behDVsave{n}/std(behDVsave{n});
        for j=1:5
        for k=1:2
        notY=data{n}.Y==8 & data{n}.cond==1;
            switch j
                case 1
                    not8=8;
                case 2
                    not8=[1,7];
                case 3
                    not8=[2,6];
                case 4
                    not8=[3,5];
                case 5
                    not8=[4];
            end
            cond=(data{n}.cond==k & logical(sum(data{n}.Y==not8,2)));
            CDbeh(n,k,:)=sum(1./(1+exp(-behDVsave{n}(cond)))>dv)/sum(cond);
            FPbeh(n,k,:)=sum(1./(1+exp(-behDVsave{n}(notY)))>dv)/sum(notY);
            AUROCbeh(n,k,j)=-trapz(squeeze(FPbeh(n,k,:)),squeeze(CDbeh(n,k,:)));
            CDLR(n,k,:)=sum(1./(1+exp(-DVLRsave{n}(cond)))>dv)/sum(cond);
            FPLR(n,k,:)=sum(1./(1+exp(-DVLRsave{n}(notY)))>dv)/sum(notY);
            AUROCLR(n,k,j)=-trapz(squeeze(FPLR(n,k,:)),squeeze(CDLR(n,k,:)));
            if(j==1 & k==1)
                AUROCbeh(n,k,j)=0.5;
                AUROCLR(n,k,j)=0.5;
            end
%             if(j==1 & k==2)
%                 AUROCbeh(n,k,j)=0.5;
%                 AUROCLR(n,k,j)=0.5;
%                 
%             end
                
        end
        end
            
%            [m,loc]=max((sum(1./(1+exp(-behDVsave{n}(Y&cond)))>dv)+sum(1./(1+exp(-behDVsave{n}(~Y&cond)))<=dv))/sum(cond));
            [m,loc]=max(squeeze(CDbeh(n,k,:)*1/8+(1-FPbeh(n,k,:))*7/8));
            dcbeh(n)=dv(loc);
            locbeh(n)=loc;
%            dcbeh(n)=0.5%2785;%mean(Y);
%            [m,loc]=max((sum(1./(1+exp(-DVLRsave{n}(Y&cond)))>dv)+sum(1./(1+exp(-DVLRsave{n}(~Y&cond)))<=dv))/sum(cond));
            [m,loc]=max(squeeze(CDLR(n,k,:)*1/8+(1-FPLR(n,k,:))*7/8));
            dcLR(n)=dv(loc);
            locLR(n)=loc;
%            dcLR(n)=0.5;%mean(Y);

            for k=1:2
                cond=data{n}.cond==k;
%                pcbeh2(n,k)=(sum(1./(1+exp(-behDVsave{n}(Y&cond)))>dcbeh(n))+sum(1./(1+exp(-behDVsave{n}(~Y&cond)))<=dcbeh(n)))/sum(cond);
%                pcLR2(n,k)=(sum(1./(1+exp(-DVLRsave{n}(Y&cond)))>dcLR(n))+sum(1./(1+exp(-DVLRsave{n}(~Y&cond)))<=dcLR(n)))/sum(cond);
                pcbeh2(n,k)=(CDbeh(n,k,locbeh(n))*1/8+(1-FPbeh(n,k,locbeh(n)))*7/8);
                pcLR2(n,k)=(CDLR(n,k,locLR(n))*1/8+(1-FPLR(n,k,locLR(n)))*7/8);
            end
            
    end
    
    
    figure
    errorbar([0:22.5:90],mean(squeeze(AUROCbeh(usedatasets,1,:))),std(squeeze(AUROCbeh(usedatasets,1,:)))/sqrt(10))
    hold on
    errorbar([0:22.5:90],mean(squeeze(AUROCbeh(usedatasets,2,:))),std(squeeze(AUROCbeh(usedatasets,2,:)))/sqrt(10))
%    title('Behavior')

%    figure
    errorbar([0:22.5:90],mean(squeeze(AUROCLR(usedatasets,1,:))),std(squeeze(AUROCLR(usedatasets,1,:)))/sqrt(10))
    hold on
    errorbar([0:22.5:90],mean(squeeze(AUROCLR(usedatasets,2,:))),std(squeeze(AUROCLR(usedatasets,2,:)))/sqrt(10))
    title('LR')
    
    figure
%    scatter(pcLR2(usedatasets,:),pcbeh2(usedatasets,:))
    figure
    hist(pcLR2(usedatasets,:))
    
figure
scatter(Ball,wall)
xlabel('LRweights')
ylabel('Behavior Fit weights')

[m,idx1] = sort(-cos(prefsall));
[m,idx2] = sort(-cos(prefsall),'descend');
len=floor(length(prefsall)/2);

idx1=find(-cos(prefsall)<-cos(22.5/180*pi));
idx2=find(-cos(prefsall)>-cos(22.5/180*pi));

[H,p]=ttest(wall(idx1))
[H,p]=ttest(wall(idx2))
[H,p]=ttest(Ball(idx1))
[H,p]=ttest(Ball(idx2))


len=min(length(idx1),length(idx2));

%[m,idx1] = sort(-cos(prefsall));
[m,idx2] = sort(-cos(prefsall),'descend');
%idx1=idx1(randi(length(idx1),len,1));
idx2=idx2(1:len);
[H,p]=ttest(wall(idx2),wall(idx1),'tail','right')
[H,p]=ttest(Ball(idx2),Ball(idx1),'tail','right')

[H,p]=ttest(wall(idx1))
[H,p]=ttest(wall(idx2))
[H,p]=ttest(Ball(idx1))
[H,p]=ttest(Ball(idx2))


% idx1=find(-cos(prefsall)>-1/2);
% idx2=find(-cos(prefsall)<-1/2);
% idx1=idx1(randi(length(idx1),length(idx2),1));

figure
[H,p]=ttest(wall(idx2),wall(idx1))
scatter(-cos(prefsall),wall)
xlabel('PV weights (-cos(pref))')
ylabel('Behavioral weights')
title(['pv<0 vs pv>0 significance = ',num2str(p)])

[H,p]=ttest(Ball(idx2),Ball(idx1))
figure
scatter(-cos(prefsall),Ball)
xlabel('PV weights (-cos(pref))')
ylabel('LR weights')
title(['pv<0 vs pv>0 significance = ',num2str(p)])

figure
scatter(Ball,wall)
xlabel('LR weights')
ylabel('Behavioral weights')


figure
scatter((prefsall)/pi*90,wall)
xlabel('PV weights (-cos(pref))')
ylabel('Behavioral weights')
%title(['pv<0 vs pv>0 significance = ',num2str(p)])

[H,p]=ttest(Ball(idx2),Ball(idx1),'tail','left')
figure
scatter(prefsall/pi*90,Ball)
xlabel('PV weights (-cos(pref))')
ylabel('LR weights')

figure
scatter(prefsall/pi*90,Ball), hold on
scatter((prefsall)/pi*90,wall), hold off
xlabel('Orientation Preference')
ylabel('Weights')
legend('LR','Beh')

%title(['pv<0 vs pv>0 significance = ',num2str(p)])

% [m,idx1] = sort(cos(prefsall).*(1-cos(ds)-sin(prefsall)*sin(ds)));
% [m,idx2] = sort(cos(prefsall).*(1-cos(ds)-sin(prefsall)*sin(ds)),'descend');
% len=floor(length(prefsall)/2);
% idx1=idx1(1:len);
% idx2=idx2(1:len);
% [H,p]=ttest(wall(idx2),wall(idx1),'tail','left')

idx1=find( -cos(prefsall)>=1/2);
idx2=find( -cos(prefsall)<=-1/2);
idx3=find( abs(-cos(prefsall))<1/2);
%idx1=find( Ball>0);
%idx2=find( Ball<0);

figure
scatter(Ball(idx1),wall(idx1)),hold on,scatter(Ball(idx2),wall(idx2)),scatter(Ball(idx3),wall(idx3)), hold off
mean(wall(idx1))
mean(wall(idx2))
mean(Ball(idx1))
mean(Ball(idx2))

%figure
%scatter(-cos(prefsall)*(1-cos(ds))+sin(prefsall)*sin(ds),wall)
%
%
%figure
%scatter(SNRall,abs(wall))
%scatter(SNRall,abs(Ball))
figure
scatter(-cos(prefsall),Ball)
hold on
scatter(-cos(prefsall),wall)
hold off



% figure
% scatter(prefsall,wall)
% 
% figure
% scatter(-cos(prefsall),Ball)

