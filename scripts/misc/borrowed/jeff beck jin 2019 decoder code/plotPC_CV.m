
clear all
close all

%not8 = [1:7];
not8 = [1,2,6,7];
%not8 = [2,3,4,5,6];
%not8 = [6,2];
%not8 = [1,7];
%not8 = [5,3]
%DV = DV.{ opt, opt8, est, JM, PV, PVemp, logit, sum }
NDC=200;
dv=[0:NDC]/NDC;

usedatasets =[2,3,4,5,6,7,8,9,11,12];
optdv=0.27*2/3;

for train =1:3
    load(['alldataAnalyzedTR',num2str(train)])
    figure
    hold on
    for j=1:3
        idxfp = (DV.Y==8 & DV.cond==j & logical(sum(DV.dataset==usedatasets,2)));
        idxcd = (logical(sum(DV.Y==not8,2)) & DV.cond==j & logical(sum(DV.dataset==usedatasets,2))) ;
        DVtemp = abs(getfield(DV,'opt'));
        
        PCtemp{train,j} = (sum(DVtemp(idxcd)>dv)+sum(DVtemp(idxfp)<=dv))/(sum(idxfp)+sum(idxcd));
        [PC(train,j),DVopt(train,j)] = max(PCtemp{train,j}');
        DVopt(train,j)=dv(DVopt(train,j));
        
        %(sum(DVtemp(idxcd)>optdv)+sum(DVtemp(idxfp)<=optdv))/(sum(idxfp)+sum(idxcd));
        chance = sum(idxcd)/(sum(idxfp)+sum(idxcd))
        plot(dv,PCtemp{train,j})
        xlabel('DV'), ylabel('PC')
    end
end




train = 3;
load(['alldataAnalyzedTR',num2str(train)])
names = fieldnames(DV);
names = names(1:8);

for n=1:length(usedatasets);
    
    figure
    for k=1:length(names)
        subplot(3,3,k)
    for j=1:3
        idxfp = (DV.Y==8 & DV.cond==j & DV.dataset==usedatasets(n));
        idxcd = (logical(sum(DV.Y==not8,2)) & DV.cond==j & DV.dataset==usedatasets(n)) ;
        DVtemp = getfield(DV,names{k});
        DVtemp=abs(DVtemp);
        CD = mean(DVtemp(idxcd)>dv);
        FP = mean(DVtemp(idxfp)>dv);
        
        scatter(dv,FP./CD), hold on
        xlabel('DV'), ylabel('FP/CD')
    end
    legend('250','750','control')
    title([names{k},'  dataset = ',num2str(usedatasets(n))])
    end
end
   
    figure
    for k=1:length(names)
        subplot(3,3,k)
    for j=1:3
        idxfp = (DV.Y==8 & DV.cond==j & logical(sum(DV.dataset==usedatasets,2)));
        idxcd = (logical(sum(DV.Y==not8,2)) & DV.cond==j & logical(sum(DV.dataset==usedatasets,2))) ;
        DVtemp = getfield(DV,names{k});
        DVtemp=abs(DVtemp);
        CD = mean(DVtemp(idxcd)>dv);
        FP = mean(DVtemp(idxfp)>dv);
        
        scatter(dv,FP./CD), hold on
        xlabel('DV'), ylabel('FP/CD')
    end
    legend('250','750','control')
    title([names{k},'  dataset = ALL'])
    end


    
        figure
    for k=1:length(names)
        subplot(3,3,k)
    for j=1:3
        idxfp = (DV.Y==8 & DV.cond==j & logical(sum(DV.dataset==usedatasets,2)));
        idxcd = (logical(sum(DV.Y==not8,2)) & DV.cond==j & logical(sum(DV.dataset==usedatasets,2))) ;
        DVtemp = getfield(DV,names{k});
        DVtemp=abs(DVtemp);
        CD = mean(DVtemp(idxcd)>dv);
        FP = mean(DVtemp(idxfp)>dv);
        
        scatter(dv,FP./CD), hold on
        xlabel('DV'), ylabel('FP/CD')
    end
    legend('250','750','control')
    title([names{k},'  dataset = ALL'])
    end

    
names = fieldnames(DV);
names = names(1:8);
    

numdatasets=12;

%usedatasets =[1:numdatasets];
usedatasets =[2,3,4,7,8,9,11,12];
figure
for k=1:length(names);
    subplot(3,3,k)
    for j=1:3
        idxfp = (DV.Y==8 & DV.cond==j & logical(sum(DV.dataset==usedatasets,2)));
        idxcd = (logical(sum(DV.Y==not8,2)) & DV.cond==j & logical(sum(DV.dataset==usedatasets,2))) ;
        DVtemp = getfield(DV,names{k});
%        if(names{k}=='est' | names{k}=='PV' | names{k}=='PVemp' | names{k}=='JM' )
            DVtemp=abs(DVtemp);
%        end
        CD = mean(DVtemp(idxcd)>dv);
        FP = mean(DVtemp(idxfp)>dv);
        scatter(FP,CD), hold on
        xlabel('FP'), ylabel('CD'), title([names{k},' control AUROC =',num2str(-trapz(FP,CD))])
    end
legend('250','750','control')
end

usedatasets =[2,3,4,7,8,9,11,12];
figure
for k=1:length(names)
    subplot(3,3,k)
    for j=1:3
        idxfp = (DV.Y==8 & DV.cond==j & logical(sum(DV.dataset==usedatasets,2)));
        idxcd = (logical(sum(DV.Y==not8,2)) & DV.cond==j & logical(sum(DV.dataset==usedatasets,2))) ;
        DVtemp = getfield(DV,names{k});
%        if(names{k}=='est' | names{k}=='PV' | names{k}=='PVemp' | names{k}=='JM' )
            DVtemp=abs(DVtemp);
%        end
        CD = mean(DVtemp(idxcd)>dv);
        FP = mean(DVtemp(idxfp)>dv);
        scatter(FP,CD), hold on
        xlabel('FP'), ylabel('CD'), title([names{k},' control AUROC =',num2str(-trapz(FP,CD))])
    end
legend('250','750','control')
end


usedatasets =[2,3,4,7,8,9,11,12];
for n=1:length(usedatasets)
figure
for k=1:length(names)
    subplot(3,3,k)
    for j=1:3
        idxfp = (DV.Y==8 & DV.cond==j & DV.dataset==usedatasets(n));
        idxcd = (logical(sum(DV.Y==not8,2)) & DV.cond==j & DV.dataset==usedatasets(n)) ;
        DVtemp = getfield(DV,names{k});
%        if(names{k}=='est' | names{k}=='PV' | names{k}=='PVemp' | names{k}=='JM' )
            DVtemp=abs(DVtemp);
%        end
        CD = mean(DVtemp(idxcd)>dv);
        FP = mean(DVtemp(idxfp)>dv);
        scatter(FP,CD), hold on
        xlabel('FP'), ylabel('CD'), title([names{k},' control AUROC =',num2str(-trapz(FP,CD))])
    end
legend('250','750','control')
end
end




figure
for k=1:length(names);
    subplot(3,3,k)
    for j=1:3
        idxfp = (DV.Y==8 & DV.cond==j & logical(sum(DV.dataset==usedatasets,2)));
        idxcd = (logical(sum(DV.Y==not8,2)) & DV.cond==j & logical(sum(DV.dataset==usedatasets,2))) ;
        DVtemp = getfield(DV,names{k});
%        if(names{k}=='est' | names{k}=='PV' | names{k}=='PVemp' | names{k}=='JM' )
            DVtemp=abs(DVtemp);
%        end
        CD = mean(DVtemp(idxcd)>dv);
        FP = mean(DVtemp(idxfp)>dv);
        scatter(dv,CD), hold on
        xlabel('DV'), ylabel('CD'), title(names{k})
    end
legend('250','750','control')
end

figure
for k=1:length(names);
    subplot(3,3,k)
    for j=1:3
        idxfp = (DV.Y==8 & DV.cond==j & logical(sum(DV.dataset==usedatasets,2)));
        idxcd = (logical(sum(DV.Y==not8,2)) & DV.cond==j & logical(sum(DV.dataset==usedatasets,2))) ;
        DVtemp = getfield(DV,names{k});
 %       if(names{k}=='est' | names{k}=='PV' | names{k}=='PVemp' | names{k}=='JM' )
            DVtemp=abs(DVtemp);
%        end
        CD = mean(DVtemp(idxcd)>dv);
        FP = mean(DVtemp(idxfp)>dv);
        scatter(dv,FP), hold on
        xlabel('DV'), ylabel('FP'), title(names{k})
    end
legend('250','750','control')
end


for l=1:numdatasets
    usedatasets=l;
    figure
for k=1:length(names);
    subplot(3,3,k)
    for j=1:3
        idxfp = (DV.Y==8 & DV.cond==j & DV.dataset==l);
        idxcd = (logical(sum(DV.Y==not8,2)) & DV.cond==j & DV.dataset==l) ;
        DVtemp = getfield(DV,names{k});
%        if(names{k}=='est' | names{k}=='PV' | names{k}=='PVemp' | names{k}=='JM' )
            DVtemp=abs(DVtemp);
%        end
        CD = mean(DVtemp(idxcd)>dv);
        FP = mean(DVtemp(idxfp)>dv);
        scatter(FP,CD), hold on
        xlabel('FP'), ylabel('CD'), title([names{k},' control AUROC =',num2str(-trapz(FP,CD))])
    end
legend('250','750','control')
end
end



for l=1:numdatasets
    usedatasets=l;
    figure
for k=1:length(names);
    subplot(3,3,k)
    for j=1:3
        idxfp = (DV.Y==8 & DV.cond==j & DV.dataset==l);
        idxcd = (logical(sum(DV.Y==not8,2)) & DV.cond==j & DV.dataset==l) ;
        DVtemp = getfield(DV,names{k});
 %       if(names{k}=='est' | names{k}=='PV' | names{k}=='PVemp' | names{k}=='JM' )
            DVtemp=abs(DVtemp);
%        end
        CD = mean(DVtemp(idxcd)>dv);
        FP = mean(DVtemp(idxfp)>dv);
        scatter(dv,CD), hold on
        xlabel('DV'), ylabel('CD'), title(names{k})
    end
legend('250','750','control')
end
end

for l=1:numdatasets
    usedatasets=l;
    figure
for k=1:length(names);
    subplot(3,3,k)
    for j=1:3
        idxfp = (DV.Y==8 & DV.cond==j & DV.dataset==l);
        idxcd = (logical(sum(DV.Y==not8,2)) & DV.cond==j & DV.dataset==l) ;
        DVtemp = getfield(DV,names{k});
%        if(names{k}=='est' | names{k}=='PV' | names{k}=='PVemp' | names{k}=='JM' )
            DVtemp=abs(DVtemp);
%        end
        CD = mean(DVtemp(idxcd)>dv);
        FP = mean(DVtemp(idxfp)>dv);
        scatter(dv,FP), hold on
        xlabel('DV'), ylabel('FP'), title(names{k})
    end
legend('250','750','control')
end
end

figure(1)
figure(2)
figure(3)

ntheta=size(DV.pr,2);
theta=[1:ntheta]/ntheta*180-90;
figure
hold on
symb{1}='o';
symb{2}='.';
symb{3}='-';
dataset=2;
for j=1:2
    idx=(DV.Y==8 & DV.cond==j & DV.dataset==dataset);
%    prbar=exp(mean(log(DV.pr(idx,:))));
    prbar=mean(DV.pr(idx,:));
    prbar=prbar/sum(prbar);
    plot(theta,prbar,strcat('b',symb{j}))
    idx=(DV.Y==7 & DV.cond==j & DV.dataset==dataset);
%    prbar=exp(mean(log(DV.pr(idx,:))));
    prbar=mean(DV.pr(idx,:));
    prbar=prbar/sum(prbar);
    plot(theta,prbar,strcat('r',symb{j}))
    idx=(DV.Y==1 & DV.cond==j & DV.dataset==dataset);
%    prbar=exp(mean(log(DV.pr(idx,:))));
    prbar=mean(DV.pr(idx,:));
    prbar=prbar/sum(prbar);
    plot(theta,prbar,strcat('g',symb{j}))
end
hold off
title('Epr:  o = 250, . = 750, - = control')

ntheta=size(DV.pr,2);
theta=[1:ntheta]/ntheta*180-90;
figure
hold on
symb{1}='o';
symb{2}='.';
symb{3}='-';
dataset=2;
for j=1:2
    idx=(DV.Y==8 & DV.cond==j & DV.dataset==dataset);
    prbar=exp(mean(log(DV.pr(idx,:))));
%    prbar=mean(DV.pr(idx,:));
    prbar=prbar/sum(prbar);
    plot(theta,prbar,strcat('b',symb{j}))
    idx=(DV.Y==7 & DV.cond==j & DV.dataset==dataset);
    prbar=exp(mean(log(DV.pr(idx,:))));
%    prbar=mean(DV.pr(idx,:));
    prbar=prbar/sum(prbar);
    plot(theta,prbar,strcat('r',symb{j}))
    idx=(DV.Y==1 & DV.cond==j & DV.dataset==dataset);
    prbar=exp(mean(log(DV.pr(idx,:))));
%    prbar=mean(DV.pr(idx,:));
    prbar=prbar/sum(prbar);
    plot(theta,prbar,strcat('g',symb{j}))
end
hold off
title('expElogpr: o = 250, . = 750, - = control')


% compute bais and variance of estimators.

clear theta
theta(1,1,:)=[0.25,0.5,0.75,-1,-0.75,-0.5,-0.25,0];
figure
gg=gcf;
gg=gg.Number;

for n=1:12
    for j=1:3
    for k=1:8
        idx=(DV.Y==k & DV.cond==j & DV.dataset==n);        
        
        bias_est(n,j,k)=angle(mean(exp(sqrt(-1)*DV.est(idx)*pi)))/pi;
        var_est(n,j,k)=1-mean(cos(pi*(DV.est(idx)-bias_est(n,j,k))));
        bias_est(n,j,k)=mod(bias_est(n,j,k)-theta(k)+1,2)-1+theta(k);
        
        bias_JM(n,j,k)=angle(mean(exp(sqrt(-1)*DV.JM(idx)*pi)))/pi;
        var_JM(n,j,k)=1-mean(cos(pi*(DV.JM(idx)-bias_JM(n,j,k))));
        bias_JM(n,j,k)=mod(bias_JM(n,j,k)-theta(k)+1,2)-1+theta(k);

        bias_PV(n,j,k)=angle(mean(exp(sqrt(-1)*DV.PV(idx)*pi)))/pi;
        var_PV(n,j,k)=1-mean(cos(pi*(DV.PV(idx)-bias_PV(n,j,k))));
        bias_PV(n,j,k)=mod(bias_PV(n,j,k)-theta(k)+1,2)-1+theta(k);

        bias_PVemp(n,j,k)=angle(mean(exp(sqrt(-1)*DV.PVemp(idx)*pi)))/pi;
        var_PVemp(n,j,k)=1-mean(cos(pi*(DV.PVemp(idx)-bias_PVemp(n,j,k))));
        bias_PVemp(n,j,k)=mod(bias_PVemp(n,j,k)-theta(k)+1,2)-1+theta(k);

%         temp=mod(DV.JM(idx)-(k/8*2-1)+1,2)-1+(k/8*2-1);        
%         bias_JM(n,j,k)=mean(temp);
%         var_JM(n,j,k)=var(temp);
% 
%         temp=mod(DV.PV(idx)-(k/8*2-1)+1,2)-1+(k/8*2-1);        
%         bias_PV(n,j,k)=mean(temp);
%         var_PV(n,j,k)=var(temp);
% 
%         temp=mod(DV.PVemp(idx)-(k/8*2-1)+1,2)-1+(k/8*2-1);        
%         bias_PVemp(n,j,k)=mean(temp);
%         var_PVemp(n,j,k)=var(temp);
    end
    dbias_est(n,j,:)=(circshift(bias_est(n,j,:)-theta,1)-circshift(bias_est(n,j,:)-theta,-1))*8;
    dbias_JM(n,j,:)=(circshift(bias_JM(n,j,:)-theta,1)-circshift(bias_JM(n,j,:)-theta,-1))*8;
    dbias_PV(n,j,:)=(circshift(bias_PV(n,j,:)-theta,1)-circshift(bias_PV(n,j,:)-theta,-1))*8;
    dbias_PVemp(n,j,:)=(circshift(bias_PVemp(n,j,:)-theta,1)-circshift(bias_PVemp(n,j,:)-theta,-1))*8;
    
    FI_est(n,j,:)=dbias_est(n,j,:).^2./var_est(n,j,:);
    FI_JM(n,j,:)=dbias_JM(n,j,:).^2./var_JM(n,j,:);
    FI_PV(n,j,:)=dbias_PV(n,j,:).^2./var_PV(n,j,:);
    FI_PVemp(n,j,:)=dbias_PVemp(n,j,:).^2./var_PVemp(n,j,:);
    
    if(j==3)
    figure(gg)
    subplot(2,2,1), scatter(theta,squeeze(bias_est(n,j,:))),title('est')
    hold on
    subplot(2,2,2), scatter(theta,squeeze(bias_JM(n,j,:))),title('JM')
    hold on
    subplot(2,2,3), scatter(theta,squeeze(bias_PV(n,j,:))),title('PV')
    hold on
    subplot(2,2,4), scatter(theta,squeeze(bias_PVemp(n,j,:))),title('PVemp')
    hold on
    
    
    figure(gg+1)
    subplot(2,2,1), scatter(theta,squeeze(var_est(n,j,:))),title('est')
    hold on
    subplot(2,2,2), scatter(theta,squeeze(var_JM(n,j,:))),title('JM')
    hold on
    subplot(2,2,3), scatter(theta,squeeze(var_PV(n,j,:))),title('PV')
    hold on
    subplot(2,2,4), scatter(theta,squeeze(var_PVemp(n,j,:))),title('PVemp')
    hold on
    
    
    figure(gg+2)
    subplot(2,2,1), scatter(theta,squeeze(FI_est(n,j,:))),title('est')
    hold on
    subplot(2,2,2), scatter(theta,squeeze(FI_JM(n,j,:))),title('JM')
    hold on
    subplot(2,2,3), scatter(theta,squeeze(FI_PV(n,j,:))),title('PV')
    hold on
    subplot(2,2,4), scatter(theta,squeeze(FI_PVemp(n,j,:))),title('PVemp')
    hold on
    
    
    
    end
    
    end
end

figure(gg)
usedatasets =[2,3,4,7,8,9,11,12];

subplot(2,2,1),refline(1)
subplot(2,2,2),refline(1)
subplot(2,2,3),refline(1)
subplot(2,2,4),refline(1)

figure
    subplot(2,2,1), scatter(theta,squeeze(mean(bias_est(usedatasets,1,:)))-squeeze(theta)),title('est'),hold on
    subplot(2,2,1), scatter(theta,squeeze(mean(bias_est(usedatasets,2,:)))-squeeze(theta)),title('est')
    subplot(2,2,1), scatter(theta,squeeze(mean(bias_est(usedatasets,3,:)))-squeeze(theta)),title('est')
%    refline(1)
    subplot(2,2,2), scatter(theta,squeeze(mean(bias_JM(usedatasets,1,:)))-squeeze(theta)),title('JM'),hold on
    subplot(2,2,2), scatter(theta,squeeze(mean(bias_JM(usedatasets,2,:)))-squeeze(theta)),title('JM')
    subplot(2,2,2), scatter(theta,squeeze(mean(bias_JM(usedatasets,3,:)))-squeeze(theta)),title('JM')
%    refline(1)
    subplot(2,2,3), scatter(theta,squeeze(mean(bias_PV(usedatasets,1,:)))-squeeze(theta)),title('PV'),hold on
    subplot(2,2,3), scatter(theta,squeeze(mean(bias_PV(usedatasets,2,:)))-squeeze(theta)),title('PV')
    subplot(2,2,3), scatter(theta,squeeze(mean(bias_PV(usedatasets,3,:)))-squeeze(theta)),title('PV')
%    refline(1)
    subplot(2,2,4), scatter(theta,squeeze(mean(bias_PVemp(usedatasets,1,:)))-squeeze(theta)),title('PVemp'),hold on
    subplot(2,2,4), scatter(theta,squeeze(mean(bias_PVemp(usedatasets,2,:)))-squeeze(theta)),title('PVemp')
    subplot(2,2,4), scatter(theta,squeeze(mean(bias_PVemp(usedatasets,3,:)))-squeeze(theta)),title('PVemp')
%    refline(1)
legend('250','750','control')

figure
    subplot(2,2,1), scatter(theta,squeeze(mean(FI_est(usedatasets,1,:)))),title('est'),hold on
    subplot(2,2,1), scatter(theta,squeeze(mean(FI_est(usedatasets,2,:)))),title('est')
    subplot(2,2,1), scatter(theta,squeeze(mean(FI_est(usedatasets,3,:)))),title('est')
%    refline(1)
    subplot(2,2,2), scatter(theta,squeeze(mean(FI_JM(usedatasets,1,:)))),title('JM'),hold on
    subplot(2,2,2), scatter(theta,squeeze(mean(FI_JM(usedatasets,2,:)))),title('JM')
    subplot(2,2,2), scatter(theta,squeeze(mean(FI_JM(usedatasets,3,:)))),title('JM')
%    refline(1)
    subplot(2,2,3), scatter(theta,squeeze(mean(FI_PV(usedatasets,1,:)))),title('PV'),hold on
    subplot(2,2,3), scatter(theta,squeeze(mean(FI_PV(usedatasets,2,:)))),title('PV')
    subplot(2,2,3), scatter(theta,squeeze(mean(FI_PV(usedatasets,3,:)))),title('PV')
%    refline(1)
    subplot(2,2,4), scatter(theta,squeeze(mean(FI_PVemp(usedatasets,1,:)))),title('PVemp'),hold on
    subplot(2,2,4), scatter(theta,squeeze(mean(FI_PVemp(usedatasets,2,:)))),title('PVemp')
    subplot(2,2,4), scatter(theta,squeeze(mean(FI_PVemp(usedatasets,3,:)))),title('PVemp')
%    refline(1)
legend('250','750','control')

% 
% 
% NDC=100;
% l=0;
% figure
% % not8 = [1,7];
% for k=1:3
% for j=1:3
%     l=l+1;
%     subplot(3,3,l) 
%     for n=1:NDC
%         dc(n)=(n-1)/(NDC-1);
%         fpDV{k,j}(n) = mean(DV{k,j}(data{j}.Y==8)>dc(n));
% %        cdDV{k,j}(n) = mean(DV{k,j}(data{j}.Y~=8)>dc(n));
% %        cdDV{k,j}(n) = mean(DV{k,j}(data{j}.Y==7 | data{j}.Y==1)>dc(n));
%         cdDV{k,j}(n) = mean(DV{k,j}(logical(sum(data{j}.Y==not8,2))) > dc(n));
%         ndist=sum(data{j}.Y==8);
%         ntarget=sum(sum(data{j}.Y==not8));
%         perfDV{k,j}(n) = ((1-fpDV{k,j}(n))*ndist + cdDV{k,j}(n)*ntarget)/(ndist+ntarget);        
%     end
%     plot(fpDV{k,j},cdDV{k,j},'-o')
%     hold on
%     
%     for n=1:NDC
%         fpDV8{k,j}(n) = mean(DV8{k,j}(data{j}.Y==8)>dc(n));
% %        cdDV8{k,j}(n) = mean(DV8{k,j}(data{j}.Y~=8)>dc(n));
%         cdDV8{k,j}(n) = mean(DV8{k,j}(logical(sum(data{j}.Y==not8,2)))>dc(n));
%         ndist=sum(data{j}.Y==8);
%         ntarget=sum(sum(data{j}.Y==not8));
%         perfDV8{k,j}(n) = ((1-fpDV8{k,j}(n))*ndist + cdDV8{k,j}(n)*ntarget)/(ndist+ntarget);        
%     end
%     plot(fpDV8{k,j},cdDV8{k,j},'-x')
%     
%     for n=1:NDC
%         fpthetaDV{k,j}(n) = mean(thetaDV{k,j}(data{j}.Y==8)>dc(n));
% %        cdthetaDV{k,j}(n) = mean(thetaDV{k,j}(data{j}.Y~=8)>dc(n)*2*pi/4);
%         cdthetaDV{k,j}(n) = mean(thetaDV{k,j}(logical(sum(data{j}.Y==not8,2)))>dc(n));
%         ndist=sum(data{j}.Y==8);
%         ntarget=sum(sum(data{j}.Y==not8));
%         perfthetaDV{k,j}(n) = ((1-fpthetaDV{k,j}(n))*ndist + cdthetaDV{k,j}(n)*ntarget)/(ndist+ntarget);        
%     end
%     plot(fpthetaDV{k,j},cdthetaDV{k,j},'-+')
% 
%     for n=1:NDC
%         fpsumDV{k,j}(n) = mean(sumDV{j}(data{j}.Y==8)>dc(n)*(max(sumDV{j}-min(sumDV{j})))+min(sumDV{j}));
% %        cdsumDV{k,j}(n) = mean(sumDV{j}(data{j}.Y~=8)>dc(n)*(max(sumDV{j}-min(sumDV{j})))+min(sumDV{j}));
%         cdsumDV{k,j}(n) = mean(sumDV{j}(logical(sum(data{j}.Y==not8,2)))>dc(n)*(max(sumDV{j}-min(sumDV{j})))+min(sumDV{j}));
%         ndist=sum(data{j}.Y==8);
%         ntarget=sum(sum(data{j}.Y==not8));
%         perfsumDV{k,j}(n) = ((1-fpsumDV{k,j}(n))*ndist + cdsumDV{k,j}(n)*ntarget)/(ndist+ntarget);        
%     end
%     plot(fpsumDV{k,j},cdsumDV{k,j},'-s')
%     
%     for n=1:NDC
%         fplogit{k,j}(n) = mean(logitDV{k,j}(data{j}.Y==8)>dc(n));
% %        cdlogit{k,j}(n) = mean(logitDV{k,j}(data{j}.Y~=8)>dc(n));
%         cdlogit{k,j}(n) = mean(logitDV{k,j}(logical(sum(data{j}.Y==not8,2)))>dc(n));        
%         ndist=sum(data{j}.Y==8);
%         ntarget=sum(sum(data{j}.Y==not8));
%         perflogit{k,j}(n) = ((1-fplogit{k,j}(n))*ndist + cdlogit{k,j}(n)*ntarget)/(ndist+ntarget);        
%     end
%     plot(fplogit{k,j},cdlogit{k,j},'-*')
%     plot(0:1,0:1,'k-')
%     xlabel('FP')
%     ylabel('CD')
%     title(['Train = ',num2str(k),' Test = ',num2str(j)])
%     legend('Optimal','8 Category opt','Theta est','Sum','logit','refline')
%     hold off
% 
% end
% end
% 
% 
% 
% NDC=100;
% l=0;
% figure
% for j=1:3
%     l=l+1;
%     subplot(3,1,l) 
%     for n=1:NDC
%         dc(n)=(n-1)/(NDC-1);
%         fpDVAll{j}(n) = mean(DVAll{j}(data{j}.Y==8)>dc(n));
% %        cdDVAllj}(n) = mean(DVAll{j}(data{j}.Y~=8)>dc(n));
% %        cdDVAll{j}(n) = mean(DVAll{j}(data{j}.Y==7 | data{j}.Y==1)>dc(n));
%         cdDVAll{j}(n) = mean(DVAll{j}(logical(sum(data{j}.Y==not8,2))) > dc(n));
%         ndist=sum(data{j}.Y==8);
%         ntarget=sum(sum(data{j}.Y==not8));
%         perfDVAll{j}(n) = ((1-fpDVAll{j}(n))*ndist + cdDVAll{j}(n)*ntarget)/(ndist+ntarget);        
%     end
%     plot(fpDVAll{j},cdDVAll{j},'-o')
%     hold on
%     
%     for n=1:NDC
%         fpDV8All{j}(n) = mean(DV8All{j}(data{j}.Y==8)>dc(n));
% %        cdDV8All{j}(n) = mean(DV8All{j}(data{j}.Y~=8)>dc(n));
%         cdDV8All{j}(n) = mean(DV8All{j}(logical(sum(data{j}.Y==not8,2)))>dc(n));
%         ndist=sum(data{j}.Y==8);
%         ntarget=sum(sum(data{j}.Y==not8));
%         perfDV8All{j}(n) = ((1-fpDV8All{j}(n))*ndist + cdDV8All{j}(n)*ntarget)/(ndist+ntarget);        
%     end
%     plot(fpDV8All{j},cdDV8All{j},'-x')
%     
%     for n=1:NDC
%         fpthetaDVAll{j}(n) = mean(thetaDVAll{j}(data{j}.Y==8)>dc(n));
% %        cdthetaDVAll{j}(n) = mean(thetaDVAll{j}(data{j}.Y~=8)>dc(n));
%         cdthetaDVAll{j}(n) = mean(thetaDVAll{j}(logical(sum(data{j}.Y==not8,2)))>dc(n));
%         ndist=sum(data{j}.Y==8);
%         ntarget=sum(sum(data{j}.Y==not8));
%         perfthetaDVAll{j}(n) = ((1-fpthetaDVAll{j}(n))*ndist + cdthetaDVAll{j}(n)*ntarget)/(ndist+ntarget);        
%     end
%     plot(fpthetaDVAll{j},cdthetaDVAll{j},'-+')
% 
%     for n=1:NDC
%         fpsumDVAll{j}(n) = mean(sumDV{j}(data{j}.Y==8)>dc(n)*(max(sumDV{j}-min(sumDV{j})))+min(sumDV{j}));
% %        cdsumDVAll{j}(n) = mean(sumDV{j}(data{j}.Y~=8)>dc(n)*(max(sumDV{j}-min(sumDV{j})))+min(sumDV{j}));
%         cdsumDVAll{j}(n) = mean(sumDV{j}(logical(sum(data{j}.Y==not8,2)))>dc(n)*(max(sumDV{j}-min(sumDV{j})))+min(sumDV{j}));
%         ndist=sum(data{j}.Y==8);
%         ntarget=sum(sum(data{j}.Y==not8));
%         perfsumDVAll{j}(n) = ((1-fpsumDVAll{j}(n))*ndist + cdsumDVAll{j}(n)*ntarget)/(ndist+ntarget);        
%     end
%     plot(fpsumDVAll{j},cdsumDVAll{j},'-s')
%     
%     for n=1:NDC
%         fplogitAll{j}(n) = mean(logitDVAll{j}(data{j}.Y==8)>dc(n));
% %        cdlogitAll{j}(n) = mean(logitDVAll{j}(data{j}.Y~=8)>dc(n));
%         cdlogitAll{j}(n) = mean(logitDVAll{j}(logical(sum(data{j}.Y==not8,2)))>dc(n));        
%         ndist=sum(data{j}.Y==8);
%         ntarget=sum(sum(data{j}.Y==not8));        
%         perflogitAll{j}(n) = ((1-fplogitAll{j}(n))*ndist + cdlogitAll{j}(n)*ntarget)/(ndist+ntarget);
%     end
%     plot(fplogitAll{j},cdlogitAll{j},'-*')
%     plot(0:1,0:1,'k-')
%     xlabel('FP')
%     ylabel('CD')
%     title(['Train = All',' Test = ',num2str(j)])
%     legend('Optimal','8 Category opt','Theta est','Sum','logit','refline')
%     hold off
% 
% end
% 
% 
% 
% 
% hand=figure;
% l=hand.Number-1;
% for k=1:3
%     figure(l+2*k-1)
%     clf
%     figure(l+2*k)
%     clf
% for j=1:3
%     idx = find(fpDV{k,j} ~= cdDV{k,j});
%     minidx=min(idx)-1;
%     maxidx=max(idx)+1;
%     idx=[max(minidx,1):min(maxidx,length(dc))];
%     figure(l+2*k-1)
%     subplot(3,2,1), scatter(dc(idx),fpDV{k,j}(idx)), hold on
%     title(['Optimal with train = ',num2str(k)])
%     %legend('FP 250','FP 500','FP 750')
%     figure(l+2*k)
%     subplot(3,2,1), scatter(dc(idx),cdDV{k,j}(idx)), hold on
%     title(['Optimal with train = ',num2str(k)])
%     %legend('CD 250','CD 500','CD 750')
% 
%     idx = find(fpDV8{k,j} ~= cdDV8{k,j});
%     minidx=min(idx)-1;
%     maxidx=max(idx)+1;
%     idx=[max(minidx,1):min(maxidx,length(dc))];
%     figure(l+2*k-1)
%     subplot(3,2,2), scatter(dc(idx),fpDV8{k,j}(idx)), hold on
%     %legend('FP 250','FP 500','FP 750')
%     title('Optimal8')
%     figure(l+2*k)
%     subplot(3,2,2), scatter(dc(idx),cdDV8{k,j}(idx)), hold on
%     title('Optimal8')
%     %legend('CD 250','CD 500','CD 750')
% 
%     idx = find(fpthetaDV{k,j} ~= cdthetaDV{k,j});
%     minidx=min(idx)-1;
%     maxidx=max(idx)+1;
%     idx=[max(minidx,1):min(maxidx,length(dc))];
%     figure(l+2*k-1)
%     subplot(3,2,3), scatter(dc(idx)*360/2/2,fpthetaDV{k,j}(idx)), hold on
%     title('theta estimated')
%     %legend('FP 250','FP 500','FP 750')
%     figure(l+2*k)
%     subplot(3,2,3), scatter(dc(idx)*360/2/2,cdthetaDV{k,j}(idx)), hold on
%     title('theta estimated')
%     %legend('CD 250','CD 500','CD 750')
% 
%     idx = find(fpsumDV{k,j} ~= cdsumDV{k,j});
%     minidx=min(idx)-1;
%     maxidx=max(idx)+1;
%     idx=[max(minidx,1):min(maxidx,length(dc))];
%     figure(l+2*k-1)    
%     subplot(3,2,4), scatter(dc(idx)*(max(sumDV{j}-min(sumDV{j})))+min(sumDV{j}),fpsumDV{k,j}(idx)), hold on
%     title('sumDV')
%     %legend('FP 250','FP 500','FP 750')
%     figure(l+2*k)
%     subplot(3,2,4), scatter(dc(idx)*(max(sumDV{j}-min(sumDV{j})))+min(sumDV{j}),cdsumDV{k,j}(idx)), hold on
%     title('sumDV')
%     %legend('CD 250','CD 500','CD 750')    
%     
%     idx = find(fplogit{k,j} ~= cdlogit{k,j});
%     minidx=min(idx)-1;
%     maxidx=max(idx)+1;
%     idx=[max(minidx,1):min(maxidx,length(dc))];
%     figure(l+2*k-1)    
%     subplot(3,2,5), scatter(dc(idx),fplogit{k,j}(idx)), hold on
%     if(j==3)
%         legend('FP 250','FP 750','FP Control')    
%     end
%     title('weighted sum')
%     xlabel('decision criterion')
%     figure(l+2*k)
%     subplot(3,2,5), scatter(dc(idx),cdlogit{k,j}(idx)), hold on
%     title('weighted sum')
%     xlabel('decision criterion')
%     if(j==3)
%         legend('CD 250','CD 750','CD Control')    
%     end
% end
% end
% 
% 
% 
% 
% 
% hand=figure;
% l=hand.Number-1;
% k=1;
% figure(l+2*k-1)
% clf
% figure(l+2*k)
% clf
% for j=1:3
%     idx = find(fpDVAll{j} ~= cdDVAll{j});
%     minidx=min(idx)-1;
%     maxidx=max(idx)+1;
%     idx=[max(minidx,1):min(maxidx,length(dc))];
%     figure(l+2*k-1)
%     subplot(3,2,1), scatter(dc(idx),fpDVAll{j}(idx)), hold on
%     title(['Optimal with train = All'])
%     %legend('FP 250','FP 500','FP 750')
%     figure(l+2*k)
%     subplot(3,2,1), scatter(dc(idx),cdDVAll{j}(idx)), hold on
%     title(['Optimal with train = All'])
%     %legend('CD 250','CD 500','CD 750')
% 
%     idx = find(fpDV8All{j} ~= cdDV8All{j});
%     minidx=min(idx)-1;
%     maxidx=max(idx)+1;
%     idx=[max(minidx,1):min(maxidx,length(dc))];
%     figure(l+2*k-1)
%     subplot(3,2,2), scatter(dc(idx),fpDV8All{j}(idx)), hold on
%     %legend('FP 250','FP 500','FP 750')
%     title('Optimal8')
%     figure(l+2*k)
%     subplot(3,2,2), scatter(dc(idx),cdDV8All{j}(idx)), hold on
%     title('Optimal8')
%     %legend('CD 250','CD 500','CD 750')
% 
%     idx = find(fpthetaDVAll{j} ~= cdthetaDVAll{j});
%     minidx=min(idx)-1;
%     maxidx=max(idx)+1;
%     idx=[max(minidx,1):min(maxidx,length(dc))];
%     figure(l+2*k-1)
%     subplot(3,2,3), scatter(dc(idx)*360/2/2,fpthetaDVAll{j}(idx)), hold on
%     title('theta estimated')
%     %legend('FP 250','FP 500','FP 750')
%     figure(l+2*k)
%     subplot(3,2,3), scatter(dc(idx)*360/2/2,cdthetaDVAll{j}(idx)), hold on
%     title('theta estimated')
%     %legend('CD 250','CD 500','CD 750')
% 
%     idx = find(fpsumDVAll{j} ~= cdsumDVAll{j});
%     minidx=min(idx)-1;
%     maxidx=max(idx)+1;
%     idx=[max(minidx,1):min(maxidx,length(dc))];
%     figure(l+2*k-1)    
%     subplot(3,2,4), scatter(dc(idx)*(max(sumDV{j}-min(sumDV{j})))+min(sumDV{j}),fpsumDVAll{j}(idx)), hold on
%     title('sumDV')
% %    legend('FP 250','FP 500','FP 750')
%     figure(l+2*k)
%     subplot(3,2,4), scatter(dc(idx)*(max(sumDV{j}-min(sumDV{j})))+min(sumDV{j}),cdsumDVAll{j}(idx)), hold on
%     title('sumDV')
% %    legend('CD 250','CD 500','CD 750')    
%     
%     idx = find(fplogitAll{j} ~= cdlogitAll{j});
%     minidx=min(idx)-1;
%     maxidx=max(idx)+1;
%     idx=[max(minidx,1):min(maxidx,length(dc))];
%     figure(l+2*k-1)    
%     subplot(3,2,5), scatter(dc(idx),fplogitAll{j}(idx)), hold on
%     if(j==3)
%         legend('FP 250','FP 750','FP Control')
%     end
%     title('weighted sum')
%     xlabel('decision criterion')
%     figure(l+2*k)
%     subplot(3,2,5), scatter(dc(idx),cdlogitAll{j}(idx)), hold on
%     title('weighted sum')
%     if(j==3)
%         legend('CD 250','CD 750','CD Control')    
%     end
%     xlabel('decision criterion')
% end
% 
% 
% 
% 
% 
% 
% 
% hand=figure;
% l=hand.Number-1;
% for k=1:3
%     figure(l+2*k-1)
%     clf
%     figure(l+2*k)
%     clf
% for j=1:3
%     idx = find(fpDV{k,j} ~= cdDV{k,j});
%     minidx=min(idx)-1;
%     maxidx=max(idx)+1;
%     idx=[max(minidx,1):min(maxidx,length(dc))];
%     figure(l+2*k-1)
%     dc=perfDV{k,j};
%     subplot(3,2,1), scatter(dc(idx),fpDV{k,j}(idx)), hold on
%     title(['Optimal with train = ',num2str(k)])
%     %legend('FP 250','FP 500','FP 750')
%     figure(l+2*k)
%     subplot(3,2,1), scatter(dc(idx),cdDV{k,j}(idx)), hold on
%     title(['Optimal with train = ',num2str(k)])
%     %legend('CD 250','CD 500','CD 750')
% 
%     idx = find(fpDV8{k,j} ~= cdDV8{k,j});
%     minidx=min(idx)-1;
%     maxidx=max(idx)+1;
%     idx=[max(minidx,1):min(maxidx,length(dc))];
%     figure(l+2*k-1)
%     dc=perfDV8{k,j};
%     subplot(3,2,2), scatter(dc(idx),fpDV8{k,j}(idx)), hold on
%     %legend('FP 250','FP 500','FP 750')
%     title('Optimal8')
%     figure(l+2*k)
%     subplot(3,2,2), scatter(dc(idx),cdDV8{k,j}(idx)), hold on
%     title('Optimal8')
%     %legend('CD 250','CD 500','CD 750')
% 
%     idx = find(fpthetaDV{k,j} ~= cdthetaDV{k,j});
%     minidx=min(idx)-1;
%     maxidx=max(idx)+1;
%     idx=[max(minidx,1):min(maxidx,length(dc))];
%     figure(l+2*k-1)
%     dc=perfthetaDV{k,j};
%     subplot(3,2,3), scatter(dc(idx),fpthetaDV{k,j}(idx)), hold on
%     title('theta estimated')
%     %legend('FP 250','FP 500','FP 750')
%     figure(l+2*k)
%     subplot(3,2,3), scatter(dc(idx),cdthetaDV{k,j}(idx)), hold on
%     title('theta estimated')
%     %legend('CD 250','CD 500','CD 750')
% 
%     idx = find(fpsumDV{k,j} ~= cdsumDV{k,j});
%     minidx=min(idx)-1;
%     maxidx=max(idx)+1;
%     idx=[max(minidx,1):min(maxidx,length(dc))];
%     figure(l+2*k-1)    
%     dc=perfsumDV{k,j};
%     subplot(3,2,4), scatter(dc(idx)*(max(sumDV{j}-min(sumDV{j})))+min(sumDV{j}),fpsumDV{k,j}(idx)), hold on
%     title('sumDV')
%     %legend('FP 250','FP 500','FP 750')
%     figure(l+2*k)
%     subplot(3,2,4), scatter(dc(idx)*(max(sumDV{j}-min(sumDV{j})))+min(sumDV{j}),cdsumDV{k,j}(idx)), hold on
%     title('sumDV')
%     %legend('CD 250','CD 500','CD 750')    
%     
%     idx = find(fplogit{k,j} ~= cdlogit{k,j});
%     minidx=min(idx)-1;
%     maxidx=max(idx)+1;
%     idx=[max(minidx,1):min(maxidx,length(dc))];
%     figure(l+2*k-1)    
%     dc=perflogit{k,j};
%     subplot(3,2,5), scatter(dc(idx),fplogit{k,j}(idx)), hold on
%     if(j==3)
%         legend('FP 250','FP 750','FP Control')    
%     end
%     title('weighted sum')
%     xlabel('Performance')
%     figure(l+2*k)
%     subplot(3,2,5), scatter(dc(idx),cdlogit{k,j}(idx)), hold on
%     title('weighted sum')
%     if(j==3)
%         legend('CD 250','CD 750','CD Control')    
%     end
%     xlabel('Performance')
% end
% end
% 
% 
% 
% 
% 
% hand=figure;
% l=hand.Number-1;
% k=1;
% figure(l+2*k-1)
% clf
% figure(l+2*k)
% clf
% for j=1:3
%     idx = find(fpDVAll{j} ~= cdDVAll{j});
%     minidx=min(idx)-1;
%     maxidx=max(idx)+1;
%     idx=[max(minidx,1):min(maxidx,length(dc))];
%     figure(l+2*k-1)
%     dc=perfDVAll{j};
%     subplot(3,2,1), scatter(dc(idx),fpDVAll{j}(idx)), hold on
%     title(['Optimal with train = All'])
%     %legend('FP 250','FP 500','FP 750')
%     figure(l+2*k)
%     subplot(3,2,1), scatter(dc(idx),cdDVAll{j}(idx)), hold on
%     title(['Optimal with train = All'])
%     %legend('CD 250','CD 500','CD 750')
% 
%     idx = find(fpDV8All{j} ~= cdDV8All{j});
%     minidx=min(idx)-1;
%     maxidx=max(idx)+1;
%     idx=[max(minidx,1):min(maxidx,length(dc))];
%     figure(l+2*k-1)
%     dc=perfDV8All{j};
%     subplot(3,2,2), scatter(dc(idx),fpDV8All{j}(idx)), hold on
%     %legend('FP 250','FP 500','FP 750')
%     title('Optimal8')
%     figure(l+2*k)
%     subplot(3,2,2), scatter(dc(idx),cdDV8All{j}(idx)), hold on
%     title('Optimal8')
%     %legend('CD 250','CD 500','CD 750')
% 
%     idx = find(fpthetaDVAll{j} ~= cdthetaDVAll{j});
%     minidx=min(idx)-1;
%     maxidx=max(idx)+1;
%     idx=[max(minidx,1):min(maxidx,length(dc))];
%     figure(l+2*k-1)
%     dc=perfthetaDVAll{j};
%     subplot(3,2,3), scatter(dc(idx),fpthetaDVAll{j}(idx)), hold on
%     title('theta estimated')
%     %legend('FP 250','FP 500','FP 750')
%     figure(l+2*k)
%     subplot(3,2,3), scatter(dc(idx),cdthetaDVAll{j}(idx)), hold on
%     title('theta estimated')
%     %legend('CD 250','CD 500','CD 750')
% 
%     idx = find(fpsumDVAll{j} ~= cdsumDVAll{j});
%     minidx=min(idx)-1;
%     maxidx=max(idx)+1;
%     idx=[max(minidx,1):min(maxidx,length(dc))];
%     figure(l+2*k-1)    
%     dc=perfsumDVAll{j};
%     subplot(3,2,4), scatter(dc(idx),fpsumDVAll{j}(idx)), hold on
%     title('sumDV')
% %    legend('FP 250','FP 500','FP 750')
%     figure(l+2*k)
%     subplot(3,2,4), scatter(dc(idx),cdsumDVAll{j}(idx)), hold on
%     title('sumDV')
% %    legend('CD 250','CD 500','CD 750')    
%     
%     idx = find(fplogitAll{j} ~= cdlogitAll{j});
%     minidx=min(idx)-1;
%     maxidx=max(idx)+1;
%     idx=[max(minidx,1):min(maxidx,length(dc))];
%     figure(l+2*k-1)
%     dc=perflogitAll{j};
%     subplot(3,2,5), scatter(dc(idx),fplogitAll{j}(idx)), hold on
%     if(j==3)
%         legend('FP 250','FP 750','FP Control')
%     end
%     xlabel('performance')
%     title('weighted sum')
%     figure(l+2*k)
%     subplot(3,2,5), scatter(dc(idx),cdlogitAll{j}(idx)), hold on
%     title('weighted sum')
%     if(j==3)
%         legend('CD 250','CD 750','CD Control')    
%     end
%     xlabel('performance')
% end
% 
% 
