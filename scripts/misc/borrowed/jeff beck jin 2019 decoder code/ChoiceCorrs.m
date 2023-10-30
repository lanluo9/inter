close all
clear all

load('D:\Dropbox\Lindsay\randIntBehavData.mat')
NE = length(expt);

prefs=[];
ccs=[];
scs=[];
B=[];
Ba=[];
SI=[];
for i=[1:9,11:NE]
   temp=expt(i).trialOutcome;
   [m,n]=size(temp);
   k=0;
   clear Y idx Ya
   for j=1:n
       if(isnan(temp{j}))
       else
           k=k+1;
           idx(k)=j;
           if(temp{j}=='s')
               Y(k,1)=1;
               Ya(k,1)=1;
           elseif(temp{j}=='f')
               Y(k,1)=1;
               Ya(k,1)=0;
           elseif(temp{j}=='m');
               Y(k,1)=0;
               Ya(k,1)=1;
           end
       end
   end
   X=expt(i).targetResp(:,idx)';
   X=[X;expt(i).lastBaseResp(:,idx)';];
   Y=[Y;zeros(length(idx),1)];
   Ya=[Ya;zeros(length(idx),1)];

   X=X(:,logical(expt(i).signifResponsiveCells));
   
   idx=find(~isnan(sum(X,2)));
   Y=Y(idx,1);   
   Ya=Ya(idx,1);   
   X=X(idx,:);

%    [m,n]=size(X);
%    if(n>20)
%        idx=randperm(n);
%        X=X(:,idx(1:20));
%    end
%    
   
   %idx=find(~isnan(expt(i).oriTuning));
   %X=X(:,idx);   
   ccs=[ccs;corr(Y,X)';];
   scs=[scs;corr(Ya,X)';];
   
   %temp=expt(i).oriTuning(idx);
   %prefs=[prefs;mod(temp+90,180)-90;];
   
   X=bsxfun(@plus,X,-mean(X));
   X=bsxfun(@times,X,1./std(X));
   C=eye(size(X,2));%C=inv(sqrtm(cov(X)));
   X=X*C;
   p=1;
   [temp1,dev1,stats1]=glmfit(X,Ya,'binomial');
   [temp2,dev2,stats2]=glmfit(X,Y,'binomial');
   idx=find(stats1.p>p|stats2.p>p);

   temp1(idx)=0;
   temp2(idx)=0;
   Ba=[Ba;C*temp1(2:end,1);];
   B=[B;C*temp2(2:end,1);];
   
   
   
   
   SI=[SI;expt(i).baseTargRatio;];
   
end
% figure(1)
% scatter(prefs,ccs)
% xlabel('Preferred Orientation')
% ylabel('Choice Correlation')
% 
% figure(2)
% scatter(prefs,B)
% xlabel('Preferred Orientation')
% ylabel('LR weights')

figure(3)

scatter(Ba(scs<0),ccs(scs<0)), hold on
scatter(Ba(scs>0),ccs(scs>0))
ylabel('Choice Correlation')
xlabel('Neural LR weights')
ax=axis;
hold on;plot([0,0],[ax(3),ax(4)]);hold off
refline(0)

figure(4)
scatter(Ba(scs<0),B(scs<0)), hold on
scatter(Ba(scs>0),B(scs>0))
ylabel('Behavioral LR weights')
xlabel('Neural LR weights')
ax=axis;
hold on;plot([0,0],[ax(3),ax(4)]);hold off
refline(0)
refline(1)

figure(5)
scatter(B(scs<0),ccs(scs<0)), hold on
scatter(B(scs>0),ccs(scs>0))
xlabel('Behavioral LR weights')
ylabel('Choice Correlations')
ax=axis;
hold on;plot([0,0],[ax(3),ax(4)]);hold off
refline(0)

figure(6)
scatter(scs(scs<0),ccs(scs<0)), hold on
scatter(scs(scs>0),ccs(scs>0))
ylabel('Choice Correlation')
xlabel('Stimulus Correlation')
ax=axis;
hold on;plot([0,0],[ax(3),ax(4)]);hold off
refline(0)
refline(1)

figure(7)
scatter(Ba,B-Ba)
ylabel('Behavioral - Actual LR weights')
xlabel('Actual')
ax=axis;
hold on;plot([0,0],[ax(3),ax(4)]);hold off
refline(0)
[h,p]=ttest(B(Ba<0)-Ba(Ba<0))
mean(B(Ba<0)-Ba(Ba<0))
[h,p]=ttest(B(Ba>0)-Ba(Ba>0))
mean(B(Ba>0)-Ba(Ba>0))

[h,p]=ttest(ccs(scs<0)-scs(scs<0))
mean(ccs(scs<0)-scs(scs<0))
[h,p]=ttest(ccs(scs>0)-scs(scs>0))
mean(ccs(scs>0)-scs(scs>0))




ndc=100;
dc=[0:ndc]/ndc*min(Ba);
for j=1:length(dc)
    [h,pdc(j)]=ttest(B(Ba<=dc(j)));
end
figure(9)
plot(dc,pdc)
mean(B(Ba<0))





