function[DV,DVAll,B,Ball]=getDVs(data,dataAll,issparse,prefs,kappa,f)

% expects data{j}.X = neural activity, data{j}.Y = condition label, data{j}.Xraw is
% raw traces in each of the j=1:3 conditions:  250,750,control

% outputs struct DV which is the decision variable for each decoder and 
% each training test condition, i.e. DV(k,j) is train k test j.  k,k case
% is leave one out.
% struct contains
% DV.opt, DV.opt8, DV.est, DV.logit, DV.sum ... 

NR=size(data{1}.X,2);

% % Build dataAll
% dataAll.X=[];
% dataAll.Y=[];
% dataAll.cond=[];
% for j=1:3
%     dataAll.X=[dataAll.X;data{j}.X];
%     dataAll.Xraw=[dataAll.X;data{j}.Xraw];        
%     dataAll.Y=[dataAll.Y;data{j}.Y];
%     dataAll.cond=[dataAll.cond;j*ones(size(data{j}.Y));];
% end

%do cross validations for self-model comparisons

%Cross validate and perform cross validation cross condition checks
maxiters=500;
ntheta=500;
theta = [1:ntheta]/ntheta*2*pi-pi;
thetaS = [1:8]/8*2*pi-pi;
dtheta=theta(2)-theta(1);
issparse=0;

for k=1:length(data)
for j=1:length(data)
    DV{k,j}.opt=[];
    DV{k,j}.opt8=[];
    DV{k,j}.est=[];
    DV{k,j}.JM=[];
    DV{k,j}.PV=[];
    DV{k,j}.PVemp=[];
    DV{k,j}.logit=[];
    DV{k,j}.sum=[];
    DV{k,j}.max=[];
    DV{k,j}.pr=[];
    DV{k,j}.PVempPrefs=[];
end
end

for j=1:length(data)
    for n=1:size(data{j}.X,1)
        m1 = VonMesisRegression(NR,issparse);
        idx = [1:n-1,n+1:size(data{j}.X,1)];
        m1.fit(data{j}.Y(idx,1)*2*pi/8,data{j}.X(idx,:),maxiters);
        [postmu{j,j}(n,1),postkappa{j,j}(n,1)]=m1.getPredictions(data{j}.X(n,:));
        pr{j,j}(n,:)=m1.getPdf(theta,data{j}.X(n,:));
        DV{j,j}.opt(n,1) = 1-sum(pr{j,j}(n,abs(theta)<pi/8),2)/sum(pr{j,j}(n,:),2);
        DV{j,j}.opt8(n,1) = 1-m1.getPdf(0,data{j}.X(n,:))/sum(m1.getPdf(thetaS,data{j}.X(n,:)));
        DV{j,j}.est(n,1) = (postmu{j,j}(n,1))/pi;

        prefs_emp=PV_emp(data{j}.Y(idx,1)*2*pi/8,data{j}.Xraw(idx,:));
        DV{j,j}.PVemp(n,1) = (angle(data{j}.Xraw(n,:)*cos(prefs_emp)+sqrt(-1)*data{j}.Xraw(n,:)*sin(prefs_emp)))/pi;

%        B1 = mnrfit(data{j}.X(idx,:),(data{j}.Y(idx,1)~=8)+1);
%        temp=mnrval(B1,data{j}.X(n,:));
%        logitDV{j,j}(n,1) = temp(2);
%        temp=mnrval(B1,data{j}.X(n,:));
%        logitDV{j,j}(n,1) = temp(2);
        B1 = glmfit(data{j}.X(idx,:),(data{j}.Y(idx,1)~=8),'binomial','link','logit','Constant','off');
        DV{j,j}.logit(n,1) = glmval(B1,data{j}.X(n,:),'logit','Constant','off');
    end
    DV{j,j}.pr=pr{j,j};
    model{j} = VonMesisRegression(NR,issparse);
    model{j}.fit(data{j}.Y*2*pi/8,data{j}.X,maxiters)
    [B{j},dev,stats] = glmfit(data{j}.X,(data{j}.Y~=8),'binomial','link','logit','Constant','off');
    DV{j,j}.PVempPrefs = PV_emp(data{j}.Y*2*pi/8,data{j}.Xraw);
end


DVsummax=-Inf(3,3);
DVsummin=Inf(3,3);
DVmaxmax=-Inf(3,3);
DVmaxmin=Inf(3,3);

for j=1:length(data)
for k=1:length(data)
    if(j~=k)
        pr{k,j}=model{k}.getPdf(theta,data{j}.X);
        DV{k,j}.opt = 1-sum(pr{k,j}(:,abs(theta)<2*pi/16),2)./sum(pr{k,j},2);        
        DV{k,j}.opt8 = 1-model{k}.getPdf(0,data{j}.X)./sum(model{k}.getPdf(thetaS,data{j}.X),2);
        [postmu{k,j},postkappa{k,j}]=model{k}.getPredictions(data{j}.X);
        DV{k,j}.est = (postmu{k,j})/pi;
        DV{k,j}.PV = [];
        DV{k,j}.logit = glmval(B{k},data{j}.X,'logit','Constant','off');
        DV{k,j}.pr=pr{k,j};
        DV{k,j}.PVemp = (angle(data{j}.Xraw*cos(DV{k,k}.PVempPrefs)+sqrt(-1)*data{j}.Xraw*sin(DV{k,k}.PVempPrefs)))/pi;    
        DV{k,j}.PVempPrefs = DV{k,k}.PVempPrefs;
    end
    DV{k,j}.sum = sum(data{j}.Xraw,2);
    DV{k,j}.max = max(data{j}.Xraw')';
%    DV{k,j}.sum=DV{k,j}.sum-min(DV{k,j}.sum);
    DVsummax(k,j)=max(DVsummax(k,j),max(DV{k,j}.sum));
    DVsummin(k,j)=min(DVsummin(k,j),min(DV{k,j}.sum));
    DVmaxmax(k,j)=max(DVmaxmax(k,j),max(DV{k,j}.max));
    DVmaxmin(k,j)=min(DVmaxmin(k,j),min(DV{k,j}.max));
    DV{k,j}.PV = angle(data{j}.Xraw*(cos(prefs').*kappa')+sqrt(-1)*data{j}.Xraw*(sin(prefs').*kappa'))/pi;    
    
    
    [m,loc] = max((data{j}.Xraw*log(f')-sum(f'))');
    
    DV{k,j}.JM = mod(loc'/180*2*pi+pi,2*pi)/pi-1;
    
end
end

for k=1:length(data)
for j=1:length(data)
    DV{k,j}.sum=(DV{k,j}.sum-min(DVsummin(k,:)))/(max(DVsummax(k,:))-min(DVsummin(k,:)));
    DV{k,j}.max=(DV{k,j}.max-min(DVmaxmin(k,:)))/(max(DVmaxmax(k,:))-min(DVmaxmin(k,:)));
end
end

prefs_empALL=PV_emp(dataAll.Y*2*pi/8,dataAll.Xraw);
for j=1:length(data)
    prALL{j}=[];
    DVAll{j}.opt=[];
    DVAll{j}.opt8=[];
    DVAll{j}.est=[];
    DVAll{j}.JM=[];
    DVAll{j}.PV=[];
    DVAll{j}.PVemp=[];
    DVall{j}.PVempPrefs=[];
    DVAll{j}.logit=[];
    DVAll{j}.sum=[];
    DVAll{j}.max=[];
    DVAll{j}.pr=[];
end

for n=1:size(dataAll.X,1)
    [Ball,dev,stats] = glmfit(dataAll.X(1:n-1:n+1:end,:),(dataAll.Y(1:n-1:n+1:end,1)~=8),'binomial','link','logit','Constant','off');
    mAll = VonMesisRegression(NR,issparse);
    mAll.fit(dataAll.Y(1:n-1:n+1:end,1)*2*pi/8,dataAll.X(1:n-1:n+1:end,:),maxiters)
    
    prAll=mAll.getPdf(theta,dataAll.X(n,:));
    for j=1:length(data)
        if(dataAll.cond(n)==j)
            DVAll{j}.opt=[DVAll{j}.opt;1-sum(prAll(abs(theta)<2*pi/16),2)/sum(prAll);];
            DVAll{j}.opt8=[DVAll{j}.opt8;1-mAll.getPdf(0,dataAll.X(n,:))/sum(mAll.getPdf(thetaS,dataAll.X(n,:)));];
            [postmuAll,postkappaAll]=mAll.getPredictions(dataAll.X(n,:));
            DVAll{j}.est=[DVAll{j}.est;abs(postmuAll)/pi;];
            
            DVAll{j}.logit = [DVAll{j}.logit;glmval(Ball,dataAll.X(n,:),'logit','Constant','off');];
            DVAll{j}.pr= [DVAll{j}.pr;prAll;];
        end
    end
end


for j=1:length(data)
        
    DVAll{j}.JM = DV{j,j}.JM;
    DVAll{j}.PV = DV{j,j}.PV;
    
    DVAll{j}.PVemp = (angle(data{j}.Xraw*cos(prefs_empALL)+sqrt(-1)*data{j}.Xraw*sin(prefs_empALL)))/pi;    
    DVAll{j}.PVempPrefs = prefs_empALL;
    
    DVAll{j}.sum=DV{j,j}.sum;
    DVAll{j}.max=DV{j,j}.max;
end

    
end