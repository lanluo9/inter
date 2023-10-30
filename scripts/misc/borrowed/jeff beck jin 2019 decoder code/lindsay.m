
data.X=[];
data.Y=[];

for n=1:8
    data.X = [data.X;testResp{n}';];
    data.Y = [data.Y;n*ones(size(testResp{n},2),1);];
end
NS=size(data.Y,1);
data.X=data.X(:,ind);
mu=mean(data.X);
stdX=std(data.X);
data.X=bsxfun(@plus,data.X,-mu);
data.X=bsxfun(@times,data.X,1./stdX);

data.X=[data.X,ones(size(data.X,1),1)];
NR = size(data.X,2);

X=[data.Y/8*2*pi,data.X];
f=@(x) fcn(x,X);
A=zeros(2*NR,1);
Opts = optimset('fminsearch');
Opts.MaxFunEvals=10000*NR*2;
Opts.MaxIter=10000*NR*2;

for n=1:10
    A=fminsearch(f,A,Opts);
end
B=A(NR+1:end,1);
A=A(1:NR,1);

ABnorm = sqrt(A.^2+B.^2);
[m,loc]=sort(ABnorm,'descend');
j=find(loc==NR);
loc(2:j)=loc(1:j-1);
loc(1)=NR;

i=sqrt(-1);
T=1000;
theta=[-T:T]/T*pi/8;
dtheta = theta(2)-theta(1);
NS=size(data.Y);
a=[];
b=[];
for n=1:NR
    
    idx=loc(1:n);
    a=[a;0;];
    b=[b;0;];
    for j=1:NS    
%    for j=1:20
%        idx=randperm(NR);

        X=data.X([1:j-1,j+1:NS],idx);
        Y=data.Y([1:j-1,j+1:NS]);
        X=[Y/8*2*pi,X];

        f=@(x) fcn(x,X);

        a=fminsearch(f,[a;b;],Opts);
        b=a(n+1:end,1);
        a=a(1:n,1);
                
        thetahat = angle(data.X(:,idx)*a+i*data.X(:,idx)*b);
        kappa = abs(data.X(:,idx)*a+i*data.X(:,idx)*b);
        
        thetahatsave(n,j,1:NS) = thetahat(1:NS);
        kappasave(n,j,1:NS) = kappa(1:NS);
        
        p = bsxfun(@times,cos(ones(length(kappa),1)*theta-thetahat*ones(1,length(theta))),kappa);
        p = sum(exp(p),2)/2/pi./besseli(0,kappa)*dtheta;
        thetahat = thetahat/pi*180;
        for k=1:8
            bias(n,j,k) = mean(mod(thetahat(data.Y==k)-k/8*360+180,360))-180;
            sigmahat(n,j,k) = std(mod(thetahat(data.Y==k)-k/8*360+180,360)-180);
            pbar{n,j,k} = p(data.Y==k);
        end
        kappabar(n,j)=mean(kappa);
        sigmabar(n,j)=mean(sigmahat(n,j,:));
        psave(n,j,1:NS)=p(1:NS);        
       thetacv(n,j)=thetahatsave(n,j,j);
       kappacv(n,j)=kappasave(n,j,j);
       Ycv(n,j)=data.Y(j);
    end
    for k=1:8
        biascv(n,k)=mean(mod(thetacv(n,Ycv(n,:)==k)/pi*180-k/8*360+180,360)-180);
        sigmacv(n,k)=std(mod(thetacv(n,Ycv(n,:)==k)/pi*180-k/8*360+180,360)-180);
    end
    save lindseysave
    n
end
figure(1)
plot(2:n,mean(sigmabar(2:n,:)')/2)
hold on
plot(2:n,mean(sigmacv(2:n,:)')/2)
hold off
legend('In sample','Out of Sample')
ylabel('Discrimination Threshold')
xlabel('Number of nurons')

for j=1:NS
    psavecv(:,j)=squeeze(psave(:,j,j));
end

nb=200;
figure(2)
for n=10:NR
    for j=1:nb
        d(j)=(j-1)/(nb-1);
        cd(n,j)=mean(psavecv(n,data.Y==8)>d(j));
        fp(n,j)=mean(psavecv(n,data.Y<8)>d(j));
    end
    plot(fp(n,:),cd(n,:))
    hold on
end
hold off
xlabel('false positives')
ylabel('correct detections')
title('Cross-validated ROC curves when using the best 10 to 20 units')
figu
figure(3)
start=10;
stop=40;
plot(mean(fp(start:stop,:)),mean(cd(start:stop,:)))
hold on
plot(mean(fp(start:stop,:))+std(fp(start:stop,:)),mean(cd(start:stop,:))+std(cd(start:stop,:)),'k')
plot(mean(fp(start:stop,:))-std(fp(start:stop,:)),mean(cd(start:stop,:))-std(cd(start:stop,:)),'k')
hold off
xlabel('false positives')
ylabel('correct detections')
axis([0 1 0 1])
title('Average +- standard error of the ROC for 10 to 20 units')


%%% build ROC for adapted cases...


data.X=[];
data.Y=[];

for n=1:8
    data.X = [data.X;testResp{n}';];
    data.Y = [data.Y;n*ones(size(testResp{n},2),1);];
end
NS=size(data.Y,1);
data.X=data.X(:,ind);
mu=mean(data.X);
stdX=std(data.X);
data.X=bsxfun(@plus,data.X,-mu);
data.X=bsxfun(@times,data.X,1./stdX);

data.X=[data.X,ones(size(data.X,1),1)];
NR = size(data.X,2);

X=[data.Y/8*2*pi,data.X];
f=@(x) fcn(x,X);
A=zeros(2*NR,1);
Opts = optimset('fminsearch');
Opts.MaxFunEvals=10000*NR*2;
Opts.MaxIter=10000*NR*2;

for n=1:10
    A=fminsearch(f,A,Opts);
end
B=A(NR+1:end,1);
A=A(1:NR,1);


ABnorm = sqrt(A.^2+B.^2);
[m,loc]=sort(ABnorm,'descend');
j=find(loc==NR);
loc(2:j)=loc(1:j-1);
loc(1)=NR;

i=sqrt(-1);
T=1000;
theta=[-T:T]/T*pi/8;
dtheta = theta(2)-theta(1);
NS=size(data.Y);


n=21;
X=[data.Y/8*2*pi,data.X(:,loc(1:n))];
f=@(x) fcn(x,X);
a=zeros(2*n,1);
for n=1:10
    a=fminsearch(f,a,Opts);
end
b=a(n+1:end,1);
a=a(1:n,1);


idx=loc(1:n);
thetahat = angle(data.X(:,idx)*a(1:length(idx),1)+i*data.X(:,idx)*b(1:length(idx),1));
kappa = abs(data.X(:,idx)*a(1:length(idx),1)+i*data.X(:,idx)*b(1:length(idx),1));

p = bsxfun(@times,cos(ones(length(kappa),1)*theta-thetahat*ones(1,length(theta))),kappa);
p = sum(exp(p),2)/2/pi./besseli(0,kappa)*dtheta;
thetahat = thetahat/pi*180;
for k=1:8
    bias(k) = mean(mod(thetahat(data.Y==k)-k/8*360+180,360))-180;
    sigmahat(k) = std(mod(thetahat(data.Y==k)-k/8*360+180,360)-180);
end
clear fp cd d
nb=200;
for j=1:nb
        d(j)=(j-1)/(nb-1);
        cd(j)=mean(p(data.Y==8)>d(j));
        fp(j)=mean(p(data.Y<8)>d(j));
end
plot(fp,cd,'b')
hold on

%%%%%%%%%%%%THE ABOVE SHOULD BE CROSSVALIDATED TO MAKE COMPARISON FAIR>>>>>



data2.X=[];
data2.Y=[];

testc=1;
for n=1:8
    data2.X = [data2.X;ppResp{testc,n}';];
    data2.Y = [data2.Y;n*ones(size(ppResp{testc,n},2),1);];
end
NS=size(data2.Y,1);
data2.X=data2.X(:,ind);
data2.X=bsxfun(@plus,data2.X,-mu);
data2.X=bsxfun(@times,data2.X,1./stdX);
data2.X=[data2.X,ones(size(data2.X,1),1)];


idx=loc(1:n);
thetahat2 = angle(data2.X(:,idx)*a(1:length(idx),1)+i*data2.X(:,idx)*b(1:length(idx),1));
kappa2 = abs(data2.X(:,idx)*a(1:length(idx),1)+i*data2.X(:,idx)*b(1:length(idx),1));

p2 = bsxfun(@times,cos(ones(length(kappa2),1)*theta-thetahat2*ones(1,length(theta))),kappa2);
p2 = sum(exp(p2),2)/2/pi./besseli(0,kappa2)*dtheta;
thetahat2 = thetahat2/pi*180;
for k=1:8
    bias2(k) = mean(mod(thetahat2(data2.Y==k)-k/8*360+180,360))-180;
    sigmahat2(k) = std(mod(thetahat2(data2.Y==k)-k/8*360+180,360)-180);
end
clear fp cd d
nb=200;
for j=1:nb
        d(j)=(j-1)/(nb-1);
        cd(j)=mean(p2(data2.Y==8)>d(j));
        fp(j)=mean(p2(data2.Y<8)>d(j));
end
plot(fp,cd,'r')
xlabel('false positives')
ylabel('correct detections')
title('ROC curves')
hold off
legend('control','adapted')
