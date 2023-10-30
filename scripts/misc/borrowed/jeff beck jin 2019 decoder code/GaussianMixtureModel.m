classdef GaussianMixtureModel < handle

    properties
        % Parameters
        Nc % number of clusters
        D  % dimension of the observations
        L % lower bound
        alpha_0
        minclustersize
        
        % Cluster Prior Parameters
        mu_0
        lambda_0
        V_0
        nu_0
        
        p  % assignment probabilities (N x Nc)
        logptilde
        NA % number of assignments to each class
        

        pi 
        NWs  % cluster parameters
    end
    
    methods
        function self = GaussianMixtureModel(Nc,D,alpha_0,mu_0,lambda_0,V_0,nu_0)
             self.minclustersize = 0.1;
             if (nargin == 3)
                 self.Nc = Nc;
                 self.D = D;
                 self.alpha_0 = alpha_0/Nc;
                 self.mu_0 = zeros(D,1);
                 self.lambda_0=1;
                 self.nu_0 = D+2;
                 self.V_0 = eye(D)/(D+2)*Nc;
             elseif(nargin==7)
                 self.Nc = Nc;
                 self.D = D;
                 self.alpha_0 = alpha_0/Nc;
                 self.mu_0 = mu_0;
                 self.lambda_0 = lambda_0;
                 self.nu_0 = nu_0;
                 self.V_0 = V_0;
             else
                 'invalid inputs'
                 stop
             end
                    
             self.L = -Inf;
             self.pi=dists.expfam.dirichlet(self.Nc,self.alpha_0*ones(1,self.Nc));
             for i=1:Nc
                 self.NWs{i}=dists.expfam.NW(self.mu_0,self.lambda_0,self.V_0,self.nu_0);
             end
             
        end
        
        function fit(self,X,tol,maxiters,Ncguess)
            if(~exist('Ncguess','var'))
                Ncguess=self.Nc;
            elseif(Ncguess>self.Nc)
                for i=self.Nc+1:Ncguess
                    self.NWs{i}=dists.expfam.NW(self.mu_0,self.lambda_0,self.V_0,self.nu_0);
                end
                self.alpha_0 = self.alpha_0*self.Nc/Ncguess;
                self.Nc = Ncguess;
                self.pi=dists.expfam.dirichlet(self.Nc,self.alpha_0*ones(1,self.Nc));
            elseif(Ncguess<self.Nc/2)
                self.Nc = ceil(self.Nc/2);
                self.pi=dists.expfam.dirichlet(self.Nc,self.alpha_0*ones(1,self.Nc));
            end
                        
            tic
            k=0;
            Llast=-Inf;
            self.smartinitialization(X,Ncguess);
            stop=0;
            while(stop<=0 & k < maxiters)
                k=k+4;
                Llast = self.L;
                self.update(X,0);
                self.update(X,0);
                self.update(X,0);
                self.update(X,0);
%                self.merge(X);
                if(self.L-Llast < abs(self.L*tol))
%                     self.merge(X);
%                     self.merge(X);
%                     self.merge(X);
%                     self.merge(X);
%                     self.merge(X);
                    if(self.L-Llast < abs(self.L*tol) )
                        stop=stop+1; 
                    end
                end
%                 figure(2)
%                 hold on
%                 plot(k,self.L,'k.')
%                 hold off
%                 if(mod(k,5)==1)
%                    self.plotclusters(X,1);
%                    title(strcat('ELBO = ',num2str(self.L)))
%                    drawnow
%                    pause
% %                    self.perturbunusedclusters;
%                 end
                 
            end
            self.plotclusters(X,1)
            if (k>=maxiters)
                fprintf('maximum iterations reached\n')
            else
                fprintf(['Discovered ',num2str(sum(self.NA>1)),' clusters after ',num2str(k),' iterations in ',num2str(toc),' seconds\n'])
            end,
            fprintf(['Final <ELBO> = ',num2str(self.L),'\n'])
        end
        
        function smartinitialization(self,X,Ncguess)
            CC = cov(X);
            C=sqrtm(CC);
            mu=mean(X)';
            [idx, muk] = kmeans(X, self.Nc);

            for i=1:self.Nc
                self.NWs{i}.mu_0 = mu;
                
                self.NWs{i}.invV_0 = CC*self.NWs{i}.nu_0;%-self.D-1);
%                 self.NWs{i}.invV_0 = cov(X(idx2,:))*(self.NWs{i}.nu_0-self.D-1);
                
                self.NWs{i}.logdetinvV_0 = self.NWs{i}.logdet(self.NWs{i}.invV_0);                
                self.NWs{i}.V_0 = inv(self.NWs{i}.invV_0);

                self.NWs{i}.nu = self.NWs{i}.nu_0 + size(X,1)/Ncguess;
                self.NWs{i}.lambda = self.NWs{i}.lambda_0 + size(X,1)/Ncguess;

                self.NWs{i}.mu = muk(i,:)'; %mu + C*randn(size(mu));
                
                idx2=find(idx==i);
                if(length(idx2)>2)
                    self.NWs{i}.invV = self.NWs{i}.invV_0 + cov(X(find(idx==i),:))*(self.NWs{i}.nu_0-self.D-1);
                else
                    self.NWs{i}.invV = self.NWs{i}.invV_0 + C*cov(randn(size(mu')))*C'*(self.NWs{i}.nu_0-self.D-1);                   
                end
                
                self.NWs{i}.V = inv(self.NWs{i}.invV);
                self.NWs{i}.logdetinvV = self.NWs{i}.logdet(self.NWs{i}.invV);

                self.NWs{i}.setUpdated(false);
            end
            self.logptilde = zeros(size(X,1),self.Nc);
            self.updateassignments(X);
        end
        
        function DL = update(self,X,iters,fighandle)
            if(~exist('fighandle','var'))
                fighandle=0;
            end
                
            for i=1:iters
                L=self.L;
                self.updateparms(X,fighandle);
                self.updateassignments(X); 
                DL = self.L - L;
            end
        end
        
        function updateassignments(self,X)

            [N,D]=size(X);
            self.logptilde=zeros(N,self.Nc);
            for i=1:self.Nc
                self.logptilde(:,i) = self.NWs{i}.Eloglikelihood(X);
            end
            self.logptilde = bsxfun(@plus,self.logptilde,self.pi.loggeomean);

            self.p = bsxfun(@minus,self.logptilde,max(self.logptilde')');            
            self.p = exp(self.p);
            self.p = bsxfun(@rdivide,self.p, sum(self.p,2));

            self.NA = sum(self.p,1);
            
            self.L = - self.pi.KLqprior;            
            for i=1:self.Nc
%                if(self.NA(i)>0)                        
                    self.L = self.L - self.NWs{i}.KLqprior;
%                end
            end
            
            self.L = self.L + sum(sum(self.p.*(self.logptilde)));
            idx = find(self.p(:)>0);
            self.L = self.L - sum(self.p(idx).*log(self.p(idx)));

        end
        
        function KLqprior(self)            

            self.L = - self.pi.KLqprior;            
            for i=1:self.Nc
                self.L = self.L - self.NWs{i}.KLqprior;
            end
                        
        end
        
        function updateparms(self,X,fighandle)
            if(isempty(self.p))
                self.updateassignments(X);
            end
            self.pi.update(self.NA);
            
            for i=1:self.Nc
                if(self.NA(i)>self.minclustersize)
                    Ex = (self.p(:,i)'*X)/self.NA(i);
                    Exx = bsxfun(@times,X,sqrt(self.p(:,i)));
                    Exx = Exx'*Exx/self.NA(i);
                    self.NWs{i}.update(Ex',Exx,self.NA(i));
                else
                    self.NWs{i}.update(0,0,0);
                end
            end
            if(fighandle>0)
                self.plotclusters(X,fighandle);
            end
        end
        
        function perturbunusedclusters(self)
            idx = find(self.NA<1);
            for i=1:length(idx)
                self.NWs{idx(i)}.mu =  randn(self.D,1);
            end
        end
        
        function merge(self,X)
            idx=find(self.NA>1);
            if(length(idx)<2) %do nothing
                fprintf('no possible merges\n')
            else                
                idx=idx(randperm(length(idx)));
                i=idx(1);
                j=idx(2);
            
                psave = self.p;
                NAsave = self.NA;
                Lsave = self.L;
                self.p(:,i) = (self.p(:,j)+self.p(:,i));
                self.p(:,j) = 0;
                self.NA(i)=self.NA(i)+self.NA(j);
                self.NA(j)=0;

                self.updateparms(X,0);
                self.updateassignments(X);

                if(self.L <= Lsave) % reject merge
                    self.p = psave;
                    self.NA = NAsave;
                    self.updateparms(X,0);
                    self.L = Lsave;
                end
            end
        end
        
        function split(self,X)
            self.update(X,1);
            [m,j]=min(self.NA);
            i=util.discretesample(self.pi.mean,1);
            [m,idx]=max(self.p');
            idxi=find(idx==i);
            self.p(idxi,:)=zeros(length(idxi),self.Nc);
            [z,c]=kmeans(X(idxi,:),2);
            
            
            self.p(idx==i,i)=real(z==1);
            self.p(idx==i,j)=real(z==2);
            
            self.updateparms(X,0);
        end
        
        function [pc,prYz]=plotclusters(self,X,fighandle,Y,pcutoff,V,D,mu)
            cc = hsv(self.Nc);
            [temp,idxc] = max(self.p');
            
            shape='xo+*sdv^ph<>.';
            if(~exist('Y','var'))
               Y=ones(size(X,1),1); 
               label=1;
               pc=NaN;
               prYz=NaN;
            elseif(isempty(Y))
               Y=ones(size(X,1),1); 
               label=1;
               pc=NaN;
               prYz=NaN;
            else
                %compute performance of MAP decoder
                [m,z]=max(self.p');                 
                zlabel=unique(z);
                label = unique(Y);
                if(~exist('pcutoff','var'))
                    pcutoff=1/length(zlabel);
                end
                pc=0;
                for i=1:length(zlabel)
                    prz(i)=mean(z==zlabel(i));
                    for j=1:length(label)
                        prYz(j,i)=mean(z'==zlabel(i) & Y==label(j));
                    end
                    [m,mapgz]=max(prYz(:,i));
                    pc = pc + prYz(mapgz,i);                     
                end
                if(length(zlabel)==2)
                    Yhat = self.p(:,zlabel(1))>pcutoff;
                    sw=0;
                    if(corr(Yhat,Y)<0)
                        Yhat=1-Yhat;
                        sw=1;
                    end
                    pc = mean(Yhat==Y);
                    if(sw==0)
                        idxc(Yhat==0)=zlabel(2);
                        idxc(Yhat==1)=zlabel(1);
                    else
                        idxc(Yhat==0)=zlabel(1);
                        idxc(Yhat==1)=zlabel(2);                        
                    end
                end
            end
            if(~exist('mu','var'))
                D=ones(size(X,2),1);
                V=eye(size(X,2));
                mu=zeros(size(X,2),1)';
            end
            
            X=X*diag(sqrt(D))*V;
            X=bsxfun(@plus,X,mu);
            mu=mu';
            for i=1:self.D
            for j=i+1:self.D    
                figure(fighandle);
                fighandle = fighandle+1;
                for n=1:length(label)
                    idx=find(Y==label(n));
                    scatter(X(idx,i),X(idx,j),50*ones(size(idx)),cc(idxc(idx),:),shape(n))
                    hold on
                end
                t=[0:1:100]/100*2*pi;
                for k=1:self.Nc
                    if(self.NA(k)>= 1)
                        C =V'*diag(sqrt(D))*(self.NWs{k}.ESigma)*diag(sqrt(D))*V;
                        C=C([i,j],[i,j]);
                        C=sqrtm(C);                
                        nwsmu = V'*diag(sqrt(D))*self.NWs{k}.mu + mu;
                        stdring = repmat(nwsmu([i,j]),1,length(t)) + 2*C*[sin(t);cos(t)];
                        plot(stdring(1,:),stdring(2,:),'color',cc(k,:))
                    end
                end 
                hold off
                title(['Observable Dimensions = ',num2str(i),' and ',num2str(j)])
                xlabel(['Dimension ',num2str(i)])
                ylabel(['Dimension ',num2str(j)])
            end
            end
            
        end
        
    end
    
end

