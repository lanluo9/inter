classdef LogisticRegression < handle
    properties
        
        NR %number of regressors
        W
        xi
        SEXX
        SEYX
        SEYX
        n
        
        L
        iters
    end
    
    methods
        function self = LogisticRegression(NR,issparse)
            if(~exist('issparse','var'))
                issparse = 0;
            end
            if(issparse==0)
                self.W = dists.expfam.MVN(zeros(NR,1),eye(NR)); %second arg is precision of prior
            else
                self.W = dists.normalsparse(NR,1/NR,1/NR);
            end
            self.iters=0;
            self.L=-Inf;
        end
        
        function fit(self,Y,X,maxiters)
            DL=Inf;
            while((DL/abs(self.L)>1e-12 | self.iters<5) & self.iters < maxiters)
                DL = self.update(Y,X,1);
            end
            if(self.iters == maxiters) 
                fprintf('maxiters reached\n')
            end
        end
        
        function DL = update(self,Y,X,niters)
            if(~exist('niters','var'))
                niters = 1;
            end
            for i=1:niters
                DL = self.update_suffstats(Y,X);
                self.updateparms;
                self.iters=self.iters+1;
                if(DL<0)
                fprintf(['Warning Lower Bound Decreasing (DL/L = ',num2str(DL/abs(self.L)),') at iteration ',num2str(self.iters),'\n'])
                end
            end
        end
        
        function DL = update_suffstats(self,Y,X,p)
            
            DL=self.L;
            Ns=size(Y,1);
            if(~exist('p','var'))
                p=ones(Ns,1);
                self.n = Ns;
            else
                self.n=sum(p);
            end
            
            self.xi = sqrt(sum(X*(self.W.secondmoment).*X,2));
            self.lambda = 0.5./self.xi.*(1./(1+exp(-self.xi))-0.5);
            self.SEXX = X'*bsxfun(@times,X,self.lambda(self.xi).*p);
            self.SEYX = ((Y.*p)'*X);

            
            self.L = - self.KLqprior + 0.5*(Y'*X*W) ...
                     - sum(X*(self.W.secondmoment).*bsxfun(@timees,X,self.lambda),2) ...
                     - sum(log(1+exp(-xi))) - sum(self.xi)/2 + sum(self.lambda.*self.xi.^2);
            DL=self.L-DL;

        end

        function updateparms(self)
            self.W.updateSS(self.SEYX'/self.n,self.SEXX/self.n,self.n);
        end
        
                
        function res = getPYeq1(self,theta,X)  %%%NEEDS WORK
            
            self.xi = sqrt(sum(X*(self.W.secondmoment).*X,2));
            self.lambda = 0.5./self.xi.*(1./(1+exp(-self.xi))-0.5);
            
        end
 
        function KL = KLqprior(self)
            KL = self.W.KLqprior;
       end
                
    end
   
end
