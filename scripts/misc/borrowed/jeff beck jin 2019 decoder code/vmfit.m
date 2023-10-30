function [Yhat, YhatTE, MEtrain, MSEtest, A, B] = vmfit( X,Y ,XTE, YTE)

    % Takes in training data X and Y and test data XTE and YTE
    % Fits von mises regression and generates predicstions for tetsting 
    % data.  Here Y is assumed to be in radians (i.e. period 2*pi).
    % Additional outputs are mean squared errors and the weights of the 
    % von-mesis regression which are composed of two vectors of weights:
    % A and B.  

     [ns,NR]=size(X);
     XTE=bsxfun(@plus,XTE,-mean(X));
     X=bsxfun(@plus,X,-mean(X));
     invD=diag(1./std(X));
     X=X*invD;
     XTE=XTE*invD;
     X=[X,ones(ns,1)];
     XTE=[XTE,ones(size(XTE,1),1)];
     
     NR=NR+1;
     model=VonMesisRegression(NR,0);
     model.fit(Y,X,50,0);

     A=model.a.mean;
     B=model.b.mean;
     A=A(1:end-1);
     B=B(1:end-1);
     A=invD*A;
     B=invD*B;
     
    [Yhat,kappa] = model.getPredictions(X);              
    [YhatTE,kappa] = model.getPredictions(XTE);         
    
    MSEtrain = mean((mod(Y-Yhat+pi,2*pi)-pi).^2);
    MSEtest = mean((mod(YTE-YhatTE+pi,2*pi)-pi).^2)
end

