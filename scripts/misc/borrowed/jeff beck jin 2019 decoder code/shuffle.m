function [X] = shuffle(X)
   [N,D]=size(X);
   for i=1:D
       X(:,i)=X(randperm(N),i);
   end
end

