function [bin_pca,rest]=designbinary_pca(num,Sigma,d)
%
% designbinary_pca.m
%
% This function calculates a matriz containing
% the binary principal components of a covariance matrix
% Sigma
%
% Copyright (2018): Jonathan Monsalve
% Universidad Industrial de Santander (2018)
% -----------------------------------------------------------------------
%
% designbinary_pca is distributed under the terms
% of the GNU General Public License 2.0.
%
% Permission to use, copy, modify, and distribute this software for
% any purpose without fee is hereby granted, provided that this entire
% notice is included in all copies of any software which is or includes
% a copy or modification of this software and in all copies of the
% supporting documentation for such software.
% This software is being provided "as is", without any express or
% implied warranty.  In particular, the authors do not make any
% representation or warranty of any kind concerning the merchantability
% of this software or its fitness for any particular purpose."
% ----------------------------------------------------------------------
%
%
%   ===== Required inputs =====
%
%   num:        Number of binary vectors to be designed
%   Sigma:      Covariance matrix
%   d:          Transmitance of the output filters (example d=1/4 for 25%)
%
% 	===== Outputs =====
%   bin_pca:    Matrix where each row is a binary PC.
%   rest:       Represent what was deflated from the covariance maatrix.
%   Such that Sigma-rest is th residual
% ========================================================
[L,~]=size(Sigma);

bin_pca=zeros(num,L);
prior_list=zeros(1,L);
to=0;
rest=0;
tr=trace(Sigma);
list=1:L;
for k=1:num
    vm=-inf;
    index=zeros(L,1);
    mfun=0;
    for i=1:L %forward propagation
        for j=list
            v=Sigma(j,j);
            v1=sum(Sigma(index(1:i-1),j),1);
            v=v+(2*v1);
            v=(mfun+v)/i;
            if v>vm
                index(i)=j;
                vm=v;
            end
        end
        if index(i)==0
            break
        end
        if(sum(index>0)>=round(L*d))
           break; 
        end
        mfun=(vm)*i;
        list=list(list~=index(i));
    end
    q=vm;
    index=index(index>0)';
    mfun=vm*length(index);
    for i=index % backward propagation
       v=Sigma(i,i);
       v1=sum(Sigma(index(index~=i),i),1);
       v=v+(2*v1);
       q1=(mfun-v)/(length(index)-1);
       if(q1>q)
          q=q1;
          mfun=(mfun-v);
          index=index(index~=i);
       end
    end
    index1=index(index>0);
    index1=index1(index>0);
    prior_list(index1)=prior_list(index1)+1;
    [~,list]=sort(prior_list);
    s=zeros(L,1);
    s(index(index>0))=1;
    bin_pca(k,:)=s';
    s=s./norm(s);
    temp=(s'*Sigma*s);
    to=to+temp;
    Q=s*s';
    rest=rest+((Sigma*Q)+(Q*Sigma)-(Q*Sigma*Q));
    
    Sigma = Sigma-(Sigma*Q)-(Q*Sigma)+(Q*Sigma*Q);
    if(temp/tr<1e-6) % if no enough variance is explaind finish
       bin_pca=bin_pca(1:k-1,:);
       break; 
    end
end

