function [bin_pca]=designbinary_pca(num,L,Sigma)
%  
% designbinary_pca.m
%
% This function calculates a matriz containing
% the binary principal components of a covariance matrix
% Sigma
% 
% Copyright (2012): Jonathan Monsalve
% Universidad Industrial de Santander (2017)
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
%   L:          Number of spectral bands (dimension of Sigma)     
%   Sigma:      Covariance matrix
%
% 	===== Outputs =====
%   bin_pca:    Matrix where each row is a binary PC.
% ========================================================
bin_pca=zeros(num,L);
for k=1:num
    s=zeros(L,1);
    vm=0;
    index=zeros(L/1,1);
    mfun=0;
    for i=1:(L/1)
        for j=1:L
            if(sum(index==j)==0)
                 v=Sigma(j,j);
                 v1=sum(Sigma(index(1:i-1),j),1);
                 v=v+(2*v1);
                 v=(mfun+v)/i;
                if abs(v)>abs(vm)
                    index(i)=j;
                    vm=v;
                end
            end
        end
        
        if index(i)==0
           break 
        end
        mfun=(vm)*i;
        
    end
    s(index(index>0))=1;
    bin_pca(k,:)=s';
    s=s./norm(s);
    Sigma=Sigma-(s*s'*(s'*Sigma*s));
end
