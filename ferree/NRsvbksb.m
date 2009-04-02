function x=NRsvbksb(u,w,v,m,n,b)

% function x=svbksb(u,w,v,m,n,b)
% Does back-substitutation for matrix inversion via SVD.
% First implementation follows exactly that in Numerical Recipes.
% Thomas Ferree
% Created 11/17/1999
% Mary Kathryn Reagor
% Revised/optimized 8/28/2007
% Last revised 8/31/2007

fast = 1; % fast = 0 for old code

if fast == 0
    tmp=zeros(n,1);

	% first compute (1/w)*transpose(u)*b when w not = 0
    for j=1:n
        s=0.0;
        if w(j)~=0 % can put higher threshold here
            for i=1:m
                s=s+u(i,j)*b(i);
            end
            s=s/w(j);
        end
        tmp(j)=s;
    end

    for j=1:n
        s=0.0;
        for jj=1:n
            s = s + v(j,jj)*tmp(jj);
        end
        x(j)=s;
    end

end % end if fast = 0

if fast == 1

	% u: 58x58    26912  double               
	% v: 58x58    26912  double               
	% w: 58x1       464  double    

    j=find(w~=0);
        %%******************* Extra Tests for W matrix *******************
        %% MK: checking to see if some w-values need to be zeroed.
        % if(rank(u)<=size(w)) 
        %      error('W is too large');
        % end
        %
        %% MK: seeing if any elements in w are w/in single fl-pt precision
        % e = and((abs(w)<10^-37),(w~=0)); 
        % if sum(e)~=0; 
        %   display(e)
        % end
        %%****************************************************************
    x = v(:,j)*((u(:,j)'*b)./w(j));

end % end if fast = 1
