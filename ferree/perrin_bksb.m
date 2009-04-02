function Cperrin = perrin_bksb(u,w,v,Vmeas)

Nmeas = length(Vmeas);
b = [Vmeas;0.0];
Cperrin = NRsvbksb(u,w,v,Nmeas+1,Nmeas+1,b);

% % compute interpolated potentials
% 
% Vinterp=zeros(Ninterp,1); % uV
% for i=1:Ninterp
% 	Vinterp(i)=c(Nmeas+1);
% 	for j=1:Nmeas
% 		x=dot(Einterp(i,:),Emeas(j,:));
%       Vinterp(i)=Vinterp(i)+c(j)*g_perrin(morder,x);
% 	end
% end	
