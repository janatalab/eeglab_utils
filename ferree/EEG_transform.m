function [EEG_interp, EEG_averef, EEG_surflap] = EEG_transform(EEG, EEG_interplocs)

% function [EEG_interp, EEG_averef, EEG_surflap] = EEG_transform(EEG, EEG_interplocs)
% Transforms data in EEGLAB structure to Perrin-based ave ref and surf Lap
% Takes as input EEG data structure for either .cnt or .eeg file.
% Optional input is different channel structure for interpolation.
% Output is interpolated raw data, spline-based average reference,
% and spline-based surface Laplacian.
% Thomas Ferree
% Created 5/23/2006 from older codes
% Last revised 11/20/2007

Nmeas = EEG.nbchan;
Ntime = EEG.pnts;
Ntrials = EEG.trials;
if size(size(EEG.data),2) == 2
    datatype = 'cnt';
elseif size(size(EEG.data),2) == 3
    datatype = 'eeg';
end

% set up Rmeas array
Rmeas = zeros(Nmeas,3);
for k = 1:Nmeas
    Rmeas(k,1) = EEG.chanlocs(k).X;
    Rmeas(k,2) = EEG.chanlocs(k).Y;
    Rmeas(k,3) = EEG.chanlocs(k).Z;
end

% normalize Rmeas to unit sphere
for k = 1:Nmeas
    m = mag(Rmeas(k,:));
    Rmeas(k,:) = Rmeas(k,:) / m;
end

% front calculations for spherical spline interpolation
morder = 3;
lambda = 0;
display('Decomposing Perrin spline matrix...');
[u w v] = perrin_gsvd(Rmeas, Rmeas, morder, lambda);

% compute spline coefficients through time
display(['Computing Perrin spline coefficients at ' num2str(Ntrials*Ntime) ' time points...']);
if strcmp(datatype,'cnt')
    Cperrin = zeros(Nmeas+1, Ntime);
    for t = 1:Ntime
        Cperrin(:,t) = perrin_bksb(u,w,v,EEG.data(:,t));
    end
elseif strcmp(datatype,'eeg')
    Cperrin = zeros(Nmeas+1, Ntime, Ntrials);
    for epoch = 1:Ntrials
        for t = 1:Ntime
            Cperrin(:,t,epoch) = perrin_bksb(u,w,v,EEG.data(:,t,epoch));
        end
    end
end

% contruct interpolation array
if nargin == 1
    Ninterp = Nmeas;
    Rinterp = Rmeas;
elseif nargin == 2
    Ninterp = length(EEG_interplocs);
    % set up Rinterp array
    Rinterp = zeros(Ninterp,3);
    for k = 1:Ninterp
        Rinterp(k,1) = EEG_interplocs(k).X;
        Rinterp(k,2) = EEG_interplocs(k).Y;
        Rinterp(k,3) = EEG_interplocs(k).Z;
    end
    % normalize Rinterp to unit sphere
    for k = 1:Ninterp
        m = mag(Rinterp(k,:));
        Rinterp(k,:) = Rinterp(k,:) / m;
    end
end

% fill G and H arrays (note updated dimensions)
g = zeros(Ninterp,Nmeas);
h = zeros(Ninterp,Nmeas);
for i = 1:Ninterp
    for j = 1:Nmeas
        xij = dot(Rinterp(i,:),Rmeas(j,:));
        g(i,j) = g_perrin(morder, xij);
        h(i,j) = h_perrin(morder,xij);
    end
end

% interpolate data onto EEG_interplocs
if nargin == 1
    Vinterp = EEG.data;
elseif nargin == 2
    display('Interpolating data onto output array...');
    if strcmp(datatype,'cnt')
        Vinterp = zeros(Ninterp,Ntime);
        for t = 1:Ntime
            Vinterp(:,t) = Cperrin(Nmeas+1,t) + g * Cperrin(1:Nmeas,t);
        end
    elseif strcmp(datatype,'eeg')
        Vinterp = zeros(Ninterp,Ntime,Ntrials);
        for epoch = 1:Ntrials
            for t = 1:Ntime
                Vinterp(:,t,epoch) = Cperrin(Nmeas+1,t,epoch) + g * Cperrin(1:Nmeas,t,epoch);
            end
        end
    end
end

% compute Perrin average reference
display('Computing Perrin average reference...');
if strcmp(datatype,'cnt')
    AveRefC0 = zeros(Ninterp,Ntime);
    for t = 1:Ntime
        AveRefC0(:,t) = Vinterp(:,t) - Cperrin(Nmeas+1,t);
    end
elseif strcmp(datatype,'eeg')
    AveRefC0 = zeros(Ninterp,Ntime,Ntrials);
    for epoch = 1:Ntrials
        for t = 1:Ntime
            AveRefC0(:,t,epoch) = Vinterp(:,t,epoch) - Cperrin(Nmeas+1,t,epoch);
        end
    end
end

% compute Perrin surface Laplacian
display('Computing Perrin surface Laplacian...');
if strcmp(datatype,'cnt')
    SLinterp = zeros(Ninterp,Ntime);
    for t = 1:Ntime
        SLinterp(:,t) = SLinterp(:,t) + h * Cperrin(1:Nmeas,t);
    end
elseif strcmp(datatype,'eeg')
    SLinterp = zeros(Ninterp,Ntime,Ntrials);
    for epoch = 1:Ntrials
        for t = 1:Ntime
            SLinterp(:,t,epoch) = SLinterp(:,t,epoch) + h * Cperrin(1:Nmeas,t,epoch);
        end
    end
end

% put surface Laplacian in proper units
sradius = 9.2/100; % standard scalp radius in meters
SLinterp = SLinterp/(sradius)^2; % uV/m^2
SLinterp = SLinterp/10000; % uV/cm^2

% append EEG stucture to acccount for electrode replacement
if nargin == 2
    EEG.nbchan = Ninterp;
    EEG.chanlocs = EEG_interplocs;
end

% EEG structure outputs
EEG_interp = EEG;
EEG_interp.data = Vinterp;
EEG_averef = EEG;
EEG_averef.data = AveRefC0;
EEG_surflap = EEG;
EEG_surflap.data = SLinterp;
