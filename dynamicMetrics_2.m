%Function that calculates the metastability, dFC as well as whole-brain
%integration and segregation (Q) for a given time series across subjects.
%Inputs are the time series arranged within a cell-array with dimensions:
%1 x numberOFSubjects and delta is the TR of the sample. 
%numberOFSubjects musd be 'ROIs x TimeSeries' format
%Victor Saenger, 2018. UPF, Barcelona.

function [meta, sync_all, dFC, dFC_cos, integ, Q] = dynamicMetrics(ts,delta)

    NSUB = length(ts);
    N = size(ts{1},1);

    %Basic filtering parameters
    %%%%%%%%%%%%%

     flp = .03;   %flp = .01;  0.03; .008;      % lowpass frequency of filter
     fhi = .07;   %fhi = .09;  0.08; .09;       % highpass
     k = 2;                  % 2nd order butterworth filter
     fnq = 1/(2*delta);       % Nyquist frequency
     Wn = [flp/fnq fhi/fnq]; % butterworth bandpass non-dimensional frequency
     [bfilt2,afilt2] = butter(k,Wn);   % construct the filter


    %Filtering method: demeaning, filtering using the band pass filter,
    %Hilbert-transforming the BIOLD signal and calculating the instantaneous
    %phase (line 32) at each node across all timepoints
    for nsub = 1:NSUB
        clear timeseriedata events Phases bpl bplth;
        xs = ts{1,nsub};
        Tmax = size(xs,2); 
        T = 1:Tmax;
        timeseriedata = zeros(N,Tmax);%length(xs)
     for seed = 1:N
          x = demean(detrend(xs(seed,:)));
          timeseriedata(seed,:) = filtfilt(bfilt2,afilt2,x);    % zero phase filter the data
          Xanalytic = hilbert(demean(timeseriedata(seed,:)));
          Phases(seed,:) = angle(Xanalytic); %%% calculating phase of each ROI for each signal
                                             %%% which will use to compute
                                             %%% metastability and other
                                             %%% parameters
     end
    % Metastability: the kuramoto order parameter (ku) is given by the
    % phase difference between the phases at two given time points in the
    % complex plain. The real part of this number is the synchronization
    % level (synch). The standard deviation of synch is in fact the
    % metastability. You could also save synch as an output across
    %subjects to see how it fluctuates over time.
    T = 1:Tmax;
    sync = zeros(1, numel(T));
    for t = T
        ku = sum(complex(cos(Phases(:,t)),sin(Phases(:,t))))/N;
        sync(t) = abs(ku);
    end

    meta(nsub) = std(sync(:));  %std(sync(:));
    sync_all (nsub,:) =sync;
    
    %fCD: This matrix is computed as the cosine difference between all phase
    %pairs.This process sis repeated at each time point so the output is s 3d
    %matrix with n x n x t dimensions.
     for t = T
      for i = 1:N
        for j = 1:N
         dM(i,j,t) = exp(-3*adif(Phases(i,t),Phases(j,t)));  % computes dynamic matrix/ dFC
         dM_cos(i,j,t)=cos(Phases(i,t)-Phases(j,t)); %dM/dFC/phasematrix all same 
         %phasematrix(i,j)=abs(x(i,t)-x(j,t));
        end
      end 
      
      %Integration: Read Gustavo's paper: Novel Intrinsic Ignition Method
     %Measuring Local-Global Integration Characterises Wakefulness and 
     %Deep Sleep. It basically measures the resilience of a given netwrok
     %to pruning. The large the resilience, the bigger the integration
     %will be
      cc = dM(:,:,t); %.*Cbin;
      cc = cc-eye(N);
      pp = 1;
      PR = 0:0.01:0.99;
      
      for p = PR
       A = abs(cc)>p;
       [~, csize] = get_components(A);
       cs(pp) = max(csize);
       pp = pp+1;
      end
      integ(nsub,t) = sum(cs)*0.01/N;
      
      %Q is the modularity index and it is basically the oposite of
      %integration. if you correlate Q ant integ in one given subject you
      %should get a value closer to -1
      [~, Q(nsub,t)] = community_louvain(abs(dM(:,:,t)));
     end
    display(nsub);
    dFC{nsub} = dM; 
    dFC_cos{nsub} = dM_cos; 
    end
end


   
   
   