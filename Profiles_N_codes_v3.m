%=========================================================================
% This Matlab code can be adopted to benchmark solvers for optimization
% problems. In particular, this code can be used to accomplish the 
% following three specific purposes:
%
%   1) Computing PERFORMANCE PROFILES (E.D.Dolan, J.More',
%      "Benchmarking optimization software with performance profiles",
%      Mathematical Programming 91, pp. 201-213, 2002)
%
%   2) Computing single/multiple DATA PROFILES (J.More', S.Wild,
%      Benchmarking Derivative-Free Optimization Algorithms, SIAM
%      J. Optimization 20, pp. 172-191, 2009)
%
%   3) Computing QUALITY PROFILES (G.Fasano, C.Piermarini, M.Roma,
%      Benchmarking Optimization Algorithms with Quality Profiles,
%      submitted to Mathemtical Programming Computations, 2025)
%
% All the above mentioned profiles can be plotted using either linear or
% logarithmic scale on the abscissa axis. Hereafter we describe the steps 
% for each of the above three tasks, including settings of the required
% parameters in the code.
%
%------------------------------------------------------------------
% TO RUN PERFORMANCE PROFILES the next quantities must be set, all
% the remaining quantities in the code can be ignored. The code will
% automatically detect both the number of test problems and the number 
% of solvers from the provided data.
%------------------------------------------------------------------
%
%   profiles = must be set equal to the value 1
%
%   limit = represents the largest value on the abscissa axis
%           used to plot the performance profiles
%   
%   logscaling = (true) for a logarithmic scale on the abscissa axis
%                (false) for a linear scale on the abscissa axis
%
%   MM = represents a PxS real matrix, being P the number of test
%        problems and S the number of solvers to be compared. The (p,s)
%        entry of MM represents the performance of the solver 's' on the 
%        test problem 'p'. In case the solver 's' has a failure on the 
%        test problem 'p', then provide a very large value (say 10 times 
%        larger than any other value in MM), in the corresponding (p,s) 
%        position of the matrix
%
% The code provides warnings and error messages in case data is
% inconsistent (e.g. some entries of MM are missing or too small, so that 
% performance profile(s) cannot be computed)
%
%
%------------------------------------------------------------------
% TO RUN DATA PROFILES the next quantities must be set, all the 
% remaining quantities in the code can be ignored. The code will
% automatically detect both the number of test problems and the number 
% of solvers from the provided data.
%------------------------------------------------------------------
%
%   profiles = must be set equal to the value 2
%
%   logscaling = (true) for a logarithmic scale on the abscissa axis
%                (false) for a linear scale on the abscissa axis
%
%   FF = represents a PxSxT real matrix, being P the number of test
%        problems, S the number of solvers to be compared and T the
%        number of precision level considered by the user (namely the 
%        number of entries of the vector 'tau'). The entry (p,s,t) of FF 
%        represents the number of function evaluations necessary to the 
%        solver 's' to find a solution to the test problem 'p', within the 
%        precision level 't'. In case the solver 's' has a failure on the 
%        test problem 'p', when using the precision level 't', then 
%        (p,s,t) must be equal to any 'negative number'
%
%  tau = represents the vector of the precision level(s) associated to the 
%        required data profile(s). It can either be a value in (0,1], or 
%        a vector (in the last case it must have an even number of 
%        entries), being each entry in (0,1]
%
%   NV = represents a P vector of positive integers. Each of its entries
%        represents the number of unknowns of test problem 'p'
%  
% The code provides warnings and error messages in case data is
% inconsistent (e.g. on a test problem 'p', the solution found by a 
% solver 's' is larger than the initial objective function value f(x_0^p), 
% so that data profiles cannot be computed)
%
%
%------------------------------------------------------------------
% TO RUN QUALITY PROFILES and TEST SET PROFILES the next quantities 
% must be set, all the remaining quantities in the code can be ignored. 
% The code will automatically detect both the number of test problems 
% and the number of solvers from the provided data.
%------------------------------------------------------------------
%
%   profiles = must be set equal to the value 3
%
%   logscaling = (true) for a logarithmic scale on the abscissa axis
%                (false) for a linear scale on the abscissa axis
%
%   TT = represents a Px(1+S) matrix, being P the number of test
%        problems and S the number of solvers to be compared. The first 
%        column (any real number) gives the value of the objective 
%        function at the initial iterate x_0, while each of the remaining 
%        S columns reports the optimal value (any real number) of the 
%        objective function found by the solver 's' on the test problems
%        in P. 
% 
%   r1 = represents a zooming parameter for the 'abscissa' axis when
%        drawing the plot of a Quality Profile, so that
%           r1 = 1 implies no zooming
%           r1 < 1 implies zooming nearby the abscissa value 1
%           r1 > 1 implies zooming nearby the abscissa value 0 
%
%   r2 = represents a zooming parameter for the 'ordinate' axis when
%        drawing the plot of a Quality Profile, so that
%           r2 = 1 implies no zoomin
%           r2 < 1 implies zooming nearby the ordinate value 1
%           r2 > 1 implies zooming nearby the ordinate value 0 
%
%   strong = must be set to decide between creating either weak or  
%            strong quality profile(s), so that
%               strong = 'true' implies that strong quality profile(s) 
%                        is/are generated
%               strong = 'false' implies that weak quality profile(s) 
%                        is/are generated
%  
%   boot = (true) when bootstrapping on the results is claimed, i.e. the
%          Test Set Profiles will be computed, too
%          (false) when bootstrapping on the results is NOT claimed, i.e. 
%          the Test Set Profiles will NOT be computed
%
%   cycles = number of cycles of bootstrapping to be (in case) considered,
%            for the computation of the Test Set Profiles
%
% The code provides warnings and error messages in case data is
% inconsistent (e.g. on a test problem 'p', the solution found by a 
% solver 's' is larger than the initial objective function value f(x_0^p), 
% so that quality profile(s) cannot be computed)
%
%
%==================================================================
%--------------------------------
% Clear screen and all variables
%--------------------------------
clc; clf;
close all;
clear all;

%-----------------------
% Start to compute time
%-----------------------
tic;

%-----------------------------
% Select the type of profiles
%-----------------------------
%profiles = 1;   % Performance profiles
%profiles = 2;   % Data Profiles
profiles = 3;   % Quality Profiles

% ------------------------------------------------------------------------
% Set the 'limit' for the abscissa parameter of Performance Profiles
% Set the required precision values for Data Profiles
% Set the required scaling parameters 'r1' and 'r2' for Quality Profiles
% Set the parameters 'boot' and 'cycles' for Quality Profiles
% ------------------------------------------------------------------------
limit = 15;         % Upper bound for the abscissa value of the 
                    % Performance Profiles (user's choice)
                   
tau = [0.1];        % Precision level in case only one data profile 
                    % is sought: the user can select any number
                    % in (0,1]
% tau = [0.1 
%        0.01 
%        0.001 
%        0.0001 
%        0.00001 
%        0.000001]; % Precision levels in case multiple (even number) 
                    % of data profiles are sought: any entry must be
                    % selected by the user in (0,1], and will be
                    % associated to a different data profile

r1 = 1; r2=1;       % Scaling parameters for Quality Profiles. Each of
                    % these parameters must be a real positive number       

%strong = true;     % The user creates 'Strong' Quality Profiles
strong = false;     % The user creates 'Weak' Quality Profiles

boot = true;        % bootstrapping required in order to compute the
                    % Test Set Profiles
%boot = false;      % bootstrapping NOT required (i.e. the picture with
                    % the Test Set Profiles must be ignored)
cycles = 1000;      % number of cycles for bootstrapping (ignored if
                    % boot = false)

%-----------------------------------------
% Decide the scaling of the abscissa axis
% for all the types of profiles
%-----------------------------------------
logscaling = true;         % log10 scaling for the abscissa axis
%logscaling = false;       % no log scaling for the abscissa axis


% -----------------------------------------------------------------------
% Matrix MM used to plot Performance Profiles. Its dimensions are PxS with 
%       P = # test Problems 
%       S = # of Solvers
% The entry (p,s) of MM represents the 'performance' of the solve 's' on 
% the problem 'p'. Hence, the entry (p,s) represents (e.g.) the number 
% of iterations, the number of function evaluations, etc. of the solve 
% 's' on the test problem 'p'. In case the solver 's' has a failure on 
% the test problem 'p', then provide a very large value for (p,s)
% -----------------------------------------------------------------------
% MM = [];
 MM = 10 + (100-10)*rand(100,3);      % Just an EXAMPLE for matrix MM
                                      % including 100 test problems
                                      % and 3 solvers. Data is randomly 
                                      % selected

% -----------------------------------------------------------------------
% Matrix FF used to plot Data Profiles. Its dimensions are PxSxT, being 
%       P = # test Problems 
%       S = # of Solvers
%       T = # of precision levels
% The entry (p,s,t) of FF represents the number of function evaluations 
% necessary to the solver 's' to find a solution of the test problem 'p', 
% within the precision level 't'. In case the solver 's' has a failure on 
% the test problem 'p', when using the precision level 't', then (p,s,t)
% must be equal to 'any NEGATIVE number'
% -----------------------------------------------------------------------
% FF = [];
 FF = [     15	    16	22
            22	    25	62
            37	    33	25
            22	    15	32
            48	    28	35
            -1	    84	38
            22	    35	55
            98	    47	98
            43	    65	29
            12	    91	63
            28	    65	61
            100	    56	134
            72	    82	53
            84	    200	142
            33	    178	165
            115	    255	200
            19	    32	28
            102	    205	175
            305	    280	391
            200	    111	124  ];      % Just an EXAMPLE of the matrix FF
                                     % when 'tau' has just one entry
                                     % (courtesy of C.Audet, W.Hare, 
                                     % "Derivative-Free and Blackbox
                                     % Optimization", Springer Series in 
                                     % Operations Research and Financial 
                                     % Engineering, 2017)

% -----------------------------------------------------------------------
% Vector NV used to plot Data Profiles. It is a P dimensional vector  
% where P is the number of test Problems. The entry 'p' of NV represents 
% the number of variables of the test problem 'p'
% -----------------------------------------------------------------------
% NV = [];
 NV = [ 2
        2
        2
        2
        3
        3
        3
        5
        5
        5
        7
        9
        9
        11
        12
        15
        16
        17
        20
        25 ];   % Just an EXAMPLE of the vector NV
                % (courtesy of C.Audet, W.Hare, 
                % "Derivative-Free and Blackbox
                % Optimization", Springer Series in 
                % Operations Research and Financial 
                % Engineering, 2017)

% -----------------------------------------------------------------
% Matrix TT used to plot Quality Profiles and Test Set Profiles. Its 
% dimensions are Px(1+S), being 
%       P = # test Problems 
%       S = # of Solvers
% The columns of TT are (respectively): 
%   1st column  = the objective function value at the initial iterate
%   s-th column = (s >= 2) the objective function value at the optimal 
%                 point detected by the solver 's'
% This matrix will be uniquely used to plot Quality Profiles and Test Set 
% Profiles, and will be ignored if Performance Profiles or Data Profiles 
% are sought. If the solver 's' has a failure on the test problem 'p', 
% then the corresponding entry (p,s) in TT must be equal to 'NaN'
% -----------------------------------------------------------------
% TT = [];
 TT = [(-1 + 2*rand(100,1)) ...
       (-2 + 2*rand(100,1)) (-2.01 + 2*rand(100,1)) (-3.0 + 2*rand(100,1))];
       % Just an EXAMPLE for matrix TT including 100 test problems 
       % and 3 solvers

%------------------------------------------------------------------
% To plot Data Profiles, assign the dimensions of the test problems. 
% To plot Quality Profiles, assig the value of the objective function 
% at the first iterate, for each test problem 
%------------------------------------------------------------------
if (profiles == 2)
    n_p = NV(:);
elseif (profiles == 3)
    f_0 = TT(:,1);
end

%------------------------------------
% Print a summary of the user's data
%------------------------------------
if (profiles == 1)
    [P,S] = size(MM);
elseif (profiles == 2)
    [P,S,entries] = size(FF);
else
    [P,S] = size(TT);
    S = S - 1;
end
if (~boot)
    cycles = 1;
end
fprintf('\n====================================================== \n');
fprintf('  Number of TEST PROBLEMS ------>  %i  \n',P);
fprintf('  Number of SOLVERS       ------>  %i  \n',S);
if (profiles ==1)
    fprintf('  Performance Profiles: user abscissa limit = %8.3f  \n',limit);
elseif (profiles ==2)
    fprintf('  Data Profiles: precision value tau = %d  \n',tau);
elseif(profiles ==3)
    fprintf('  Quality Profiles: scaling parameter (x) r1 =%6.2f  \n',r1);
    fprintf('  Quality Profiles: scaling parameter (y) r2 =%6.2f  \n',r2);
    fprintf('  Quality Profiles: bootstrapping = %i  \n',boot);
    fprintf('  Quality Profiles: cycles = %i  \n',cycles);
end
fprintf('====================================================== \n \n');

% ---------------------------
% Check for data consistency
% ---------------------------
if (profiles == 2)
    [dim_n_p,dummy2] = size(n_p);
elseif (profiles == 3)
    [dim_f_0,dummy1] = size(f_0);
end
if (profiles == 1)
    for i=1:P
        for j=1:S
            if (MM(i,j) <= 0) 
                fprintf('  ERROR: wrong data for Performance Profiles !!!');
                return;
            end
        end
    end
elseif ((profiles == 2) & (min(n_p) < 1))
    fprintf('  ERROR: Wrong problems dimension !!!');
    return;
end
if (profiles == 2)
    [a,b] = size(tau);
    if (b ~= 1) & (mod(b,2) ~= 0) 
        fprintf('  ERROR: odd number of required Data Profiles !!!');
        return;
    end
end

% -----------------------------------------------------------
% In case of Quality Profiles then
%       - select the stepsize for the parameter representing 
%         the precision 'tau'
%       - update the vector f_0 if Strong Quality Profiles 
%         are sought
% -----------------------------------------------------------
if (profiles ==3)
    intervals = floor(max([1000 1000*r1 1000/r1]));
    stepsize = 1/intervals;
    tau = zeros(intervals+1,1);
    if (strong == true)
        for i=1:P
            f_0(i) = max(TT(i,2:(2+S-1)));
        end
    end
end

% -----------------------------------------------------------------
% Compute the best value among the solvers, for each test problem. 
% The best value is only used for Quality Profiles. 
% -----------------------------------------------------------------
if (profiles == 3)
    f_star = zeros(P,1);
    for p=1:P
        f_star(p) = min(TT(p,2:(2+S-1)));
        if (f_star(p) > f_0(p))
            fprintf('  ERROR: wrong data for Quality Profiles !!!');
            return;
        end
    end
end

% ------------------------------
% Compute the required profiles
% ------------------------------
% ----------------------
%  PERFORMANCE PROFILES
% ----------------------
if (profiles == 1)  
    % Here we create the Performance Profiles "RHO_s(tau)"
    % for any solver 's'. The parameter 'tau' has no dimension and
    % represents the abscissa parameter in Performance Profiles
    RHO = zeros(S,1000);
    tps = MM;
    rps = zeros(P,S);
    r_M = 0;        % initialize the maximum performance ratio
    min_tps = zeros(P,1);
    for p=1:P
        min_tps(p) =  min(MM(p,:));
        for s=1:S
           if (abs(tps(p,s)) <= 1.0e-10)
                fprintf('WARNING: very small data in Performance Profiles !!!');
                pause;
           end
           rps(p,s) =  tps(p,s) / min_tps(p);
           r_M = max(r_M , rps(p,s));
        end    
    end
    r_M = 1.1*min(r_M,limit);   % compute an upper bound to 'r_M' so that 
                                % the plot of the Performance Profiles
                                % will correctly fit the user's interval 
    tau = linspace(1, r_M, 1000);
    for t=1:1000
        for s=1:S
            counter = 0;
            for p=1:P
                if (rps(p,s) <= tau(t)) 
                    counter = counter + 1;
                end 
            end
            RHO(s,t) = counter/P;
        end
    end

% ---------------
%  DATA PROFILES
% ---------------
elseif (profiles == 2)  
    % Here we create the Data Profiles "d_s(alpha)" for any 
    % solver 's'. The parameter 'alpha' has no dimension and
    % represents a multiple of the cost of one simplex gradient 
    % evaluation, i.e. n_p + 1, being 'n_p' the number of 
    % variables of the test problem 'p'
    [entries,dummy] = size(tau);
    d = zeros(S,5000,entries);
    tps = FF;
    i_M = ones(entries,1);
    
    % Compute the maximum performance ratio and the smallest interval for 
    % each precision level
    delta = zeros(entries);
    r_max = zeros(entries);
    for t=1:entries
        for p=1:P
             for s=1:S
                 r_max(t) = max(r_max(t) , tps(p,s,t)/(n_p(p)+1));
             end
        end 
    end
    for t=1:entries
        r_max(t) = 1.05*r_max(t);     % Compute an upper bound on 'r_M'
        delta(t) = r_max(t)/5000;     % Compute an elementary interval
    end
    for t=1:entries         % Compute Data Profiles for each precision value
        for s=1:S           % For any solver 's'
            for i=1:5000    % Increase the abscissa parameter 'alpha' of
                            % a small interval 'delta(t)'
                counter = 0;
                for p=1:P
                    if ( (tps(p,s,t) > 0) & (tps(p,s,t)/(n_p(p)+1) <= i*delta(t)) )
                         counter = counter + 1;
                    end
                end
                d(s,i,t) = counter/P;
            end    
        end
    end
    
    % Compute the largest abscissa such that, for a given entry of "tau", 
    % the plot of Data Profiles for the solvers does not change
    for t=1:entries
        for i=1:4999
            i_found = false;
            for s=1:S
                if ( abs(d(s,5000-i+1,t) - d(s,5000-i,t)) > 1.0e-5 )
                    i_M(t) = floor( min((5000-i+1)*1.1,5000) );
                    i_found = true;
                    break
                end
            end
            if i_found
                break
            end
        end
    end
   
% ------------------
%  QUALITY PROFILES
% ------------------
else
    % Here we create the Quality Profiles "Qs(tau)". Each row of the 
    % matrix 'Q' is associated with a different precision level.
    % Each column of the matrix 'Q' is associated with a different solver.
    % Moreover, we create the matrix "QQs(tau,c)" such that: in case of 
    % bootstrapping, for any cycle 'c' we store QQs(tau,c)=Qs(tau). 
    Q = zeros(intervals+1,S);
    QQ = zeros(intervals+1,S,cycles);
    f_star_w = zeros(size(f_star)); % create a dummy vector to work
                                    % within bootstrapping cycles
    f_0_w = zeros(size(f_0));       % create another a dummy vector to work
                                    % within bootstrapping cycles
    TT_work = zeros(size(TT));      % create a dummy matrix to work
                                    % within bootstrapping cycles
    % Perform the first cycle of bootstrapping
    % fprintf('  Cycle = %i / %i  \n',1,cycles);
    for s=1:S
            for t=1:(intervals+1)
                tau(t) = ((t-1)*stepsize)^(r1);
                counter = 0;
                for w=1:P
                    if (~isnan(f_0(w)) & ~isnan(TT(w,1+s)))
                        if ((TT(w,1+s) - f_star(w)) <= tau(t)*(f_0(w) - f_star(w)))
                            counter = counter + 1;
                        end
                    end
                end
                % Scaling, through the parameter 'r2', of the Quality Profiles 
                Q(t,s) = (counter/P)^(r2);
            end
    end
    QQ(:,:,1) = Q(:,:); 

    % Perform (if required) the other cycles of bootstrapping
    if (boot)
        for c=2:cycles
            % fprintf('  Cycle = %i / %i  \n',c,cycles);
            % Shuffle the test problems to generate the novel test set in the
            % current bootstrapping cycle
            % u = floor(1 + (P-1).*rand(P,1));    % this formula does NOT allow
                                                  % to select the p-th problem
            u = floor(1 + P.*rand(P,1));          % this formula allows to
                                                  % select also the p-th problem
            for w=1:P
                TT_work(w,:) = TT(u(w),:);
                f_star_w(w) = f_star(u(w));
                f_0_w(w) = f_0(u(w));
            end
    
            for s=1:S
                for t=1:(intervals+1)
                    tau(t) = ((t-1)*stepsize)^(r1);
                    counter = 0;
                    for w=1:P
                        if (~isnan(f_0_w(w)) & ~isnan(TT_work(w,1+s)))
                            if ((TT_work(w,1+s) - f_star_w(w)) <= tau(t)*(f_0_w(w) - f_star_w(w)))
                                counter = counter + 1;
                            end
                        end
                    end
                    % Scaling, through the parameter 'r2', of the Quality Profiles 
                    Q(t,s) = (counter/P)^(r2);
                end
            end
            QQ(:,:,c) = Q(:,:);
        end
    end
end

%----------------------------%
% Plot the required profiles %
%----------------------------%

%---------------------------%
% Plot PERFORMANCE PROFILES %
%---------------------------%
if (profiles ==1)     
    figure1 = figure;
        % colors = [linspace(1, 0, S)', rand(S,1), linspace(0, 1, S)'];
        colors = [linspace(1, 0, S)', zeros(S,1), linspace(0, 1, S)'];
        if (logscaling)
            for s=1:S
                semilogx(tau,RHO(s,:),':','Color',colors(s,:),'LineWidth',1.5);
                hold on
            end
        else
            for s=1:S
                plot(tau,RHO(s,:),':','Color',colors(s,:),'LineWidth',1.5);
                hold on
            end
        end
        title(['{\color{red}Performance Profiles:} |S| = ',num2str(S), ...
               '; |P| = ',num2str(P), ...
               '; \tau--range = [1 , ',num2str(r_M),']'],'FontName','Baskerville Old Face', ...
               'Fontsize',35);
        grid on
        xlabel('\it \tau','FontWeight','bold','FontSize',35); 
        ylabel(' \it \rho_s(\tau)','FontWeight','bold','FontSize',30);
        set(gca, 'FontSize', 15);
        % Here the user can specify the name of the solvers in the vector 
        % "SolversNames": the following example is provided 
        SolversNames = strseq(' Solver - ',1:S);
        legend1 = legend(SolversNames);
        set(legend1,'Location','southeast','FontSize',20);
        hold off
        
% --------------------
%  Plot DATA PROFILES
% --------------------
elseif (profiles ==2)  
    figure2 = figure;
        % colors = [linspace(1, 0, S)', rand(S,1), linspace(0, 1, S)'];
        colors = [linspace(1, 0, S)', zeros(S,1), linspace(0, 1, S)'];
        % ---------------------------
        %  create 'entries' subplots
        % ---------------------------
        if (entries ==1)
                if (logscaling)
                    for s=1:S
                        semilogx((1:i_M(1))*delta(1),d(s,1:i_M(1),1),':','Color',colors(s,:),'LineWidth',1.5);
                        hold on
                    end
                else
                    for s=1:S
                        plot((1:i_M(1))*delta(1),d(s,1:i_M(1),1),':','Color',colors(s,:),'LineWidth',1.5);
                        hold on
                    end
                end
                title(['{\color{red}Data Profile:} |S| = ', ...
                        num2str(S),';  |P| = ', ...
                        num2str(P),';  \tau = ',num2str(tau(1))],'FontName','Baskerville Old Face', ...
                        'Fontsize',35);
                grid on
                xlabel('\it \alpha','FontWeight','bold','FontSize',35); 
                ylabel(' \it d_s(\alpha)','FontWeight','bold','FontSize',30);
                set(gca, 'FontSize', 15);
                % Here the user can specify the name of the solvers in the 
                % vector "SolversNames": the following example is provided
                SolversNames = strseq(' Solver - ',1:S);
                legend2 = legend(SolversNames);
                set(legend2,'Location','southeast','FontSize',10);
        else
            for j=1:entries
                subplot(entries/2,2,j)
                if (logscaling)
                    for s=1:S
                        semilogx((1:i_M(j))*delta(j),d(s,1:i_M(j),j),':','Color',colors(s,:),'LineWidth',1.5);
                        hold on
                    end
                else
                    for s=1:S
                        plot((1:i_M(j))*delta(j),d(s,1:i_M(j),j),':','Color',colors(s,:),'LineWidth',1.5);
                        hold on
                    end
                end
                title(['{\color{red}Data Profile:} |S| = ', ...
                        num2str(S),';  |P| = ', ...
                        num2str(P),';  \tau = ',num2str(tau(j))],'FontName','Baskerville Old Face', ...
                        'Fontsize',20);
                grid on
                xlabel('\it \alpha','FontWeight','bold','FontSize',20); 
                ylabel(' \it d_s(\alpha)','FontWeight','bold','FontSize',20);
                set(gca, 'FontSize', 13);
                % Here the user can specify the name of the solvers in the 
                % vector "SolversNames": the following example is provided
                SolversNames = strseq(' Solver - ',1:S);
                legend2 = legend(SolversNames);
                set(legend2,'Location','southeast','FontSize',10);
            end
        end

% -----------------------
%  Plot QUALITY PROFILES
% -----------------------
else
    figure3 = figure;   % PLOT of the QUALITY PROFILE
        % colors = [linspace(1, 0, S)', rand(S,1), linspace(0, 1, S)'];
        colors = [linspace(1, 0, S)', zeros(S,1), linspace(0, 1, S)'];
        if (logscaling)
            idx = 2:(intervals+1);
            for s=1:S
                semilogx(tau(idx),Q(idx,s),':','Color',colors(s,:),'LineWidth',2.5);
                %loglog(tau(idx),Q(idx,s),':','Color',colors(s,:),'LineWidth',2.5);
                hold on
            end
        else
            for s=1:S
                plot(tau(1:intervals+1),Q(:,s),':','Color',colors(s,:),'LineWidth',2.5);
                hold on
            end
            [rows,cols] = size((0:stepsize:1)');
            plot(0.001*ones(rows,1),(0:stepsize:1)','-k',...
                 0.01*ones(rows,1),(0:stepsize:1)','-k',...
                 0.1*ones(rows,1),(0:stepsize:1)','-k','LineWidth',0.5);
        end
        if (strong)
            title([' {\color{red}(Strong) Quality Profiles:} |S| = ', ...
                            num2str(S),';  |P| = ', ...
                            num2str(P),';  x-scaling (r_1 = ',num2str(r1), ...
                   ');  y-scaling (r_2 = ',num2str(r2),')'],'FontName','Baskerville Old Face', ...
                   'Fontsize',35);
        else
            title([' {\color{red}(Weak) Quality Profiles:} |S| = ', ...
                            num2str(S),';  |P| = ', ...
                            num2str(P),';  x-scaling (r_1 = ',num2str(r1), ...
                   ');  y-scaling (r_2 = ',num2str(r2),')'],'FontName','Baskerville Old Face', ...
                   'Fontsize',35);
        end
        grid on
        xlabel('0 \leq \tau^{r_1} \leq 1','FontWeight','bold','FontSize',35); 
        ylabel(' [ Q_s(\tau) ]^{r_2}','FontWeight','bold','FontSize',35);
        set(gca, 'FontSize', 15);
        % Here you can specify the name of the solvers in the 
        % vector "SolversNames": the following example is provided
        SolversNames = strseq(' Solver - ',1:S);
        legend3 = legend(SolversNames);
        set(legend3,'Location','southeast','FontSize',15);
        hold off
      
        if (boot)
          figure4 = figure;         % PLOT of the TEST SET PROFILE
             % compute the matrix of variances
             QQ_var = zeros(intervals+1,S);
             for t=1:(intervals+1)
                 for s=1:S
                     QQ_var(t,s) = var(QQ(t,s,:));
                 end
             end
    
            % colors = [linspace(1, 0, S)', rand(S,1), linspace(0, 1, S)'];
            colors = [linspace(1, 0, S)', zeros(S,1), linspace(0, 1, S)'];
            if (logscaling)
                idx = 2:(intervals+1);
                for s=1:S
                    semilogx(tau(idx),QQ_var(idx,s),':','Color',colors(s,:),'LineWidth',2.5);
                    %loglog(idx),QQ_var(idx,s),':','Color',colors(s,:),'LineWidth',2.5);
                    hold on
                end
            else
                for s=1:S
                    plot(tau(1:intervals+1),QQ_var(:,s),':','Color',colors(s,:),'LineWidth',2.5);
                    hold on
                end
                [rows,cols] = size((0:stepsize:1)');
                plot(0.001*ones(rows,1),(0:stepsize:1)','-k',...
                     0.01*ones(rows,1),(0:stepsize:1)','-k',...
                     0.1*ones(rows,1),(0:stepsize:1)','-k','LineWidth',0.5);
            end
            if (strong)
                title([' {\color{red}    Var (Strong) Quality Profiles:} |S| = ', ...
                                num2str(S),'; |P| = ', ...
                                num2str(P),'; x-scal. (r_1 = ',num2str(r1), ...
                       '); y-scal. (r_2 = ',num2str(r2),')'],'FontName','Baskerville Old Face', ...
                       'Fontsize',35);
            else
                title([' {\color{red}    Var (Weak) Quality Profiles:} |S| = ', ...
                                num2str(S),'; |P| = ', ...
                                num2str(P),'; x-scal. (r_1 = ',num2str(r1), ...
                       '); y-scal. (r_2 = ',num2str(r2),')'],'FontName','Baskerville Old Face', ...
                       'Fontsize',35);
            end
            grid on
            xlabel('0 \leq \tau^{r_1} \leq 1','FontWeight','bold','FontSize',35); 
            ylabel(' \sigma^2_s(\tau)','FontWeight','bold','FontSize',35);
            set(gca, 'FontSize', 15);
            % Here you can specify the name of the solvers in the 
            % vector "SolversNames": the following example is provided
            SolversNames = strseq(' \sigma^2\_Solver - ',1:S);
            legend4 = legend(SolversNames);
            set(legend4,'Location','northwest','FontSize',15);
            hold off
        end
end

   
%----------------------
% Stop to compute time
%----------------------
toc;
