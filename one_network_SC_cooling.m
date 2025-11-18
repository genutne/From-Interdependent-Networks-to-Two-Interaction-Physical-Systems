clear
% Evaluation of the effective resistance in a superconducting square-grid network 
% under constant applied current, as a function of varying heat-bath temperature.


% general structure of the network
L = 100; % length size
W = L; % width size
NR = (2*W-1)*L; % number of links (without entrance)
NN = W*L; % number of nodes
epsilon = 10^(-5); % threshold for epsilon to determine that the system converged
gamma = 0.1*10^9; %gamma describes the relationship between heat and temperature in a system, depending on the medium
Rhot = 500; % the resistance in the Normal phase
RSc = 10^(-5); % the resistance in the SC (Super-Conducting) phase
Ic0 = 58*10^(-6); % avg critical current
sigma = 0.1; % variance of the criticality of the system resistors
Dt = 1000; % Dt sets the level of heat diffusion. Dt~L^2 = mean-filed aproximation
Itotal = 400*10^(-6); % fixed current through the system
Ti = 0.1; % lowest temperature
Tf = 3; % highest temperature
Tjumps = 1000;  
heat = 0; % chose 0 for cooling the system and 1 for heating.
plateau_vec = nan;

%if a plateau mesurment wanted, as Tc is known, activate the lines:
% Tc = 1;
% dT = [0:150]*0.01;
% if heat
%     plateau_vec = Tc+dT;
% elseif ~heat
%     plateau_vec = Tc-dT;
% end

if isnan(plateau_vec)
    plateau = 0;
else
    plateau = 1;
end  

if plateau
    pl = struct('T',[],'iter',[],'R',[]);
    plateau_structor(length(plateau_vec)) = pl;
end

for p = 1:length(plateau_vec)

    % initial properties of the network

    I = zeros(1,NN+1); % seting the general currents vector
    I(end) = Itotal;

    cin = zeros(1,W); % vector currents for each link
    currents = zeros(1,NR);

    if heat % vector resistors for each link
        Rin = RSc*ones(1,W);
        R = RSc*ones(1,NR);
    elseif ~heat
        Rin = Rhot*ones(1,W);
        R = Rhot*ones(1,NR);
    end

    states_in = zeros(1,W); % vector states for each link. 1=SC , 2=intermidiate , 3=normal.
    states = zeros(1,NR);

    Ic_in = Ic0*(1+sigma*randn(1,W)); % Ic for each resistor - creating the network
    Ic = Ic0*(1+sigma*randn(1,NR));

    % Tc for each resistor
    factor = 100;  % the 100 here can be modified but should be close to it to match the experiment
    Tc = factor*Rhot*Ic;
    Tcin = factor*Rhot*Ic_in;

    %temperature properties
    T_range = linspace(Ti,Tf,Tjumps);
    if ~heat
        flip(T_range);
    end
    if plateau %for plateau
        if heat
            T_range = [Ti plateau_vec(p)]; 
        elseif ~heat
            T_range = [Tf plateau_vec(p)];
        end
    end
    R_all = zeros(1,length(T_range));

    for i=1:length(T_range)

        T0 = T_range(i); % heat bath temperature T0
        finish = 0;
        iter = 0;
        flag = 0;                
        iter_per_T = [];
        R_per_T = [];
        cc1 = [];
        cc2 = [];

        while flag==0 
            iter = iter+1;

            %saving 2 last voltege results for exit condition
            if iter>2
                Vold2 = Vold;
            end

            if iter>1
                Vold = V;
            end

            % first, calculate the voltage at each node in the two networks    
            V = Voltage_nodes(L, W, NN, NR,Rin, R,I);

            % checking exit condition
            if iter>2
                if norm(V-Vold2)/norm(V)<epsilon || norm(V-Vold)/norm(V)<epsilon
                    flag=1;
                end
                cc1 = [cc1 norm(V-Vold2)/norm(V)];
                cc2 = [cc2 norm(V-Vold)/norm(V)];
                if flag == 0 && iter>10000 % to avoide numeric-loops, checking the exit condition for the mean error
                    if mean(cc1(end-1000:end-500))-mean(cc1(end-500:end-2))<10^(-5) || mean(cc2(end-1000:end-500))-mean(cc2(end-500:end))<10^(-5)
                        flag=1;
                        ["Number of iterations is very high"]
                    end
                end
            end

            % now, calculate the currents in each resistor and its heat production
            [currents,cin,heat,heat_in] = find_currents(V, L, W, NN, NR, Rin, R,1);              

            % let the heat diffuse (only in the bulk, no exit/entrence nodes)       
            Q = heat_transfer_spectral(L,W,Dt,heat(1:NR-W)*gamma);
            heat(1:NR-W) = Q;

            % find the new heated resistance and states
            [R , Rin, states ,states_in] = update_resistors_SC(NR, W, R, Rin, heat, heat_in, Ic, Ic_in, Tc, Tcin, Rhot, RSc, currents, cin,T0, states ,states_in);

            iter_per_T = [iter_per_T iter];
            R_per_T = [R_per_T V(end)/Itotal];
        end

        if finish
            break
        end

        %%%% resistence
        R_all(i) = V(end)/Itotal;
        [T0 R_all(i) iter] %true-time-state
    end

    if plateau
        plateau_structor(p).T = T0;
        plateau_structor(p).iter = iter_per_T;
        plateau_structor(p).R = R_per_T;
    end
end

if heat
    description ='heating';
elseif ~heat
    description ='cooling';
end

if plateau
    %ploting plateau results
    figure; hold on;
    for iii = 1:length(plateau_vec)
        legendlist{iii} = ['T = ', num2str(plateau_structor(iii).T)];
        plot(plateau_structor(iii).iter, plateau_structor(iii).R);    
    end   
    title(['the plateau effect, for: Ib=',num2str(Itotal*(10^6)), ', L=' ,num2str(L),' one layer']);
    legend(legendlist, 'Location', 'best');
    ylabel('R (\Omega)');
    xlabel('iterations');
else    
    %ploting regular results
    figure;
    hold on;
    plot(T_range, R_all,'-o');
    title(['Resistence by ',description,': ', 'I0 =', num2str(Itotal*10^6), ', Dt=', num2str(Dt), ', L=' , num2str(L), ', kappa=',num2str(gamma), ', one layer']);
    ylabel('R (\Omega)');
    xlabel('T(K)'); 
end

%save results
if plateau
    %saving plateau results
    filename =  [description,'_sigma' ,num2str(sigma),'Iin' ,num2str(Itotal*10^6),'Dt', num2str(Dt), 'L',num2str(L), 'k',num2str(gamma), '_onelayer_plateau.mat'];
    save(filename, 'plateau_structor');
else
    %saveing regular results
    filename = [description,'_sigma' ,num2str(sigma),'Iin' ,num2str(Itotal*10^6),'Dt', num2str(Dt), 'L',num2str(L), 'k',num2str(gamma), '_onelayer.mat'];
    save(filename, 'T_range','R_all','iterations');
end
