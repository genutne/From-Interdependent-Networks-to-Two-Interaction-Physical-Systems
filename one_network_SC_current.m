clear
% Evaluation of the effective resistance in a superconducting square-grid network 
% under constant applied heat-bath temperature, as a function of varying current.

% general structure of the network
L = 100; % length size
W = L; % width size
NR = (2*W-1)*L; % number of links (without entrance)
NN = W*L; % number of nodes
epsilon = 10^(-5); % threshold for epsilon to determine that the system converged
gamma = 2.5*10^9; % gamma describes the relationship between heat and temperature in a system, depending on the medium
Rhot = 500; % the resistance in the Normal phase
RSc = 10^(-5); % the resistance in the SC (Super-Conducting) phase
Ic0 = 58*10^(-6); % avg critical current
sigma = 0.1; % variance of the criticality of the system resistors
Dt = 1000; % Dt sets the level of heat diffusion
T0 = 2; %fixed hea-bath temperature
lower = 1; % (1) for lowering the current each step (0) for increasing
Ii = 10; % lowest current
If = 300; % highet current
Ijumps = 1000;
plateau_vec = nan;

I_vec = linspace(Ii,If,Ijumps)*10^(-6);
if lower
    I_vec = flip(I_vec);
end

% % if a plateau mesurment wanted, as Ic is known, activate the lines:
% Ic = 1;
% x = 1:100;
% f = @(x) (10.^(-4))*exp(0.12*x);
% dI = f(x);
% if lower
%     plateau_vec = Ic-dI;
% elseif ~lower
%     plateau_vec = Ic+dI;
% end
    
if isnan(plateau_vec)
    plateau = 0;
else
    plateau = 1;
end  

if plateau
    pl = struct('I',[],'iter',[],'R',[]);
    plateau_structor(length(plateau_vec)) = pl;
end

for p = 1:length(plateau_vec)

    %current properties
    if plateau %for plateau
        I_vec = [Ii plateau_vec(p)]*10^(-6); 
        if lower
            I_vec = [If plateau_vec(p)]*10^(-6); 
        end
    end
    LI = length(I_vec);
    iterations = zeros(1,LI);
    R_all = zeros(1,LI);

    % vector currents for each link
    cin = zeros(1,W);
    currents = zeros(1,NR);

    % vector resistors for each link. we start with all at the Normal phase
    if lower
        Rin = Rhot*ones(1,W);
        R = Rhot*ones(1,NR);
    else
        Rin = RSc*ones(1,W);
        R = RSc*ones(1,NR);
    end

    % vector states for each link. 1=SC , 2=intermidiate , 3=normal.
    states_in = zeros(1,W);
    states = zeros(1,NR);

    % Ic for each resistor
    Ic_in = Ic0*(1+sigma*randn(1,W));
    Ic = Ic0*(1+sigma*randn(1,NR));

    % Tc for each resistor
    factor = 100;  % the 100 here can be modified but should be close to it to match the experiment
    Tc = factor*Rhot*Ic;
    Tcin = factor*Rhot*Ic_in;

    for i=1:length(I_vec)

        Itotal = I_vec(i); % heat bath temperature T0
        finish = 0;
        iter = 0;
        flag = 0;                
        iter_per_I = [];
        R_per_I = [];
        cc1 = [];
        cc2 = [];            

        %seting the general currents vector
        I = zeros(1,NN+1);
        I(end) = Itotal;

        while flag==0 
            iter = iter+1;

            %saving 2 last voltege results for exit condition
            if iter>2
                Vold2 = Vold;
            end

            if iter>1
                Vold = V;
            end

            % first, calculate the voltage at each node in the network    
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

            %let the heat diffuse (only in the bulk, no exit/entrence nodes)       
            Q = heat_transfer_spectral(L,W,Dt,heat(1:NR-W)*gamma);
            heat(1:NR-W) = Q;

            %Updates the states and resistance of the resistors following the updated temperature
            [R , Rin, states ,states_in] = update_resistors_SC(NR, W, R, Rin, heat, heat_in, Ic, Ic_in, Tc, Tcin, Rhot, RSc, currents, cin,T0, states ,states_in);

            iter_per_I = [iter_per_I iter];
            R_per_I = [R_per_I V(end)/Itotal];

        end

        if finish
            break
        end

        %%%% resistence
        R_all(i) = V(end)/Itotal;
        iterations(i) = iter;

        [Itotal*10^6 R_all(i) iter]

    end
    if plateau
        plateau_structor(p).I = Itotal;
        plateau_structor(p).iter = iter_per_I;
        plateau_structor(p).R = R_per_I;
    end
end

if lower
    PT_direction ='decrease_current';
else
    PT_direction ='increase_current';
end
if plateau
    %ploting plateau results
    figure; hold on;
    for iii = 1:length(plateau_vec)
        legendlist{iii} = ['I = ', num2str(plateau_structor(iii).I.*10^6)];
        plot(plateau_structor(iii).iter, plateau_structor(iii).R,'-o');    
    end   
    title(['the plateau effect for: T0=',num2str(T0), ', L=' ,num2str(L),' one layer ',PT_direction]);
    legend(legendlist, 'Location', 'best');
    ylabel('R (\Omega)');
    xlabel('iterations');
else    
    %ploting regular results
    figure;
    hold on;
    plot(I_vec.*10^6, R_all,'-o');
    title(['Resistence  T0 =', num2str(T0), ', Dt=', num2str(Dt), ', L=' , num2str(L), ', kappa=',num2str(gamma), ', one layer ',PT_direction]);
    ylabel('R (\Omega)');
    xlabel('I(\muA)'); 
end

%save results
if plateau
    %saving plateau results
    filename =  [PT_direction ,'_sigma' ,num2str(sigma),'T' ,num2str(T0),'Dt', num2str(Dt), 'L',num2str(L), 'k',num2str(gamma), '_onelayer_plateau.mat'];
    save(filename, 'plateau_structor');
else
    %saveing regular results
    filename = [PT_direction ,'_sigma' ,num2str(sigma),'T' ,num2str(T0),'Dt', num2str(Dt), 'L',num2str(L), 'k',num2str(gamma), '_onelayer.mat'];
    save(filename, 'I_vec','R_all','iterations');
end
