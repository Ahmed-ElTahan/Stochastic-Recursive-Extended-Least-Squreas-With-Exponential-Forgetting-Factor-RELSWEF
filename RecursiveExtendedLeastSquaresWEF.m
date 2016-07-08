% This function is made by Ahmed ElTahan

%{
        This function is intended to estimate the parameters of a dynamic
        system of unknown parameters using the Recursive Extended Least Squares With Exponential Forgetting Factor
        Method (RELSWEF) for time varying parameter system which has an noise addition.
        After an experiment, we get the inputs, the outputs of the system.
        The experiment is operated with sample time Ts seconds.
                                    
        The model is given by           A(z) y(t) = B(z)sys u(t) + C(z) eps(t)
        which can be written in
                                                                     z^(-d) B(z)                C(z)
                                                         y(t) = ------------------- u  + ------------ e = L*u + M*e
                                                                         A(z)                       A(z)

    where:
    -- y : output of the system.
    -- u : control action (input to the system).
    -- e : color guassian noise (noise with non zero mean).
    -- Asys = 1 + a_1 z^-1 + a_2 z^-2 + ... + a_na z^(-na). [denominator polynomail]
    -- Bsys = b_0 + b_1 z^-1 + b_2 z^-2 + ... + b_nb z^(-nb). [numerator polynomail]
    -- C = 1 + c_1 z^-1 + c_2 z^-2 + ... + c_nc z^(-nc). [noise characteristics]
    -- d : delay in the system.
    A and C are monic polynomials. (in output estimation of the stochastic
    system as C is monic, we add e(t) to the estimation i.e. not starting from c1*e(t-1))


    Function inputs
    u : input to the system in column vector form
    y : input of the system in column vector form
    na : order of the denominator polynomail
    nb : order of the numerator polynomail
    nc : order of the characteristics of the noise (usually <=2 for max)
    d : number represents the delay between the input and the output
    lambda : forgetting factor -->>>     1>lambda>0
    
    Function Output
    Theta_final : final estimated parameters.
    Gz_estm : pulse (discrete) transfer function of the estimated parameters
    1 figure for the history of the parameters that are being estimated
    2 figure to validate the estimated parameters on the given output
    using the instantaneous estimated parameters.
    3 figure to plot the input versus time.

        Note: the noise added shall not to be with a magnitude close to the
        system output, it should be smaller, this is in simulation such as
        here or the algorithm will go crazy that can't distinguish between
        the main and the noisy signal (This can be measured in practical
        case finding noise to signal ratio).

    An example is added to illustrate how to use the funcrtion
%}

function [ Theta_final Gz_estm] = RecursiveExtendedLeastSquaresWEF( u, y, na, nb, d, nc, Ts, lambda )

nu = na + nb + nc + 1; % number of unkowns
N = length(y); % length of the vectors
n = max(na+1, nb+d+1); % used to avoid zeros at the begining and added 1 to avoid zero index in matlab
ne = N - n + 1; % number of equations that will used in estimation

L = zeros(1, N);
M = zeros(1, N);

%% Initialization
% Covariance matrix Initialization
alpha = 10E6;
for i = 1:(n-1)
    P{1,i} = alpha*eye(nu, nu);
end

% Initialization of estimated parameters
theta = cell(1, N);
for i = 1:(n-1)
    theta{1, i} = zeros(nu,1);
end

% initialization of the epsilon (residual)
for i = 1:n
    eps(i) =0;
end
C = 0; % as if nc = 0, as in RLS

% Initialization of cell for legend the estimated paramters
Name = cell(nu, 1);

% Initialization of y_estm
y_estm = zeros(N, 1);


%% Recursive Extended Least Squares Algorithm
% Filling the Phi matrix
for j = n : N
    
    for i = 1 : na % this for loop used to fill parts in the same line that pertains to the output "denomenator"
        if ((j-i)<=0)
            Phi(j, i) = 0;
        else
            Phi(j, i) = -y((j-i));
        end
    end
    
    
    for i = na+1 : na+nb+1 % this for loop used to fill parts in the same line that pertains to the input "numerator" starts from the last postion column of "na +1
        if ((j-d-(i-(na+1)))<=0) % add na as we left the output going to the input and we start index after "na"
            Phi(j, i) = 0;
        else
            Phi(j, i) = u(j-d-(i-(na+1)));
        end
    end
    
    
    for i = (na+nb+1)+1 : nu % this for loop used to fill parts in the same line that pertains to the noise polynomial starts from the last postion column of "na+nb +1
        if ((j-(i-(na+nb+1)))<=0) % add na+nb as we left the input going to the noise and we start index after "na"
            Phi(j, i) = 0;
        else
            Phi(j, i) = eps(j-(i-(na+nb+1)));
        end
    end
    
    % Phi = [phi(i) phi(i+1) ... phi(N)], Phi is the matrix N*nu but
    % the phi is the vector nu*1
    phi = Phi(j, :)'; % to get the j th row and then convert into column and that what we need in RLS not Phi
    
    % RELS algorithm
    P{1, j} = 1/lambda*(P{1, j-1} - P{1, j-1}*phi*inv(lambda+phi'*P{1, j-1}*phi)*phi'*P{1, j-1}); % Covariance Matrix
    K = P{1, j}*phi; % Gain
    eps(j) = y(j) - phi'*theta{1, j-1}; % residual calculation
    theta{1, j} = theta{1, j-1} + K*(eps(j)); % estimated parameters
    eps(j) = y(j) - phi'*theta{1, j-1}; % residual calculation
    
    A = 1;
    for i = 1: na
        A(i+1) = theta{1, j}(i);
    end
    for i = na+1:na+nb+1
        B(i-na) = theta{1, j}(i);
    end
    for i = (na+nb+1)+1 : nu
        C(i-(na+nb+1)) = theta{1, j}(i);
    end
    
    L(j) = outputestimation( A, B, d, u, L, j );
    M(j) = outputestimation( A, [1 C], 0, eps, M, j );
    y_estm(j) = L(j)+M(j);
end

Theta = cell2mat(theta); % converting the Total_loop_estimate from cell array to matrix, the final column is the same as the final estimate
Theta_final = Theta(:, N); % final estimation of parameters

%% Parameters and transfer function preparation (final Parameters)
z = tf('z');
A = z^(d+nb);
B = 0;
C = 1;
for i = 1:na
    a(i) = Theta_final(i);
    A = A + a(i)*z^(d-na+nb+na-i);
end
for i = na+1:na+nb+1
    b(i-na) = Theta_final(i);
    B = B + b(i-na)*z^(nb - (i-(na+1)));
end
for i = (na+nb+1)+1:nu
    c(i-(na+nb+1)) = Theta_final(i);
    C = C + c(i-(na+nb+1))*z^(-(i-(na+nb+1)));
end

Gz_estm = B/A;
[num, den] = tfdata(Gz_estm);
num = cell2mat(num);
den = cell2mat(den);
Gz_estm = tf(num, den, Ts);
C
%% Validation
final_time = (N-1)*Ts;
t = 0:Ts:final_time;

%original system
h = figure;
hold all
plot(t, y, 'LineWidth', 2)

% estimated system
plot(t, y_estm, 'LineWidth', 2);
legend(['Ordinary System'], ['Instantaneous Estimated Parameter System Output', num2str(na)]);
xlabel('Time (sec.)');
ylabel('Output');
title('Output from both Real and Estimated Systems')
grid on
print(h,'-dpng','-r500','Output')

% Parameter History Plot
for i = 1 : nu-nc
    h = figure(2);
    hold all
    plot(t, Theta(i, :), 'LineWidth', 2);
    if (i>=1 & i<=na)
        Name{i, 1} = ['a', num2str(i)]; % Polynomal A parameters
    elseif (i>=(na+1) & i<=(na+nb+1) )
        Name{i, 1} = ['b', num2str(i-na-1)]; % Polynomal B parameters
        
    end
end

grid on
xlabel('Time (sec.)')
ylabel('Estimated Parameters')
title('Estimated Parameters History')
legend(Name{1:nu-nc})
print(h,'-dpng','-r500','Estimated Parameter History');


h = figure();
for i = (na+nb+1)+1:nu
    hold all
    plot(t, Theta(i, :), 'LineWidth', 2);
    Name{i, 1} = ['c', num2str(i-(na+nb+1))]; % Polynomal B parameters
end
grid on
xlabel('Time (sec.)')
ylabel('Estimated Noise Parameters')
title('Estimated Noise Parameters History')
legend(Name{nu-nc+1:nu})
print(h,'-dpng','-r500','Estimated Noise Parameter History');

% Input plot over time
h = figure;
plot(t, u, 'LineWidth', 2);
xlabel('Time (sec.)');
ylabel('Input');
title('Input Versus Time')
grid on
print(h,'-dpng','-r500','Input')
end