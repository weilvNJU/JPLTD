function [G, L,S, Loss] = JPLTD(X, index, gt, param, k)
%% Normalized data
num_view = length(X);
N = numel(gt);
dk=cell(1, num_view);
for i=1:num_view
    X{i} = NormalizeData(X{i});
    dk{i} = size(X{i},1);
    num{i} = size(X{i}, 2);
end
%%
O = cell(1, num_view); 
for i=1:num_view
    pos = find(index(:,i)==1);
    T= zeros(N, length(pos));
    for j=1:length(pos)
        T(pos(j),j) = 1;
    end
    O{i} = T';
end
%% parameters
alpha  = param.lambda;
beta   =  param.theta;
mu = 1;
%% initialization
W = cell(num_view,1);

for i = 1:num_view
    E{i} = zeros(k, num{i});
    L{i} = zeros(N);
    Z{i} = zeros(N);
    F1{i} = zeros(N);
    ki = size(X{i},1);
    W{i} = zeros(k,ki);
end
Y = E;
S = F1;
MAX_iter = 30;

rho = 0.001;
delta = 1.3;
%% update
for iter = 1:MAX_iter
    if mod(iter, 10)==0
        fprintf('%d..',iter);
    end

    % W_i
    for i=1:num_view
        kt = X{i}-X{i}*O{i}*Z{i}*O{i}';
        tempC = kt * (E{i} - Y{i}/rho)';
        [U,~,V] = svd(tempC,'econ');
        W{i} = V*U';
    end


    %E
    for i=1:num_view
        E{i} = prox_l21(W{i} * X{i} - W{i} * X{i}*O{i}*Z{i}*O{i}' + Y{i}/rho, mu/rho);
    end


    % Zi
    for i=1:num_view
        OiO=  O{i}'*O{i}+1e-10*eye(N);%对应求解的G'*G
        temp1 = W{i} * X{i}*O{i};%对应求解的W * Xo*G
        tempQ = W{i} * X{i}-E{i}+Y{i}/rho;%对应Q
        tempP = L{i} + S{i} -  F1{i}/rho;%对应P
        tempL = (temp1'*temp1);%第一个Z左乘的部分
        tempR = pinv(OiO);%第二个Z的右乘部分
        tempM = (tempP + temp1' * tempQ * O{i})*pinv(OiO);
        Z{i} = lyap(tempL,tempR, -tempM);
    end
    clear OiO temp1 tempL tempR tempM


    % update L
    Z_tensor = cat(3, Z{ : , : });
    F1_tensor = cat(3,F1{ : , : });
    S_tensor = cat(3, S{ : , : });
    Zv = Z_tensor(:);
    F1v = F1_tensor(:);
    Sv = S_tensor(:);
    [Lv, ~] = wshrinkObj(Zv - Sv + 1/rho*F1v, alpha/rho, [N, N, num_view], 0, 3);
    L_tensor = reshape(Lv, [N, N, num_view]);
    for i=1:num_view
        L{i} = L_tensor(:,:,i);
    end


    % S
    for i=1:num_view
        S{i} = prox_l1(Z{i}-L{i}+F1{i}/rho, beta/rho);
    end


    RR1 = []; RR2 = [];
    for i=1:num_view
        res1 = W{i} *(X{i} - X{i}*O{i}*Z{i}*O{i}') - E{i};
        res2 = Z{i}-L{i} -S{i};
        Y{i} = Y{i} + rho * res1;
        F1{i} = F1{i} + rho * res2;
        RR1 =[RR1, norm(res1,'inf')];
        RR2 =[RR2, norm(res2,'inf')];
    end

    rho = min(1e6, delta*rho);

    thrsh = 1e-5;
    if(norm(RR1, inf)<thrsh && norm(RR2, inf)<thrsh)
        break;
    end

    loss1(iter) = norm(RR1, inf); loss2(iter) = norm(RR2, inf);
    RR1 = []; RR2 = [];
end

temp=0;
for i=1:num_view
    temp = temp + (abs(L{i})+(abs(L{i}))');
end
G = temp/2/num_view;
Loss= [loss1; loss2];

end
