function [sols,minerr] = tdoa_offset(D,options,offset)
%% [sols,minerr] = tdoa_offset(D,options,offset)
% Unknown offsets for 3D TDOA problems with rank constraints or linear
% techniques with minimal settings
% Solvable cases ( m receivers / n senders )
% (6/8) , (9/5), (10/5), (7/6)
% INPUT:
%     D    : TDOA measurement matrix
%  options :
%            normalize  : normalize the measurement matrix so that sum_i(Dij) = 0
%            eqsf       : templates for generate rank constraint equations
%   offset : for testing, if offset is provided, use to select solutions

% m and n
m = size(D,1);
n = size(D,2);

if m == 6 && n == 8
    dim = 14;  basis_size  = 30; mdegree  = 0; action_variable = 2;
    dim2 = 12; basis_size2 = 30; mdegree2 = 5; action_variable2 = 1;
elseif m == 9 && n == 5
    dim = 1;   basis_size  = 10; mdegree  = 0; action_variable = 1;
elseif m == 7 && n == 6
    dim = 5;   basis_size = 16;  mdegree  = 0;   action_variable =2 ; % action variable seems to be important , 2nd
    dim2 = 36; basis_size2 = 60; mdegree2 = 5; action_variable2 = 1;
elseif m == 10 && n == 5
else
    error('(6,8), (9,5) , (10,5) and (7,6) are the only working cases');
end

if ~isfield(options,'normalize');options.normalize = 0;end

% normalize the measurement matrix
if options.normalize
    Dorg = D;
    mD = mean(D);
    D = D - repmat(mD,m,1);
end

if ~( m == 9) && ~(m==10)
    
    %% rank-3 constraints on sub-determinants
    if ~isfield(options,'eqsf');
        V = eye(n);
        for i = 1:n;
            tv(i) = multipol(1,V(:,i));
        end
        one = multipol(1,zeros(n,1));
        zero =multipol(0,zeros(n,1));
        
        Dv = D;
        clear Tv;
        for i = 2:m
            for j = 2:n
                Tv(i-1,j-1) = Dv(i,j)^2 - 2*Dv(i,j)*tv(j) - ...
                    Dv(1,j)^2 + 2*Dv(1,j)*tv(j)     - ...
                    Dv(i,1)^2 + 2*Dv(i,1)*tv(1)     + ...
                    Dv(1,1)^2 - 2*Dv(1,1)*tv(1) ;
            end
        end
        %% pick out all sub-matrices of 4 by 4
        rp = nchoosek([1:m-1],4);
        cp = nchoosek([1:n-1],4);
        
        keyboard;
        
        clear eqs;
        k = 1;
        for i = 1:size(rp,1)
            for j = 1:size(cp,1)
                CC = Tv(rp(i,:),cp(j,:));
                eqs(k) = detv(CC);
                k = k+1;
            end
        end
        
    else
        % use pre-compute template to generate required equations
        eqsf = options.eqsf;
        DD = reshape(D,m*n,1);
        
        for j = 1:length(eqsf)
            
            coeff = eqsf(j).coeff;
            mm1 = eqsf(j).mm1;            
            mm = eqsf(j).mm;            
            
            ccsub = (mm(n+1:end,:));
            Dr = repmat(DD,1,size(ccsub,2));
            Drr = prod(Dr.^ccsub);
            cc = ( coeff .* Drr);
            for i = 1:size(mm1,2)
                nnn = i;
                n1 = eqsf(j).ext(i).monsid;
                ccc(i) = sum(cc(n1));
            end
            eqss = matrix2polynomials(ccc,mm1);
            eqs(j) = eqss{1};
        end
        
    end
    
    %% solve for the offset
    [eqs2] = generate_equations(eqs,mdegree);
    settings.dim = dim;
    settings.basis_size = basis_size;
    settings.action_variable = action_variable;
    try
        [sols,stats] = polysolve (eqs2,settings);
    catch
        keyboard
    end
    
elseif (m==9 && n==5)
    % linear techniques with compaction from the right D_hat = D*Cn
    D = D';
    n00 = n;
    n   = m;
    m   = n00;
    TT = [D(2:end,:).^2 - repmat(D(1,:).^2,m-1,1) ; - 2*D(2:end,:); 2*D(1,:)];
    tt = ones(n,1);
    oo = tt'/TT;
    sols = [oo(end)/sum(oo(1:m-1))   oo(m:n-1)./oo(1:m-1)]';
    settings.relative = 1;
elseif (m == 10 && n == 5)
    % Pollefeys and Nister 's 
    TT = [D.^2 -2*D];
    tt = ones(m,1);
    oo = TT\tt;
    sols = oo(6:end)./oo(1:5);
    settings.relative = 1;
end

% push back the mean
if options.normalize
    sols = sols + repmat(mD',1,size(sols,2));
end

%% TODO: for loop
if nargin == 3
    [minerr,err] = evaluate_solutions(sols,offset,settings);
    oo = sols(:,find(err==minerr));
end
if ~exist('minerr','var')
    minerr = nan;
end
