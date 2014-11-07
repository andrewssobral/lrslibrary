function [ Z, H, dnorm ] = deep_seminmf ( X, layers, varargin )

% Process optional arguments
pnames = { ...
    'z0' 'h0' 'bUpdateH' 'bUpdateLastH' 'maxiter' 'TolFun', 'verbose', 'bUpdateZ', 'cache' ...
};

num_of_layers = numel(layers);

Z = cell(1, num_of_layers);
H = cell(1, num_of_layers);

dflts  = {0, 0, 1, 1, 500, 1e-5, 1, 1, 1};

[z0, h0, bUpdateH, bUpdateLastH, maxiter, tolfun, verbose, bUpdateZ, cache] = ...
        internal.stats.parseArgs(pnames,dflts,varargin{:});

if  ~iscell(h0)
    for i_layer = 1:length(layers)
        if i_layer == 1
            % For the first layer we go linear from X to Z*H, so we use id
            V = X;
        else 
            V = H{i_layer-1};
        end
        
        if verbose
            display(sprintf('Initialising Layer #%d with k=%d with size(V)=%s...', i_layer, layers(i_layer), mat2str(size(V))));
        end
        
        % For the later layers we use nonlinearities as we go from
        % g(H_{k-1}) to Z*H_k
        [Z{i_layer}, H{i_layer}, init_err{i_layer}] = ...
             seminmf(V, ...
                 layers(i_layer), ...
                 'maxiter', maxiter, ...
                 'bUpdateH', true, 'bUpdateZ', bUpdateZ, 'verbose', verbose); 
                 %'bUpdateH', true, 'bUpdateZ', bUpdateZ, 'verbose', verbose, 'save', cache); 
    end

else
    Z=z0;
    H=h0;
    
    if verbose
        display('Skipping initialization, using provided init matrices...');
    end
end


dnorm0 = cost_function(X, Z, H);
dnorm = dnorm0;

if verbose
    display(sprintf('#%d error: %f', 0, dnorm0));
end

%%% Error Propagation
if verbose
    display('Finetuning...');
end
H_err = cell(1, num_of_layers);


for iter = 1:maxiter  
    H_err{numel(layers)} = H{numel(layers)};
    for i_layer = numel(layers)-1:-1:1
        H_err{i_layer} = Z{i_layer+1} * H_err{i_layer+1};
    end
    
    for i = 1:numel(layers)
        if bUpdateZ
            try
                if i == 1
                    Z{i} = X  * pinv(H_err{1});
                else
                    Z{i} = pinv(D') * X * pinv(H_err{i});
                end
            catch 
                display(sprintf('Convergance error %f. min Z{i}: %f. max %f', norm(Z{i}, 'fro'), min(min(Z{i})), max(max(Z{i})))); 
            end
        end
        
        if i == 1
            D = Z{1}';
        else
            D = Z{i}' * D;
        end
       
        if bUpdateH && (i < numel(layers) || (i == numel(layers) && bUpdateLastH))
            A = D * X;
            Ap = (abs(A)+A)./2;
            An = (abs(A)-A)./2;

            B = D * D';
            
            Bp = (abs(B)+B)./2;
            Bn = (abs(B)-B)./2;
    
       
            H{i} = H{i} .* sqrt((Ap + Bn * H{i}) ./ max(An + Bp * H{i}, 1e-10));
        end
    end
    
    assert(i == numel(layers));
    
    dnorm = cost_function(X, Z, H);
    
    if verbose
        display(sprintf('#%d error: %f', iter, dnorm));
    end
    
    assert(dnorm <= dnorm0 + 1, ...
        sprintf('Rec. error increasing! From %f to %f. (%d)', ...
        dnorm0, dnorm, iter) ...
    );
    
    if dnorm0-dnorm <= tolfun*max(1,dnorm0) 
        if verbose
            display( ...
                sprintf('Stopped at %d: dnorm: %f, dnorm0: %f', ...
                    iter, dnorm, dnorm0 ...
                ) ...
            );
        end
        break;
    end
    
    dnorm0 = dnorm;
end
end

function error = cost_function(X, Z, H)
    error = norm(X - reconstruction(Z, H), 'fro');
end

function [ out ] = reconstruction( Z, H )

    out = H{numel(H)};

    for k = numel(H) : -1 : 1;
        out =  Z{k} * out;
    end

end
