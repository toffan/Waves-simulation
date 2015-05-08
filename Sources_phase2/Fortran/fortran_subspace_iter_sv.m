function [U d] = fortran_subspace_iter_sv(A, l, p, percentage, eps, maxit)
    
    % INPUT arguments:
    %
    % version    = the subspace iteration version to be used to
    %              compute the eigenvalues
    %             
    % A          = the input matrix
    %             
    % p          = the number of products per iteration in the block
    %              algorithm
    %             
    % m          = the size of the iteration subspace
    %
    % percentage = the percentage of the trace to be retained
    %
    % eps        = threshold value to stop the subspace iteration
    %
    % maxit      = the maximum number of iterations
    %
    %
    %
    % OUTPUT arguments:
    %
    % E          = a matrix containing the computed eigenvectors by
    %              column
    % d          = a vector containing the computed eigenvalues
    

    if((nargin ~= 6) || (nargout~=2))
        fprintf('Wrong number of input/output arguments\n\n\nUsage:             \n');  
        fprintf('    INPUT arguments:                                               \n');  
%        fprintf('                                                                   \n');
%        fprintf('    version    = the subspace iteration version to be used to      \n');
%        fprintf('                 compute the eigenvalues                           \n');
        fprintf('                                                                   \n');
        fprintf('    A          = the input matrix                                  \n');
        fprintf('                                                                   \n');
        fprintf('    p          = the number of products per iteration in the block \n');
        fprintf('                 algorithm                                         \n');
        fprintf('                                                                   \n');
        fprintf('    m          = the size of the iteration subspace                \n');
        fprintf('                                                                   \n');
        fprintf('    percentage = the percentage of the trace to be retained        \n');
        fprintf('                                                                   \n');
        fprintf('    eps        = threshold value to stop the subspace iteration    \n');
        fprintf('                                                                   \n');
        fprintf('    maxit      = the maximum number of iterations                  \n');
        fprintf('                                                                   \n');
        fprintf('                                                                   \n');
        fprintf('                                                                   \n');
        fprintf('    OUTPUT arguments:                                              \n');
        fprintf('                                                                   \n');
        fprintf('    E          = a matrix containing the computed eigenvectors by  \n');
        fprintf('                 column                                            \n');
        fprintf('    d          = a vector containing the computed eigenvalues      \n');
        U=0;
        d=0;
        return
    else
        [U, d]  = mex_subspace_iter_sv(A, l, p, percentage, eps, maxit);
    end
    
    
    
    
    
