function [x, r_norms] = tgcr(M,b,tol,maxiters,x0)
% Generalized conjugate residual method for solving Mx = b
% INPUTS
% M - matrix
% b - right hand side
% tol - convergence tolerance, terminate on norm(b - Mx) < tol * norm(b)
% maxiters - maximum number of iterations before giving up
% OUTPUTS
% x - computed solution, returns null if no convergence
% r_norms - the scaled norm of the residual at each iteration (r_norms(1) = 1)

if norm(b)<tol
    x=zeros(size(b));
    r_norms=1;
else
    % Generate the initial guess for x (zero)
    x = x0;
%     x=zeros(size(b));
    % Set the initial residual to b - Mx^0 = b
    r = b;
    
    % Determine the norm of the initial residual
    r_norms(1) = norm(r,2);
    
    for iter = 1:maxiters
        % Use the residual as the first guess for the new
        % search direction and multiply by M
        p(:,iter) = r;
        Mp(:, iter) = M * p(:,iter);
        
        % Make the new Mp vector orthogonal to the previous Mp vectors,
        % and the p vectors M^TM orthogonal to the previous p vectors
        for j=1:iter-1
            beta = Mp(:,iter)' * Mp(:,j);
            p(:,iter) = p(:,iter) - beta * p(:,j);
            Mp(:,iter) = Mp(:,iter) - beta * Mp(:,j);
        end
        
        % Make the orthogonal Mp vector of unit length, and scale the
        % p vector so that M * p  is of unit length
        norm_Mp = norm(Mp(:,iter),2);
        Mp(:,iter) = Mp(:,iter)/norm_Mp;
        p(:,iter) = p(:,iter)/norm_Mp;
        
        % Determine the optimal amount to change x in the p direction
        % by projecting r onto Mp
        alpha = r' * Mp(:,iter);
        
        % Update x and r
        x = x + alpha * p(:,iter);
        r = r - alpha * Mp(:,iter);
        
        % Save the norm of r
        r_norms(iter+1) = norm(r,2);
        
        % Print the norm during the iteration
        % fprintf('||r||=%g i=%d\n', norms(iter+1), iter+1);
        
        % Check convergence.
        if r_norms(iter+1) < (tol * r_norms(1))
%             fprintf(num2str(iter))
            break;
        end
    end
    
    % Notify user of convergence
    if r_norms(iter+1) > (tol * r_norms(1))
        fprintf(1, 'GCR NONCONVERGENCE!!!\n');
        x = [];
    else
        %fprintf(1, 'GCR converged in %d iterations\n', iter);
    end
    
    % Scale the r_norms with respect to the initial residual norm
    r_norms = r_norms / r_norms(1);
end