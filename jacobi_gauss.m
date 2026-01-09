function [nodes, weights] = jacobi_gauss(n, alpha, beta)

    if n < 0
        error('需要非负的积分点数');
    end
    
    if n == 0
        nodes = [];
        weights = [];
        return;
    end

    ab = alpha + beta;
    lnmuzero = (ab + 1) * log(2) + gammaln(alpha + 1) + gammaln(beta + 1) - gammaln(ab + 2);
    

    a = zeros(n, 1);
    b = zeros(n, 1);
    

    a(1) = (beta - alpha) / (ab + 2);
    if n > 1
        i2 = 2:n;
        abi = ab + 2 * i2;
        a(i2) = (beta^2 - alpha^2) ./ ((abi - 2) .* abi);
    end
    

    if n > 1
        b(1) = sqrt(4 * (alpha + 1) * (beta + 1) / (ab + 2 + 1e-10)^2 / (ab + 3));
    end
    if n > 2
        i2 = 2:n-1;
        abi = ab + 2 * i2;
        b(i2) = sqrt(4 * i2 .* (i2 + alpha) .* (i2 + beta) .* (i2 + ab) ./ (abi.^2 - 1 + 1e-10) ./ abi.^2);
    end
    
    if n == 1
        J = a(1);
    else
        J = diag(a) + diag(b(1:n-1), 1) + diag(b(1:n-1), -1);
    end
    

    [V, D] = eig(J);
    nodes = diag(D);
    

    weights = (V(1, :).^2) * exp(lnmuzero);
    weights = weights(:);
    

    weights = abs(weights);

    [nodes, idx] = sort(nodes);
    weights = weights(idx);
end
