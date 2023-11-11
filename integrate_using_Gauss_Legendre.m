function [out] = integrate_using_Gauss_Legendre(a, b, n, f)
% function that calculates integral from a to b of function f using
% Gauss_Legendre rule
% a - double
% b - double
% n - int
% f - function hendler with elementwise operations(np. "./", ".*")
%out - double
if a == inf || b == inf || a >= b  || n < 3 || mod(n,2)==0 || mod (n,1) ~= 0
    error('Wrong input !!!!!')
end
    [coef, node] = coef_nodes(n);
    out = ((b-a)/2)*sum(coef .* f((((b-a)/2)*node) + ((b+a)/2)));
%     out = 0;
%     for zs = 1:n
%         out = out + ((b-a)/2)*(coef(zs) * f(((b-a)/2)*node(zs)) + ((b+a)/2));
%     end

    function [coef, nodes] = coef_nodes(n)
    % function that calculates nodes(roots of Legendre polynomia) and their
    % coefficient
    % n - int
    % coef - array of doubles
    % nodes - array of doubles
        
        nodes = Legendre_roots(n);
        coef = Legendre_weight(nodes,n);
    
        function [out] = Legendre_poly(x, n)
        % Function which calculate the Legendre polynomial of degree n at x
        % Function will use recursive formula to find Legendre polynomia.
        % To optymize complexity of this function we will use 
        % Dynamic Programing (recursion+memoziation) to
        % have time complexity = O(n)
        % x - double
        % n - int
        % out - double
    
        % Array for memoziation
        memo = zeros(1,n+1)*nan;
        memo(1) = 1;
        memo(2) = x;
        out = Legendre(n, x);
            %recursive function
            function l = Legendre(n , x)
                if ~isnan(memo(n+1))
                    l = memo(n+1);
                else
                    l = ((2*n-1)*x*Legendre(n-1,x) - (n-1)*Legendre(n-2,x))/n;
                    memo(n+1) = l;
                end
            end
        end
    
    
        function [out] = Legendre_integral(x, n)
        % Function that calculates first order derivative of
        % Legendre polynomial of degree n at x
        % x - double
        % n - int
        % out - double 
            out = n/(x^2-1)*(x*Legendre_poly(x,n) - Legendre_poly(x,n-1));
        end
    
    
        
        function[out] = Legendre_weight(x ,n)
        % Function that calculates weights for roots of Legendre polynomial of
        % degree n at x
        % x - array of doubles
        % n - int
        % out - array of doubles
            real_n = floor(n/2);
            out = zeros(1,n);
            for k = 1:(real_n  + 2)
                out(k) = 2/((1-x(k)^2)*Legendre_integral_test(x(k),n)^2);
            end
            for k = real_n + 2:n
                out(k) = out(n + 1 - k);
            end
        end
    
        
        function [out] = Legendre_roots(n)
        % finding the roots of Legendre polynomial of degree n
        % n - int
        % out - array of doubles
    
            function [out] = Initial_guess(k, n)
            % initial guess for Newton’s method
            % k - k-th guess (int)
            % n - degree of Legendre polynomial (int)
            % out - double
                out = (1 - 1/(8*n^2) + 1/(8*n^3))*cos(pi*(4*k - 1)/(4*n + 2));
            end
    
            % we only need to find half of the roots other will be have changed
            % sign
            real_n = floor(n/2);
            out = zeros(1,n);
            for k = 1:(real_n  + 1)
                x_approx = Initial_guess(k,n);
                % Newton’s method
                x_real = x_approx - Legendre_poly(x_approx, n)/Legendre_integral(x_approx, n);
                % i is a backup measurement to prevent infinite loop
                i=0;
                % 999 is constant number which allows to go into first iteration
                % of while loop
                x_prev = 999;
                % 0.000000001 - is maximum error between itertions
                while abs(x_prev - x_real) > 0.000000001 && i < 1000000
                    i = i + 1;
                    x_prev = x_real;
                    % Newton’s method
                    x_real = x_real - Legendre_poly(x_real, n)/Legendre_integral(x_real, n);
                end
                out(k) = x_real;
            end
            % middle node will always be 0 for odd numbers of nodes
            out(real_n + 1) = 0;
            % changing sign symetrically
            for k = real_n + 2:n
                out(k) = out(n + 1 - k)*(-1);
            end
        end
    end
end