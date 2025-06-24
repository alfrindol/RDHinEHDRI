function shares=ShamirSecretSharing(t,n,s)
    % System parameters
    %n = 5; % Number of shares
    %t = 3; % Minimum shares to reconstruct the secret
    F_num = 2^8; % Finite field GF(256)

    % Generate and distribute shares
    shares = distribute_shares(s, n, t, F_num);

    % Display shares
    % Assume we randomly select t shares for reconstruction
    %collected_shares = shares(randperm(size(shares, 1), t), :);

    % Reconstruct the secret
    %reconstructed_secret = reconstruct_secret(collected_shares, t, F_num);

    % Display the reconstructed secret
    %fprintf('Reconstructed Secret: %d\n', reconstructed_secret);
end

function shares = distribute_shares(s, n, t, F_num)
    % Generate polynomial coefficients with the secret as the constant term
    coeffs = [randi([0 F_num-1], 1, t-1), s];

    % Generate random x coordinates for shares
    %x_coords = randperm(F_num-1, n);
    x_coords=1:n;
    % Evaluate polynomial at each x coordinate
    shares = zeros(n, 2);
    for i = 1:n
        shares(i, 1) = x_coords(i);
        shares(i, 2) = mod(polyval(coeffs, x_coords(i)), F_num);
    end
end

function secret = reconstruct_secret(shares, t, F_num)
    x = shares(:, 1);
    y = shares(:, 2);
    secret = 0;

    % Use Lagrange interpolation to reconstruct the secret
    for k = 1:t
        num = 1;
        den = 1;
        for j = 1:t
            if k ~= j
                num = mod(num * (0 - x(j)), F_num);
                den = mod(den * (x(k) - x(j)), F_num);
            end
        end
        % Here we need the multiplicative inverse of den in GF(256)
        inv_den = multiplicative_inverse(den, F_num);
        term = mod(num * inv_den, F_num) * y(k);
        secret = mod(secret + term, F_num);
    end
end

function inv = multiplicative_inverse(value, F_num)
    % Compute the multiplicative inverse of a value in GF(2^8)
    % using the Extended Euclidean algorithm
    g = gcd(value, F_num); % Assuming gcd outputs the steps
    inv = mod(g(2,1), F_num);
end
