function shares = ShamirSecretSharing(t, n, s)
    % t: threshold
    % n: number of shares
    % s: secret (integer in 0~255)
    
    m = 8; % GF(2^8)
    field = gf(0, m); % initialize field
    prim_poly = field.prim_poly; % use default irreducible polynomial

    % Convert secret to GF element
    s_gf = gf(s, m, prim_poly);

    % Generate random polynomial: a(x) = a_{t-1}x^{t-1} + ... + a_1 x + s
    coeffs = [gf(randi([0 255], 1, t-1), m, prim_poly), s_gf];

    % Generate x values (must be non-zero and unique in GF)
    x_coords = gf(1:n, m, prim_poly);  % x = 1 to n in GF(256)

    % Evaluate the polynomial at each x
    shares = zeros(n, 2);  % store in (x, y) form
    for i = 1:n
        x = x_coords(i);
        y = polyval(coeffs, x);  % evaluate in GF
        shares(i, :) = [double(x.x), double(y.x)]; % convert GF back to decimal
    end
end