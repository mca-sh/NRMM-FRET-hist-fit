function n_ic = noisedistrib(n_pc,p)
b = p.b; % bias offset to prevent negative values (ic)
q = p.q; % detector quantum efficiency (ec/pc)
g = p.g; % EM register gain
f = p.f; % analog-to-digital factor (ec/ic) s

switch p.noisetype
    case 'none'
        n_ic = q*n_pc*g/f + b;

    case 'PGN'
        % model parameters (notation as in Hirsh et. al. 2011)
        % I light intensity (pc)
        r = p.r; % readout noise (ec) s_d
        c = p.c; % spurious charge (ec) CIC
        
        % Poisson noise of photo-electrons + poisson noise of supurious charge 
        % (PC-->EC: nb. of input e- in the EM register)
        n_ie = random('poiss', q*n_pc+c); % mean=q*n_pc+c
        
        % Gamma amplification noise, composition 
        % (EC: nb. of output e- from the EM register)
        n_oe = random('gamma', n_ie, g);  % mean=(q*n_pc+c)*g
        
        % Gausian read-out noise, composition 
        % (EC-->IC)
        n_ic = random('norm', n_oe/f, r*g/f) + b; % mean=(q*n_pc+c)*g/f + b

    case 'N'
        r = p.r; % readout noise (ec) s_d

        % (PC-->IC)
        mu_ec = phtn2ele(n_pc, g/f, q);

        % (PC-->EC) 
        sig_pe = sqrt(q*n_pc); 
        
        % (EC-->IC)
        sig_ec = sqrt((r*g/f)^2 + (p.s_q^2) + (sig_pe*g/f).^2);

        % (EC-->IC)
        n_ic = random('norm', mu_ec, sig_ec) + b;

    case 'P'
        % (PC-->EC)
        n_pe = random('poiss', q*n_pc);

        % (EC-->IC)
        n_ic = n_pe*g/f + b;

end