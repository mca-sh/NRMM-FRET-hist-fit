function n_ic = pc2ic(n_pc,p)

switch p.noisetype
    case {'none','N','P'}
        n_ic = p.q*n_pc*p.g/p.f + p.b;

    case 'PGN'
        n_ic = (p.q*n_pc+p.c)*p.g/p.f + p.b;
end