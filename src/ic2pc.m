function n_pc = ic2pc(n_ic,p)

switch p.noisetype
    case {'none','N','P'}
        n_pc = p.f*(n_ic-p.b)/(p.g*p.q);
    case 'PGN'
        n_pc = (p.f*(n_ic-p.b)/p.g - p.c)/p.q;
end