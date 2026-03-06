function [n_max,str_fmt1,str_fmt2] = get_nb_char_to_print(prm)

ndigit_iter = nbdigit(prm.niter); % iteration field width (max nb of char)
ndigit_run = nbdigit(prm.nruns); % run field width

str_fmt1 = ['run %',num2str(ndigit_run),'i, iter %',num2str(ndigit_iter),...
    'i/%',num2str(ndigit_iter),'i:, K=%2i, Kopt=%2i, alpha=%d\n'];
str_fmt2 = ['run %',num2str(ndigit_run),'i, iter %',num2str(ndigit_iter),...
    'i: burnin phase %3i/100\n'];

n = [length(sprintf(str_fmt1,0,0,0,0,0,prm.alpha0)),...
    length(sprintf(str_fmt2,0,0,0))];
n_max = max(n);
n_blank = n_max-n;

str_fmt1 = [str_fmt1(1:end-2),repmat(' ',1,n_blank(1)),'\n'];
str_fmt2 = [str_fmt2(1:end-2),repmat(' ',1,n_blank(2)),'\n'];