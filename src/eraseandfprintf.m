function nb = eraseandfprintf(str,nb)

fprintf(repmat('\b',1,nb));
nb = fprintf(str);