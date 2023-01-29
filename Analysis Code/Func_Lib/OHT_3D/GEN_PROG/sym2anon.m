function [anon_funcs] = sym2anon(symfunc,symvars)

sym_str = ['(', char(symfunc), ')'];
sym_str = strrep(sym_str,'*','.*');
sym_str = strrep(sym_str,'/','./');
sym_str = strrep(sym_str,'^','.^');

sym_str = ['@(', symvars, ') ', sym_str, ';'];
anon_funcs = eval(sym_str);