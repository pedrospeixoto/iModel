%% Analyse engen values
%%
    B=dlmread('swm_tc11_HC_vrecperhx_crecpered_sintbary_grdtrsk_linearmatrixsparse_HR95JT_004.txt')
    A=sparse(B(:,1),B(:,2), B(:,3))
    e=eigs(A, 10240)   %%sizeofmatriz -2 to be able to use sparse matrix - matlab limitation
    e=eigs(A, 1, 'lr') %% Find the largest real eigenvaleu - there should not be any
    e=eigs(A, 1, 'li') %% Find the largest imaginary eigenvalue
    e=eigs(A, numof, 'si') %% Find the largest imaginary eigenvalue
    ie=imag(e)   %%imaginary part - of interest
    re=real(e)   %% real part - must be zero to be stable
    ies=sort(ie) %% sort eigenvalues 
    plot(ies)
