function[K] = unwrap_K(K)
%unroll K
    dimens = size(K);
    siz = dimens(1)*dimens(2);
    K= reshape(K,siz,siz);
end