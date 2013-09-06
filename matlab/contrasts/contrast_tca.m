function [score,details]=contrast_tca(x,W,param);
switch (param.contraste)
    case 'kde'
        [score,details]=contrast_tca_kde(x,W,param);
    case 'kgv'
        [score,details]=contrast_tca_kgv(x,W,param);
    case 'cum'
        [score,details]=contrast_tca_cum(x,W,param);
    case 'statio'
        [score,details]=contrast_tca_statio(x,W,param);
end
