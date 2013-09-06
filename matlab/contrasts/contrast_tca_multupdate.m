function [score,details]=contrast_tca_multupdate(x,W,param,is,details);
switch (param.contraste)
    case 'kde'
        [score,details]=contrast_tca_kde_multupdate(x,W,param,is,details);
    case 'kgv'
        [score,details]=contrast_tca_kgv_multupdate(x,W,param,is,details);
    case 'cum'
        [score,details]=contrast_tca_cum_multupdate(x,W,param,is,details);
    case 'statio'
        [score,details]=contrast_tca_statio_multupdate(x,W,param,is,details);
end
