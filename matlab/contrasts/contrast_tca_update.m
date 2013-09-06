function [score,details]=contrast_tca_update(x,W,param,i,details);

switch (param.type)
case 'ica'
   switch (param.contrast)
   case 'kde'
      [score,details]=contrast_ica_kde_update(x,W,param,i,details);  
   case 'kgv'
      [score,details]=contrast_ica_kgv_update(x,W,param,i,details);
   case 'cum'
      [score,details]=contrast_ica_cum_update(x,W,param,i,details);
   end
   
case 'tca'
   switch (param.contrast)
   case 'kde'
      [score,details]=contrast_tca_kde_update(x,W,param,i,details);
   case 'kgv'
      [score,details]=contrast_tca_kgv_update(x,W,param,i,details);
   case 'cum'
      [score,details]=contrast_tca_cum_update(x,W,param,i,details);
   end
end
