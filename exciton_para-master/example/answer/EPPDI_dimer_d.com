 %RWF=EPPDI_dimer_d
 %NoSave
 %Chk=EPPDI_dimer_d
 %NProcShared=24
 %mem=20GB
 #P SP PBE1PBE/Def2SV IOP(3/33=1) IOP(6/7=3) IOP(5/33=3) IOP(3/32=2)  Nosymm

 EPPDI_dimer_d

 0 1
   O       2.268088   -5.646484   -0.504904                                                                                                                                                              
   O      -2.265045   -5.655290   -0.411633                                                                                                                                                              
   N       0.000213   -5.656939   -0.522290                                                                                                                                                              
   C       1.245225   -0.730988   -0.087953                                                                                                                                                              
   C       0.000235   -1.422842   -0.122618                                                                                                                                                              
   C      -1.244943   -0.735814   -0.036708                                                                                                                                                              
   C      -2.413905   -1.473149   -0.072226                                                                                                                                                              
   C      -2.404173   -2.864554   -0.190545                                                                                                                                                              
   C      -1.216781   -3.544268   -0.277862                                                                                                                                                              
   C       0.000446   -2.834931   -0.245045                                                                                                                                                              
   C       1.218029   -3.539547   -0.327963                                                                                                                                                              
   C       2.405368   -2.855237   -0.289509                                                                                                                                                              
   C       2.414583   -1.463795   -0.171583                                                                                                                                                              
   C       1.235882   -5.015117   -0.455915                                                                                                                                                              
   C      -1.234157   -5.019908   -0.405091                                                                                                                                                              
   C       0.000017   -7.116557   -0.668317                                                                                                                                                              
   C       0.028798   -7.831765    0.683931                                                                                                                                                              
   C       0.027682   -9.328293    0.509005                                                                                                                                                              
   C       1.223680  -10.030516    0.383094                                                                                                                                                              
   C       1.224158  -11.403934    0.185709                                                                                                                                                              
   C       0.024015  -12.096735    0.109602                                                                                                                                                              
   C      -1.174340  -11.407338    0.231322                                                                                                                                                              
   C      -1.170255  -10.033929    0.428639                                                                                                                                                              
   H      -3.374275   -0.980954   -0.007308                                                                                                                                                              
   H      -3.329705   -3.426159   -0.215534                                                                                                                                                              
   H       3.331252   -3.413263   -0.352588                                                                                                                                                              
   H       3.374910   -0.967898   -0.146176                                                                                                                                                              
   H       0.875231   -7.383519   -1.255540                                                                                                                                                              
   H      -0.897283   -7.386880   -1.219574                                                                                                                                                              
   H      -0.841594   -7.520468    1.266066                                                                                                                                                              
   H       0.921153   -7.517564    1.230143                                                                                                                                                              
   H       2.165261   -9.493283    0.444495                                                                                                                                                              
   H       2.165292  -11.935303    0.094569                                                                                                                                                              
   H       0.022667  -13.170564   -0.041518                                                                                                                                                              
   H      -2.116745  -11.941386    0.176003                                                                                                                                                              
   H      -2.110335   -9.499355    0.525859                                                                                                                                                              
   O      -2.268088    5.646484    0.504904                                                                                                                                                              
   O       2.265045    5.655290    0.411633                                                                                                                                                              
   N      -0.000213    5.656939    0.522290                                                                                                                                                              
   C      -1.245224    0.730988    0.087952                                                                                                                                                              
   C      -0.000235    1.422842    0.122618                                                                                                                                                              
   C       1.244943    0.735814    0.036708                                                                                                                                                              
   C       2.413905    1.473149    0.072226                                                                                                                                                              
   C       2.404174    2.864554    0.190545                                                                                                                                                              
   C       1.216782    3.544268    0.277862                                                                                                                                                              
   C      -0.000446    2.834931    0.245045                                                                                                                                                              
   C      -1.218029    3.539547    0.327963                                                                                                                                                              
   C      -2.405368    2.855237    0.289509                                                                                                                                                              
   C      -2.414583    1.463795    0.171583                                                                                                                                                              
   C      -1.235882    5.015117    0.455915                                                                                                                                                              
   C       1.234157    5.019908    0.405090                                                                                                                                                              
   C      -0.000017    7.116557    0.668317                                                                                                                                                              
   C      -0.028798    7.831765   -0.683931                                                                                                                                                              
   C      -0.027682    9.328293   -0.509005                                                                                                                                                              
   C      -1.223680   10.030516   -0.383093                                                                                                                                                              
   C      -1.224158   11.403934   -0.185708                                                                                                                                                              
   C      -0.024015   12.096735   -0.109602                                                                                                                                                              
   C       1.174340   11.407338   -0.231322                                                                                                                                                              
   C       1.170254   10.033929   -0.428639                                                                                                                                                              
   H       3.374275    0.980954    0.007308                                                                                                                                                              
   H       3.329705    3.426159    0.215534                                                                                                                                                              
   H      -3.331252    3.413263    0.352588                                                                                                                                                              
   H      -3.374910    0.967898    0.146176                                                                                                                                                              
   H      -0.875231    7.383519    1.255540                                                                                                                                                              
   H       0.897283    7.386880    1.219574                                                                                                                                                              
   H       0.841594    7.520468   -1.266067                                                                                                                                                              
   H      -0.921153    7.517564   -1.230143                                                                                                                                                              
   H      -2.165261    9.493283   -0.444494                                                                                                                                                              
   H      -2.165292   11.935303   -0.094569                                                                                                                                                              
   H      -0.022667   13.170564    0.041518                                                                                                                                                              
   H       2.116745   11.941386   -0.176003                                                                                                                                                              
   H       2.110335    9.499355   -0.52585                                                                                                                                                               
    O       2.268088      -5.646484      2.775096                                                                                                                                                        
    O      -2.265045      -5.655290      2.868367                                                                                                                                                        
    N       0.000213      -5.656939      2.757710                                                                                                                                                        
    C       1.245225      -0.730988      3.192047                                                                                                                                                        
    C       0.000235      -1.422842      3.157382                                                                                                                                                        
    C      -1.244943      -0.735814      3.243292                                                                                                                                                        
    C      -2.413905      -1.473149      3.207774                                                                                                                                                        
    C      -2.404173      -2.864554      3.089455                                                                                                                                                        
    C      -1.216781      -3.544268       3.002138                                                                                                                                                       
    C       0.000446      -2.834931      3.034955                                                                                                                                                        
    C       1.218029      -3.539547       2.952037                                                                                                                                                       
    C       2.405368      -2.855237      2.990491                                                                                                                                                        
    C       2.414583      -1.463795      3.108417                                                                                                                                                        
    C       1.235882      -5.015117      2.824085                                                                                                                                                        
    C       -1.234157     -5.019908      2.874909                                                                                                                                                        
    C       0.000017      -7.116557      2.611683                                                                                                                                                        
    C       0.028798      -7.831765      3.963931                                                                                                                                                        
    C       0.027682      -9.328293      3.789005                                                                                                                                                        
    C       1.223680     -10.030516      3.663094                                                                                                                                                        
    C       1.224158     -11.403934      3.465709                                                                                                                                                        
    C       0.024015     -12.096735      3.389602                                                                                                                                                        
    C      -1.174340     -11.407338      3.511322                                                                                                                                                        
    C      -1.170255     -10.033929      3.708639                                                                                                                                                        
    H      -3.374275      -0.980954      3.272692                                                                                                                                                        
    H      -3.329705      -3.426159      3.064466                                                                                                                                                        
    H       3.331252      -3.413263      2.927412                                                                                                                                                        
    H       3.374910      -0.967898      3.133824                                                                                                                                                        
    H       0.875231      -7.383519      2.024460                                                                                                                                                        
    H       -0.897283     -7.386880      2.060426                                                                                                                                                        
    H       -0.841594     -7.520468      4.546066                                                                                                                                                        
    H       0.921153      -7.517564      4.510143                                                                                                                                                        
    H       2.165261      -9.493283      3.724495                                                                                                                                                        
    H       2.165292     -11.935303      3.374569                                                                                                                                                        
    H       0.022667     -13.170564      3.238482                                                                                                                                                        
    H      -2.116745     -11.941386      3.456003                                                                                                                                                        
    H      -2.110335      -9.499355      3.805859                                                                                                                                                        
    O      -2.268088       5.646484      3.784904                                                                                                                                                        
    O       2.265045       5.655290      3.691633                                                                                                                                                        
    N      -0.000213       5.656939      3.802290                                                                                                                                                        
    C      -1.245224       0.730988      3.367952                                                                                                                                                        
    C      -0.000235       1.422842      3.402618                                                                                                                                                        
    C       1.244943       0.735814      3.316708                                                                                                                                                        
    C       2.413905       1.473149      3.352226                                                                                                                                                        
    C       2.404174       2.864554      3.470545                                                                                                                                                        
    C       1.216782       3.544268      3.557862                                                                                                                                                        
    C      -0.000446       2.834931      3.525045                                                                                                                                                        
    C      -1.218029       3.539547      3.607963                                                                                                                                                        
    C      -2.405368       2.855237      3.569509                                                                                                                                                        
    C      -2.414583       1.463795      3.451583                                                                                                                                                        
    C      -1.235882       5.015117      3.735915                                                                                                                                                        
    C       1.234157       5.019908      3.685090                                                                                                                                                        
    C      -0.000017       7.116557       3.948317                                                                                                                                                       
    C      -0.028798       7.831765       2.596069                                                                                                                                                       
    C      -0.027682       9.328293       2.770995                                                                                                                                                       
    C      -1.223680      10.030516       2.896907                                                                                                                                                       
    C      -1.224158      11.403934       3.094292                                                                                                                                                       
    C      -0.024015      12.096735       3.170398                                                                                                                                                       
    C       1.174340      11.407338       3.048678                                                                                                                                                       
    C       1.170254      10.033929       2.851361                                                                                                                                                       
    H       3.374275       0.980954      3.287308                                                                                                                                                        
    H       3.329705       3.426159      3.495534                                                                                                                                                        
    H      -3.331252       3.413263      3.632588                                                                                                                                                        
    H      -3.374910       0.967898      3.426176                                                                                                                                                        
    H      -0.875231       7.383519       4.535540                                                                                                                                                       
    H       0.897283       7.386880       4.499574                                                                                                                                                       
    H       0.841594       7.520468       2.013933                                                                                                                                                       
    H      -0.921153       7.517564       2.049857                                                                                                                                                       
    H      -2.165261       9.493283       2.835506                                                                                                                                                       
    H      -2.165292      11.935303       3.185431                                                                                                                                                       
    H      -0.022667      13.170564       3.321518                                                                                                                                                       
    H       2.116745      11.941386       3.103997                                                                                                                                                       
    H       2.110335       9.499355       2.754150                                                                                                                                                       
                                                                                                                                                                                                         
                                                                                                                                                                                                         

