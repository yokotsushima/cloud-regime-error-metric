# Introduction #
This is a set of metrics based on RMS error which indicate the overall performance of a model to simulate the cloud radiative effects (CREs) of cloud regimes in annual mean climatology as well as in natural variations. The cloud regimes are defined as grid-boxes sharing similar daily mean cloud top pressures, optical depths and cloud covers from the ISCCP data.
It is documented in [[Williams and Webb 2009](http://rd.springer.com/article/10.1007%2Fs00382-008-0443-1)] and [[Tsushima et al.,2012](http://rd.springer.com/article/10.1007/s00382-012-1609-4)].

Results from CFMIP1 are in [[Williams and Webb 2009](http://rd.springer.com/article/10.1007%2Fs00382-008-0443-1)].

Results from CMIP5 are in [[Tsushima et al.,2012](http://rd.springer.com/article/10.1007/s00382-012-1609-4)]

## Cloud Regime Error Metric for annual mean climatology ##

This is a single scalar metric for the annual mean climatology of net CRE, a CREM for the present day (CREMpd). It is RMS error (RMSE) of Net CRE and can be broken down into two components of errors for each regimes; the error from the relative frequency of occurrence (RFO) and the error from Net CRE when the regime occurs (NCRE).

**Figure 1.** Changes between CFMIP-1 and CMIP5 of Error components in Net CRE RMSE within daily ISCCP simulator cloud regimes in the tropics. Cloud regimes are in the order of larger albedo.
![http://cfmip.metoffice.com/cremr.png](http://cfmip.metoffice.com/cremr.png)



Improvements are seen in mainly in the radiative properties of optically thicker low and high topped cloud regimes rather than frequency of occurrence.[[Tsushima et al.,2012](http://rd.springer.com/article/10.1007/s00382-012-1609-4)]

## References ##
  * Williams KD, Webb MJ (2009) A quantitative performance assessment of cloud regimes in climate models. Clim Dyn 33(1):141–157
  * Tsushima Y, Ringer MA, Webb MJ, Williams KD (2012) Quantitative evaluation of the seasonal variations in climate model cloud regimes. Clim Dyn 41(9-10):2679-2696