% t_readme.m
%
%                    Tidal Analysis Toolbox
%                              by        
%               R. Pawlowicz, R. Beardsley, S. Lentz
%
%                 Version 1.3b   December 1, 2010
%
% The tidal analysis toolbox uses harmonic analysis to estimate 
% tidal consitutents and their uncertainities in scalar and vector 
% time series. 
%
% A description of the theoretical basis of the toolbox and some
% implementation details can be found in:
%
% Pawlowicz, R., B. Beardsley, and S. Lentz, "Classical Tidal 
%   "Harmonic Analysis Including Error Estimates in MATLAB 
%    using T_TIDE", Computers and Geosciences, 28 (2002), 929-937.
%
% (citation of this article would be appreciated if you find the
%  toolbox useful).
%
% A demonstration program is included which tests many of the 
% capabilities of this toolbox. To run the demonstration, type 't_demo'.
%
% This toolbox may call functions PSD and/or CSD in the Signal 
% Processing Toolbox (see note below).
%
% ---------------------------------------------------------------------
% This package began as an attempt to translate a FORTRAN package 
% developed by M. G. G. Foreman and coworkers at the Institute of Ocean 
% Sciences (IOS). (The IOS tidal package and user manuals are available 
% from Foreman at http://www.pac.dfo-mpo.gc.ca/sci/osap/projects/tidpack/tidpack_e.htm)
% S. Lentz and R. Beardsley (WHOI) began translating the code into 
% MATLAB, and wrote a linear error estimation algorithm using MATLAB 
% spectrum code.  R. Pawlowicz (UBC) then completely rewrote it all, 
% using complex (rather than real) math, adding inference, a user 
% interface, and lots of other goodies, including an nonlinear error 
% analysis. If you want to make use of the error analysis, it is 
% strongly recommended that you read T_ERRORS.
%
% Caveats: - Presently the nodal corrections are done in such a way 
% that they may not be too accurate for time series longer than a year
% or so (corrections are based on the middle time of the input series). 
% The hand-me-down advice for this case is "break up your series into 
% yearly chunks and do successive analyses", unless you have 19 years 
% or more data, in which case nodal corrections aren't necessary, but 
% in THAT case, you should rewrite the code so that all constituents 
% are analyzed without the need for nodal corrections.
%
% Tidal predictions are also based on nodal corrections at the center 
% time (the center of the time series being analyzed).
%
% Shallow water constituents are not used automatically. They can
% be used but you must specify them manually using 'shallow' input 
% option to T_TIDE.  You can get the names and frequencies of 
% shallow water constituents from the CONST structure returned by 
% T_GETCONSTS):
%
%   CONST=t_getconsts;
%   CONST.name(isfinite(CONST.ishallow))
% 
% Note that T_TIDE has options for pretty much anything you could 
% possibly want to do -  type 'help t_tide' for more info (also look 
% at the example in t_demo).
%
%----------------------------------------------------------------------
%
% PWELCH and CPSD: Currently the functions pwelch.m and cpsd.m from the SIGNAL 
% PROCESSING toolbox are called when the default confidence interval 
% calculation is used. If you don't have these functions, you can still
% do things by specifying another algorithm that doesn't use the spectral
% estimators, e.g.,
%
%    [...]=t_tide(...'error','wboot')
%
%----------------------------------------------------------------------
%
% The toolbox presently contains the following mfiles:
%
% ---FOR ANALYSIS
% t_tide.m       - computes the tidal analysis of the real or complex 
%                  time series.
%
% t_predic.m     - computes a tidal prediction using the results of 
%                  t_tide.
%
% t_vuf.m        - computes nodal corrections. 
% 
% t_astron.m     - computes astronomical arguments.
%
% t_getconsts.m  - loads constituent data of all kinds (based on the 
%                  data file from the fortran package).
%
% t_xtide.m      - worldwide tidal predictions using the constituent data
%                  from the XTIDE package.
%
% ---FOR DOCUMENTATION
%
% t_readme.m     - this file.
%
% t_errors.m     - a long discussion of confidence intervals and how 
%                  they can be generated.
%
% ---FOR DEMONSTRATION
%
% t_synth.m      - synthesizes noisy data to test real as opposed to
%                  estimated uncertainties (crufty code).
%
% t_demo.m       - a short example using the Tuktoyuktok elevation data.
% 
% ---FOR FUN
%
% t_equilib.m    - computes the equilibrium amplitudes of main
%                  constituents at a given latitude.
%
%
% Various data files are also included:
%
% tide3.dat      - standard constituent data file from the IOS analysis
%                  package. (read once and then results stored in 
%                  various data structures in t_constituents.mat)
%
% t_equilib.dat  - equilibrium amplitude A and B factors. (read once 
%                  and results then stored in t_constituents.mat)
%
% t_constituents.mat - constituent data structures.
%
% t_example.mat  - sample Tuktoyuktuk elevation data set used as an 
%                  example in PMSR-77/10 (check example from IOS tide
%                  package).
% 
%
% Questions or comments to:
% R. Pawlowicz (rich@eos.ubc.ca) 12/1/01
% Version 1.3b

help t_readme
