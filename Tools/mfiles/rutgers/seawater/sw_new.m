%
% SW_NEW    What's new in this version of seawater.
%
% 19 April 2006  release 3.2 (For version 5.x and onwards of Matlab)
% ********************************************************
% Corrected sign of potential vorticity in sw_bfrq.
%
% 24 November 2005  release 3.1 (For version 5.x and onwards of Matlab)
% ********************************************************
% Added sw_swvel to compute surface wave velocity.
%
% 12 December 2003  release 3.0 (For version 5.x and onwards of Matlab)
% ********************************************************
% Coverted code so that temperature is now ITS-90 throughout.
%
% 25 June 1999  release 2.0.2 (For version 5.x of Matlab)
% ********************************************************
% Coding changes to enable functions to return the same shape vector as
% the input arguments.  In previous releases, some functions returned
% column vectors for row vector input.  Also some other tidying up.
%
% 22 April 1998  release 2.0.1 (For version 5.x of Matlab)
% ********************************************************
% This version is not optimised but will run under Matlab 5.x
% sw_satAr    New routine. Solubility of Ar in seawater
% sw_satN2    New routine. Solubility of N2 in seawater
% sw_satO2    New routine. Solubility of O2 in seawater
% sw_test     Updated to include tests for above
%
% April 1998  release 1.2e (For version 4.x of Matlab)
% ************************
% sw_alpha    Fixed bug where temp used in calculations regardless of
%             whether 'temp' or 'pmpt' was passed as keyword.
%
% sw_info     Shorter version. Refer users to web pages
%             http://www.marine.csiro.au
%
% sw_ver      New routine. Returns version number of SEAWATER
%
% sw_test     New Routine. Run a test on the SEAWATER routines
%             and compare results with literature values
%
% 94/11/15 release 1.2d
% **********************
% sw_bfrq.m   Now also returns potential vorticity.
%             Thanks to Greg Johnson (gjohnson@pmel.noaa.gov)
%
% sw_gvel.m   OMEGA=7.29e-5 changed to OMEGA=7.292e-5 to be
%             consistent with sw_f.m
%
%             IMPORTANT CHANGE: The usage of the following
%             routines has changed!
%
% sw_alpha.m |    All these routines expect (S,T,P) to
% sw_beta.m  |--  be passed instead of (S,PTMP,P) as in
% sw_aonb.m  |    previous releases of seawater.
%                 Fast execution can still be obtained by passing
%                 ptmp with a string flag 'ptmp' see help.
%
% 94/10/19 release 1.2c
% **********************
% Added routine sw_new.m to inform of updates and new features.
% sw_bfrq.m   Fixed bug where LAT = [] was needed as argument when
%             no latitude values are being passed.
%             Now pass PRESSURE instead of DEPTH -> more consistent
%             though only a negligible change is answers.
%
% sw_info.m   Updated to include a registration section.
%             Noted that software is FREE.
%             Noted best email address is seawater@ml.csiro.au
%             Requests for Report also via email to library@ml.csiro.au
%
% 94/10/12 release 1.2b
% ********************
% First official release and announcement on the networks.
%

more on
help sw_new
more off
