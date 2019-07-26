% Status of error estimates (Version 1.0, 1/Sep/01)
%
%  To generate a harmonic analysis, we
%
%     a) fit the series to either sines and cosines at specified
%        (positive) frequencies, or to complex exponentials at both
%        positive and negative frequencies. Traditionally sine/cosines
%        have been used, in t_tide I deal with the complex exponentials
%        as the math is conceptually (and practically) simpler. The
%        use of complex exponentials also unifies the treatment of
%        real time series (e.g. tidal height) and complex time
%        series (e.g. currents, u+i*v).
%     b) The complex amplitudes of constituents are corrected for
%        various factors (nodal modulations, inference, etc.)...
%     c) ...and then converted to ellipse parameters (semi-major axis,
%        Greenwich phase, etc.).
%
%  Since it is well-known that some of the less-important constituents
%  can be below the "geophysical noise" level, it is important to have 
%  confidence intervals for these estimates, giving some idea about
%  their trustworthiness for predictive purposes.
%
%  Currently two different methods of estimating confidence intervals
%  are implemented in t_tide.
%
%  1) "Linear Analysis"
%
%  Munk and Cartwright in their paper on the "response method" (Phil. 
%  Trans Roy. Soc. Lon. A, vol 259, 1966, pg 533-581) outlined the
%  analysis for estimating noise levels in spectra. Although this was 
%  used for "transfer function" estimates, the formalism was adapted by 
%  W. S. Brown and J. D. Irish (unpublished notes, 1991) for the case
%  when the transfer function is the output of a harmonic analysis. 
%  A conversion from errors in the cos/sine amplitudes to errors
%  in ellipse parameters (axis lengths, inclination, etc.) can be done
%  through a linearized analysis (B. Beardsley, unpublished notes, 1999,
%  checking earlier work by R. Signell).
%
%  In essence, 9 frequency bands are chosen, bracketing M0,M1...M8. In 
%  each band the amplitude of the residual power spectrum is estimated.
%  It is then assumed that this noise contaminates both sin and cosine
%  components of the harmonic fit equally. Errors in ellipse parameters
%  are determined through a linearized analysis in which variances are
%  summed, weighted by analytically calculated sensitivity terms.
%
%  The results appears to be adequate for real time series (e.g., 
%  tidal height), as long as the SNR (amp/error)^2 <10, and is probably 
%  not bad for SNR as low as 2 or 3.
%
%  The formalism has been extended to complex time series, with the
%  assumption that noise in the real and imaginary parts is 
%  uncorrelated, although the levels of noise in both directions are 
%  allowed to be different. Unfortunately this means that error bars 
%  are not rotationally independent, i.e., the size of error bars may 
%  change depending on whether you input your time series in N/E 
%  coordinates, or along/across bathymetry coordinates. 
%
%  I visualize the bivariate error model by seeing it as an ellipse on
%  the plane - if the ellipse is not pointing along the x or y axes,
%  this error analysis will fail to produce the correct results. Of 
%  course, you can always rotate the coordinate system so that it is 
%  aligned with the axes.
%
%  Thus I currently recommend that if you really care about your CI,
%  you submit your time series in a coordinate system in which the
%  ellipse semi-major axis lies roughly along one of the coordinate 
%  axes.
%  
%
% 2) "Nonlinear analysis"
%
% The linear analysis described above has two problems - first, it 
% linearizes a non-linear transformation, and second, it doesn't really
% account for complex noise which might be correlated.
%
% The first problem is probably not much of an issue, since it is of
% concern only at very low SNR where the results are probably not very
% useful anyway. However, it would be nice to do something about this.
% The second problem is more serious (in my opinion), since it 
% potentially affects any kind of current analysis.
%
% The modern way of dealing with nonlinearity is through resampling
% techniques. Here I use a 'parametric bootstrap'. The idea is to
% estimate an uncertainty in the complex amplitudes of the constituents.
% Assuming the original noise to be roughly white, so will the noise in
% the complex coefficients. However, there are interesting correlations
% between the real and imaginary parts of both positive and negative
% frequencies. These can be written in terms of a 4x4 matrix in which
% the variance and covariance of the real and imaginary parts of the 
% time series appear. An eigenvalue decomposition can be used to 
% generate a transformation matrix that will take 4 uncorrelated white
% noise series and give us noise of the correct characteristics; this 
% can then be use to generate a series of 'constituent replicates'.
%
% The replicates are then nonlinearly transformed into ellipse 
% parameters, and confidence intervals estimated directly from the 
% results, taking into account the nonlinearity.
% 
% The computational time required is essentially minimal using a modern
% PC.
%
% Thus, this can account for both the nonlinear transformation, and for
% correlated noise! So, what's the problem?
%
% The first is that, so far, I have done the math for just two cases:
%  a) white, correlated continuum spectrum (suitable for complex time
%     series with a flat background spectrum)
%  b) coloured, uncorrelated continuum spectrum (suitable for scalar
%     time series with a sloped background spectrum).
%
% Ideally, of course, I would like to handle coloured, correlated noise.
% I'm working on it...
%
%
% SUMMARY
%
% For scalar time series - either method is OK. The nonlinear method
% handles low SNR cases slightly better. This shows how well both 
% methods work for a relatively low noise level (for a detailed 
% explanation of the plots, see 'help t_synth'):
%
% >t_synth('nrun',40,'error','.1*colrand(SY,-1)','time',[0:24*30],'tidecon',[1 0 0 60]);
%
% And this is for higher noise levels where the transformations is more
% nonlinear.
%
% >t_synth('nrun',40,'error','20*colrand(SY,-1)','time',[0:24*30],'tidecon',[1 0 0 60]);
% 
% For vector time series. If your noise is isotropic (in space), but
% coloured (in time), both methods work fine. Again, the nonlinear 
% method handles low SNR cases slightly better:
%
% >t_synth('nrun',40,'error','.01*(colrand(SY,-1)+i*colrand(SY,-1))','time',[0:24*30]);
%
% >t_synth('nrun',20,'error','10*(colrand(SY,-1)+i*colrand(SY,-1))','time',[0:24*30]);
%
% If your noise is non-isotropic (in space) and spectrally flat (in time), 
% the nonlinear analysis wins hands down. Here is an extreme case, where the
% noise is at 45 degrees to the axes:
%
% >t_synth('nrun',20,'error','.1*(1+i)*randn(SY)','time',[0:24*30],'boota','w');
%
% However, noise is hardly ever spectrally flat in real life.
%
% If your noise is non-isotropic (in space) and spectrally coloured,
% then AS LONG AS YOU ROTATE THE TIME SERIES so that the noise in real
% and imaginary parts is uncorrelated, which is always possible for
% bivariate noise, then both methods works fine.
%
%
% R. Pawlowicz (rich@ocgy.ubc.ca)
% 1/Sep/01

help t_errors










