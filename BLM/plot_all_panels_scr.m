%% quick script to plot some itneresting versions of the figure plots

tic;
plot_BL_curves([40 20 60]/1000,[100 75 125]/1000,1./[-0.2:.02:-0.02,0.02:0.02:0.2]);
toc;

%% this set looks at the effects of rise time changes in the same cases as the exponential filters above
% with and without spatial RF

tic;
plot_BL_curves_risetime(40/1000,100/1000,[30 20 40]/1000,[40 20 60]/1000,1./[-0.2:.02:-0.02,0.02:0.02:0.2]);
toc;

% conclusions: delta function delays are only sensible with spatial RFs,
% and have the predicted dependency on the delay once you take PD-ND.

% The excitatory decay time
% has a strong influence on the tuning, while inhibitory is not predicted
% to (except in PD-ND, since it mostly affects ND tuning). This is a little
% different from what we see, though I guess we could attribute shifts to
% both the speed up of rise and of fall in the delay lines.

% once we include rise times, then risetimes also have some influence
% on tuning on the excitatory arm, as well as the inhibitory arm. always as
% expected, I guess, that inh and exc speed ups move things in different
% directions.

% given what we see in these curves, we could plot the delay only version
% of things in early intuition plots and save rise times for later? or just
% point out that to test these interactions, one would need to have
% non-global manipulations, then launch into it? If we set up the delta
% function delay, we could still use that as a foil, and point out that
% there are many aspects of these filters that matter.