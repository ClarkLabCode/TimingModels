function [ stimArray ] = MovingBars(params, barParam, v)
%Moving bars of a given width, repeated with a given period. 

pdBar = double(mod(v*params.t - params.x, barParam.barPeriod) < barParam.barWidth);
ndBar = double(mod(v*params.t + params.x, barParam.barPeriod) < barParam.barWidth);

stimArray = params.mask .* cat(3,...
    barParam.mlum + barParam.c * pdBar,...
    barParam.mlum + barParam.c * ndBar,...
    barParam.mlum - barParam.c * pdBar,...
    barParam.mlum - barParam.c * ndBar);

end

