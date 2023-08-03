function model = FlexibilizeConstraints(model, flex, minInf, maxInf, print)
%
% Function for flexibiling upper and lower bound constraints in a GEM.
%
% model        GEM
% flex1        Flexibilization factor (Default: 0.01).
% minInf       Value that is considered as -Inf (or desired minimum cutoff value)
% maxInf       Value that is considered as +Inf (or desired maximum cutoff value)
% print        Whether new constraits should be printed or not (Default: false)
%
% Based on the code for the Cobra Toolbox function printConstraints

if nargin < 5
    print = false;
    if nargin < 4
        maxInf = 1000;
        if nargin < 3
            minInf = -1000;
            if nargin < 2
                flex = 0.01;
            end
        end
    end
end

if print == true
    fprintf('Constraints before flexibilization: \n');
    printConstraints(model, minInf, maxInf);
end

minConstraints = intersect(find(model.lb > minInf), find(model.lb));
for i = 1:length(minConstraints)
    model.lb(minConstraints(i)) = model.lb(minConstraints(i)) * (1-flex);
end

maxConstraints = intersect(find(model.ub < maxInf), find(model.ub));
for i = 1:length(maxConstraints)
    model.ub(maxConstraints(i)) = model.ub(maxConstraints(i)) * (1+flex);
end

if print == true
    fprintf('Constraints after flexibilization: \n');
    printConstraints(model, minInf, maxInf);
end
end