function ft = time_function3(type, tm, ts, tref)
% return value of time function f(t) at tref
%
% 2011-OCT-17 Kurt Feigl
%
%    inputs:
%          tepochs - me x 1 vector of epochs in years
%          tref  - scalar reference epoch in years
%    output
%          ft      - me x 1 vector containing value of time function
%                    evaluated at each epoch
%          type == 'heaviside' Heaviside step function
%             ft(t)    = 0 if t <  tref
%             ft(t)    = 1 if t >= tref
%


if nargin ~= 4
    error(sprintf('wrong number of arguments %d. Need 4\n',nargin));
end

ft = zeros(size(tref));

switch(lower(type))
    case {'rate','secular'}
        %itime=find(isfinite(tepochs));
        %ft = tepochs-tref;
        error('A');
    case {'step','heaviside'}
        %itime=find(tepochs >= tref);
        %ft(tepochs >= tref) = 1.0;
        error('B');
    case {'pwl'}
        if numel([tm ts]) == 2
            t1 = min([tm ts]);
            t2 = max([tm ts]);
            for j=1:numel(tref)-1
                if t1 >= tref(j) && t1 <= tref(j+1)
                    if t2 >= tref(j) && t2 <= tref(j+1)
                        ft(j) = t2 - t1; % interferogram starts and ends in interval
                    elseif t2 >= tref(j+1)
                        ft(j) = tref(j+1) - t1; % interferogram starts during interval and ends after interval
                     end
                 elseif t1 <= tref(j)
                    if t2 >= tref(j) && t2 <= tref(j+1)
                        ft(j) = t2 - tref(j); % interferogram starts before interval and ends during interval
                    elseif t2 >= tref(j+1)
                        ft(j) = tref(j+1)-tref(j); % interferogram spans interval j
                    end
                end
            end
        else
            error('C');
        end
    otherwise
        warning(sprintf('undefined type %s',type));
end
return
end

