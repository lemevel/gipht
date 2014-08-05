function ft = time_function4(type, tm, ts, tbreaks)
% return value of time function f(t) at tbreaks
%
% 2011-OCT-17 Kurt Feigl
%
%    inputs:
%          tepochs - me x 1 vector of epochs in years
%          tbreaks  - scalar reference epoch in years
%    output
%          ft      - me x 1 vector containing value of time function
%                    evaluated at each epoch
%          type == 'heaviside' Heaviside step function
%             ft(t)    = 0 if t <  tbreaks
%             ft(t)    = 1 if t >= tbreaks
%


if nargin ~= 4
    error(sprintf('wrong number of arguments %d. Need 4\n',nargin));
end

%ft = zeros(size(tbreaks));

switch(lower(type))
    case {'rate','secular'}
        if numel(tbreaks) > 1
            warning('ignoring tbreaks in secular parameterization');
        end
        mparam = 1;
        ft = zeros(mparam,1);
        if numel([tm ts]) == 2 
            t1 = min([tm ts]);
            t2 = max([tm ts]);
            ft = t2 - t1;
        else
            error('A');
        end
   case {'step','heaviside'}
       mparam = numel(tbreaks); % number of breaks
       ft = zeros(mparam,1);
        if numel([tm ts]) == 2
            t1 = min([tm ts]);
            t2 = max([tm ts]);
            for j=1:numel(tbreaks)
                if t1 < tbreaks(j) && t2 >= tbreaks(j)
                    ft(j) = 1.0;
                else
                    ft(j) = 0.0;
                end
            end
         else
            error('B');
        end       
    case {'pwl'}
       mparam = numel(tbreaks)-1; % number of intervals
       ft = zeros(mparam,1);
       if numel([tm ts]) == 2
            t1 = min([tm ts]);
            t2 = max([tm ts]);
            for j=1:numel(tbreaks)-1
                if t1 >= tbreaks(j) && t1 <= tbreaks(j+1)  % interferogram starts during interval
                    if t2 >= tbreaks(j) && t2 <= tbreaks(j+1)   % and ends during interval
                        ft(j) = t2 - t1;                      
                    elseif t2 >= tbreaks(j+1)                   % and ends after interval
                        ft(j) = tbreaks(j+1) - t1;             
                     end
                 elseif t1 <= tbreaks(j)                  % interferogram starts before interval
                    if t2 >= tbreaks(j) && t2 <= tbreaks(j+1)    % and ends during interval
                        ft(j) = t2 - tbreaks(j);                
                    elseif t2 >= tbreaks(j+1)                    % and ends after interval
                        ft(j) = tbreaks(j+1)-tbreaks(j);         
                    end
                end
            end
        else
            error('C');
       end
    case 'step-pwl'
        mparam = 2*numel(tbreaks)-3; % number of intervals + number of breaks
        ft = zeros(mparam,1);
        if numel([tm ts]) == 2
            t1 = min([tm ts]);
            t2 = max([tm ts]);            
            fms = time_function4('pwl',t1,t2,tbreaks); % call this function to avoid repeating code above
            for j=1:numel(fms)
                ft(j) = fms(j);
            end
            j0 = numel(fms);
            fms = time_function4('step',t1,t2,tbreaks(2:end-1)); % call this function to avoid repeating code above
            for j=1:numel(fms)
                ft(j0+j) = fms(j);
            end
        else
            error('D');
        end
    case 'step-sec'
        mparam = numel(tbreaks)-2; % number number of breaks
        ft = zeros(mparam,1);
        if numel([tm ts]) == 2
            t1 = min([tm ts]);
            t2 = max([tm ts]);            
            fms = time_function4('secular',t1,t2,[]); % call this function to avoid repeating code above
            for j=1:numel(fms)
                ft(j) = fms(j);
            end
            j0 = numel(fms);
            fms = time_function4('step',t1,t2,tbreaks(2:end-1)); % call this function to avoid repeating code above
            for j=1:numel(fms)
                ft(j0+j) = fms(j);
            end
        else
                  error('E');
        end
    otherwise
        warning(sprintf('undefined type %s',type));
end
return
end

