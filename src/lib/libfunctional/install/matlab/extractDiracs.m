function [S] = extractDiracs(old)

% Transform diracs to ternary statements to avoid NaN
% Nesting is NOT allowed

[dirac_start dirac_end] = regexp(old, 'dirac', 'start', 'end');

% Argument first
for k = 1:length(dirac_start)
    pos = dirac_end(k);
    opens = 0;    

    started = false;
    while true
        if (old(pos) == '(')
            opens = opens + 1;
            if (~started)
                dirac_open(k) = pos;
            end
            started = true;
        end
        if (old(pos) == ')')
            opens = opens - 1;
            if (started && opens == 0)
                dirac_close(k) = pos;
                break;
            end
        end
        pos = pos + 1;
    end
end

% Left side
for k = 1:length(dirac_start)
    pos = dirac_start(k) - 1;
    opens = 0;  
    done = false;
    operator = false;
     
    while true
        switch old(pos)
            case ')'
                opens = opens + 1;
            case '('
                opens = opens - 1;
                if (opens < 0)
                    done = true;
                end
            case '+'
                if (opens == 0)
                    done = true;
                end
            case '-'
                if (opens == 0)
                    done = true;
                end
            case '=' % assignment
                done = true;
        end
        if done
            dirac_left(k) = pos;
            break;
        end
        pos = pos - 1; 
    end
end

% Right side
for k = 1:length(dirac_start)
    pos = dirac_close(k) + 1;
    opens = 0;  
    done = false;
    operator = false;
     
    while true
        switch old(pos)
            case '('
                opens = opens + 1;
            case ')'
                opens = opens - 1;
                if (opens < 0)
                    done = true;
                end
            case '+'
                if (opens == 0)
                    done = true;
                end
            case '-'
                if (opens == 0)
                    done = true;
                end
            case ';' % EOL 
                done = true;
        end
        if done
            dirac_right(k) = pos;
            break
        end
        pos = pos + 1; 
    end
end

if (length(dirac_start) == 0)
    S = old;
else
    S = old(1:dirac_left(1));
    for k = 1:length(dirac_start)
        argument = old(dirac_open(k)+1:dirac_close(k)-1);
        left = old(dirac_left(k)+1:dirac_start(k)-1);
        right = old(dirac_close(k)+1:dirac_right(k)-1);
        S = [S ' 0.0 '];
        if (k ~= length(dirac_start))
            S = [S old(dirac_right(k):dirac_left(k+1))];
        end    
    end
    S = [S old(dirac_right(end):end)];
end

