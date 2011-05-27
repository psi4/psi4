function [S] = extractHeavisides(old)

% Transform heavisides to ternary statements to avoid NaN
% Nesting is NOT allowed

[heavi_start heavi_end] = regexp(old, 'heaviside', 'start', 'end');

% Argument first
for k = 1:length(heavi_start)
    pos = heavi_end(k);
    opens = 0;    

    started = false;
    while true
        if (old(pos) == '(')
            opens = opens + 1;
            if (~started)
                heavi_open(k) = pos;
            end
            started = true;
        end
        if (old(pos) == ')')
            opens = opens - 1;
            if (started && opens == 0)
                heavi_close(k) = pos;
                break;
            end
        end
        pos = pos + 1;
    end
end

% Left side
for k = 1:length(heavi_start)
    pos = heavi_start(k) - 1;
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
            heavi_left(k) = pos;
            break;
        end
        pos = pos - 1; 
    end
end

% Right side
for k = 1:length(heavi_start)
    pos = heavi_close(k) + 1;
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
            heavi_right(k) = pos;
            break
        end
        pos = pos + 1; 
    end
end

if (length(heavi_start) == 0)
    S = old;
else
    S = old(1:heavi_left(1));
    for k = 1:length(heavi_start)
        argument = old(heavi_open(k)+1:heavi_close(k)-1);
        left = old(heavi_left(k)+1:heavi_start(k)-1);
        right = old(heavi_close(k)+1:heavi_right(k)-1);
        S = [S ' ( (' argument ' > 0.0) ? ' left '1.0' right ' : 0.0 ) '];
        if (k ~= length(heavi_start))
            S = [S old(heavi_right(k):heavi_left(k+1))];
        end    
    end
    S = [S old(heavi_right(end):end)];
end

