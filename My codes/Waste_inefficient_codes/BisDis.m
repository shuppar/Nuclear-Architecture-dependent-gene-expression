%% BisDis means finding distance by bisection. It is related to crosshair function.
%% To find distance of a point from an object boundary.
%% In fact, it uses crosshair function to start out with.
%% Anyway in this form it takes longer than if you were just to 
%% do everything without any bisection. Just use crosshair.


function [dist, c_p] = BisDis(tuple, Set, varargin)

% If you want to use bisection-based distance measurement say 
% BisDis(tuple, perim_set, 'Bis', 'y'). Else just end with 
% perim_set like BisDis(tuple, perim_set). Use bisection-based
% measurement only when perim_set is > 2000 points.

%% Reading the input correctly

% tic
p = inputParser;

if nargin < 2
    error('MATLAB:narginchk:notEnoughInputs', 'Not enough input arguments.');
else
    addParameter(p, 'Bis', 'n', ...
        @(t) (ischar(t) && ismember(t, {'n', 'y'})));
end
p.KeepUnmatched = true;

parse(p, varargin{:});
Want = p.Results;
Want = Want.Bis;

%%
a = tuple(1,1); b = tuple(1,2);
p1 = crosshair([a, b], Set);


switch Want
    case 'n'   
        
        r1 = p1(1,1); c1 = p1(1,2);
        r2 = p1(2,1); c2 = p1(2,2);
        qs = Set;

        if r1 < r2
            if c1 < c2
                qs(qs(:,1)>r2, :) = [];
                qs(qs(:,1)<r1, :) = [];
                qs(qs(:,2)<c1, :) = [];
                qs(qs(:,2)>c2, :) = [];
            else
                qs(qs(:,1)>r2, :) = [];
                qs(qs(:,1)<r1, :) = [];
                qs(qs(:,2)>c1, :) = [];
                qs(qs(2,:)<c2, :) = [];
            end
        else
            if c1 < c2
                qs(qs(:,1)<r2, :) = [];
                qs(qs(:,1)>r1, :) = [];
                qs(qs(2,:)<c1, :) = [];
                qs(qs(2,:)>c2, :) = [];
            else
                qs(qs(:,1)<r2, :) = [];
                qs(qs(:,1)>r1, :) = [];
                qs(qs(2,:)>c1, :) = [];
                qs(qs(2,:)<c2, :) = [];
            end
        end
        
        d3 = round(sqrt((a-qs(:,1)).^2 + (b-qs(:,2)).^2));
        [dist, t] = min(d3); % t is index of the entry with minimum d
        r3 = qs(t, 1); c3 = qs(t, 2);
        c_p = [r3, c3];
        
%%
    case 'y'
        
        d1 = sqrt((a-p1(1,1))^2 + (b-p1(1,2))^2);
        d2 = sqrt((a-p1(2,1))^2 + (b-p1(2,2))^2);

        if d1 < d2
            if p1(2,1) == a
                drxn = 'r';
                wch = 1;
            else
                drxn = 'c';
                wch = 2;
            end
        else
            if p1(1,1) == a
                drxn = 'r';
                wch = 3;
            else
                drxn = 'c';
                wch = 4;
            end
    
        end

        switch drxn
            case 'r'
                if wch == 1
                c = round((b+p1(2,2))/2);
                else
                c = round((b+p1(1,2))/2);
                end
                p2 = crosshair([a, c], Set);
        
            case 'c'
                if wch == 2
                    c = round((a+p1(2,1))/2);
                else
                    c = round((a+p1(1,1))/2);
                end
                p2 = crosshair([c, b], Set);
        end

    d3 = round(sqrt((a-p2(:,1)).^2 + (b-p2(:,2)).^2));
    [~, c_ind] = min(d3); [~, ma_ind] = max(d3); clear d3;
    r2 = p2(c_ind, 1); c2 = p2(c_ind, 2);

    if p1(1, 1) == p2(ma_ind, 1)
        r1 = p1(2,1); c1 = p1(2,2);
    else
        r1 = p1(1,1); c1 = p1(1,2);
    end

    qs = Set;


    if r1 < r2
        if c1 < c2
            qs(qs(:,1)>r2, :) = [];
            qs(qs(:,1)<r1, :) = [];
            qs(qs(2,:)<c1, :) = [];
            qs(qs(2,:)>c2, :) = [];
        else
            qs(qs(:,1)>r2, :) = [];
            qs(qs(:,1)<r1, :) = [];
            qs(qs(2,:)>c1, :) = [];
            qs(qs(2,:)<c2, :) = [];
        end
    else
        if c1 < c2
            qs(qs(:,1)<r2, :) = [];
            qs(qs(:,1)>r1, :) = [];
            qs(qs(2,:)<c1, :) = [];
            qs(qs(2,:)>c2, :) = [];
        else
            qs(qs(:,1)<r2, :) = [];
            qs(qs(:,1)>r1, :) = [];
            qs(qs(2,:)>c1, :) = [];
            qs(qs(2,:)<c2, :) = [];
        end
    end
        
    d3 = round(sqrt((a-qs(:,1)).^2 + (b-qs(:,2)).^2));
    [dist, t] = min(d3); % t is index of the entry with minimum d
    r3 = qs(t, 1); c3 = qs(t, 2);
    c_p = [r3, c3];

end

% toc
end