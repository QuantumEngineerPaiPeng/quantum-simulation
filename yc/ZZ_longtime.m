function ZZ_long = ZZ_longtime( ...,
    Ut, p_interval, opers, extra ...,
)
% 2019 pi
% p_interval = [p_start, p_end], p means period, not t, p=t/tau
% average from t = t_start to t_end
% or ED
% opers={ \sum S^x, \sum S^y, ... }
% Chao Yin
% last change 5/31, all_t

%% extra input
all_t = false;
dn = 1;
cross = false;
if nargin == 4
for extra_index = 1: length(extra)
    extra_item = extra{extra_index};
    if strcmp( extra_item{1}, 'all_t')
        all_t = extra_item{2};
    end
    if strcmp( extra_item{1}, 'dn')
        dn = extra_item{2};
    end
    if strcmp( extra_item{1}, 'cross')
        cross = extra_item{2};
    end
end
end % if (nargin==4)

%% start
oper_N = length(opers);
if cross
    oper_sum = zeros(oper_N);
else
    oper_sum = zeros(oper_N, 1);
end

if strcmp(p_interval, 'ED') | p_interval<0
    %% ED
    Ut.diagonalize(true); % if diagonalized, dont diagonalize again
    blockN = length(Ut.eigsys);
    for block_index = 1: blockN
        thisV = Ut.eigsys{block_index}.V;
        for line_index = 1: size(thisV, 1)
            if cross
                for oper_index = 1: oper_N
                    for oper_index1 = 1: oper_N
                        oper_sum(oper_index, oper_index1) = oper_sum(oper_index, oper_index1) + ( thisV(:, line_index)' ...
                            * opers{oper_index}.matrix{block_index}* thisV(:, line_index) ) * conj( thisV(:, line_index)' ...
                            * opers{oper_index1}.matrix{block_index}* thisV(:, line_index) );
                    end
                end
            else
                %'ED without cross'
                for oper_index = 1: oper_N
                    oper_sum(oper_index) = oper_sum(oper_index) + ( thisV(:, line_index)' ...
                        * opers{oper_index}.matrix{block_index}* thisV(:, line_index) )^2;
                end
            end % if cross
        end
    end
else
    %% dynamical
    Ut_t = Ut^ (p_interval(1)-1) ;
    Ut = Ut^dn;
    t_N = p_interval(2) - p_interval(1)+1;
    oper_sum = zeros(oper_N, floor((t_N-1)/dn)+1 );
    for t = 1:dn:t_N
        Ut_t = Ut* Ut_t;
        for oper_index = 1: oper_N
            oper_sum(oper_index, (t-1)/dn+1) = trace(Ut_t' * opers{oper_index} * Ut_t * opers{oper_index});
        end
    end % for t = 1:t_N
    if ~all_t
        oper_sum = mean( oper_sum, 2);
    end
    
end % if isstr(T_N) && strcmp(T_N, 'ED')

ZZ_long = real( oper_sum) / 2^Ut.L ;

end

