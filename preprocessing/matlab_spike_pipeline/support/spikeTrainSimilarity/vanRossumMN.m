% Calculates mulit-neuron van Rossum distance using kernels and markage vectors:
% x_spikes and y_spikes are matrices with the same number of spike trains (rows)
% rsd_tc: exponential decay (time scale parameter)
% cosalpha: cos (angle), 0: LL(labeled line),1: SP (summed population)
%
% For a detailed description of the methods please refer to:
%
% Houghton C, Kreuz T:
% On the efficient calculation of van Rossum distances.
% Network: Computation in Neural Systems, submitted (2012).
%
% Copyright:  Thomas Kreuz, Conor Houghton, Charles Dillon
%
%

function MVRD = vanRossumMN(x_spikes,y_spikes,rsd_tc,cosalpha)

x_num_trains = size(x_spikes,1);
x_num_spikes = zeros(1,x_num_trains);
for trac=1:x_num_trains
    x_num_spikes(trac) = find(x_spikes(trac,:),1,'last');
end
y_num_trains = size(y_spikes,1);
y_num_spikes = zeros(1,y_num_trains);
for trac=1:y_num_trains
    y_num_spikes(trac) = find(y_spikes(trac,:),1,'last');
end
if x_num_trains~=y_num_trains
    error('Number of trains do not match');
else
    num_trains = x_num_trains;
end

if cosalpha == 1 && rsd_tc == Inf
    
    MVRD = sqrt(sum(x_num_spikes)*(sum(x_num_spikes)-sum(y_num_spikes))+sum(y_num_spikes)*(sum(y_num_spikes)-sum(x_num_spikes)));
    
else
    if rsd_tc ~= Inf
        
        exp_x_spikes = exp(x_spikes/rsd_tc);
        exp_y_spikes = exp(y_spikes/rsd_tc);
        inv_exp_x_spikes = 1./exp_x_spikes;
        inv_exp_y_spikes = 1./exp_y_spikes;
        
        D=0;
        x_markage=ones(num_trains,max(x_num_spikes));
        y_markage=ones(num_trains,max(y_num_spikes));
        for trac=1:num_trains

            for spc=2:x_num_spikes(trac)
                x_markage(trac,spc)=1+x_markage(trac,spc-1)*exp_x_spikes(trac,spc-1)*inv_exp_x_spikes(trac,spc);
            end
            for spc=2:y_num_spikes(trac)
                y_markage(trac,spc)=1+y_markage(trac,spc-1)*exp_y_spikes(trac,spc-1)*inv_exp_y_spikes(trac,spc);
            end

            xmat=bsxfun(@rdivide,exp_x_spikes(trac,1:x_num_spikes(trac)),exp_x_spikes(trac,1:x_num_spikes(trac))');
            Dxx=x_num_spikes(trac)+2*sum(sum(tril(xmat,-1)));

            ymat=bsxfun(@rdivide,exp_y_spikes(trac,1:y_num_spikes(trac)),exp_y_spikes(trac,1:y_num_spikes(trac))');
            Dyy=y_num_spikes(trac)+2*sum(sum(tril(ymat,-1)));
            
            Dxy=f_altcor_exp2(exp_x_spikes(trac,1:x_num_spikes(trac)),exp_y_spikes(trac,1:y_num_spikes(trac)),...
                inv_exp_x_spikes(trac,1:x_num_spikes(trac)),inv_exp_y_spikes(trac,1:y_num_spikes(trac)),...
                x_markage(trac,1:x_num_spikes(trac)),y_markage(trac,1:y_num_spikes(trac)));
            
            D = D + (Dxx+Dyy)/2 - Dxy;
            
        end
        D = 2/rsd_tc * D;
        
    else                                                                   % rsd_tc = Inf --- pure rate code
        
        D = sum(x_num_spikes.*(x_num_spikes-y_num_spikes)) + ...
            sum(y_num_spikes.*(y_num_spikes-x_num_spikes));
        
    end

    if cosalpha > 0
        
        MD = 0;
        for trac1=1:num_trains-1
            for trac2=trac1+1:num_trains
                if rsd_tc ~= Inf
                    
                    MD = MD - f_altcor_exp2(exp_x_spikes(trac1,1:x_num_spikes(trac1)),exp_y_spikes(trac2,1:y_num_spikes(trac2)),...
                        inv_exp_x_spikes(trac1,1:x_num_spikes(trac1)),inv_exp_y_spikes(trac2,1:y_num_spikes(trac2)),...
                        x_markage(trac1,1:x_num_spikes(trac1)),y_markage(trac2,1:y_num_spikes(trac2)));
                    MD = MD + f_altcor_exp2(exp_x_spikes(trac1,1:x_num_spikes(trac1)),exp_x_spikes(trac2,1:x_num_spikes(trac2)),...
                        inv_exp_x_spikes(trac1,1:x_num_spikes(trac1)),inv_exp_x_spikes(trac2,1:x_num_spikes(trac2)),...
                        x_markage(trac1,1:x_num_spikes(trac1)),x_markage(trac2,1:x_num_spikes(trac2)));
                    MD = MD - f_altcor_exp2(exp_y_spikes(trac1,1:y_num_spikes(trac1)),exp_x_spikes(trac2,1:x_num_spikes(trac2)),...
                        inv_exp_y_spikes(trac1,1:y_num_spikes(trac1)),inv_exp_x_spikes(trac2,1:x_num_spikes(trac2)),...
                        y_markage(trac1,1:y_num_spikes(trac1)),x_markage(trac2,1:x_num_spikes(trac2)));
                    MD = MD + f_altcor_exp2(exp_y_spikes(trac1,1:y_num_spikes(trac1)),exp_y_spikes(trac2,1:y_num_spikes(trac2)),...
                        inv_exp_y_spikes(trac1,1:y_num_spikes(trac1)),inv_exp_y_spikes(trac2,1:y_num_spikes(trac2)),...
                        y_markage(trac1,1:y_num_spikes(trac1)),y_markage(trac2,1:y_num_spikes(trac2)));
                    
                else                                                                   % rsd_tc = Inf --- pure rate code
                              
                    MD = MD + x_num_spikes(trac1)*(x_num_spikes(trac2)-y_num_spikes(trac2)) + ...
                              x_num_spikes(trac2)*(x_num_spikes(trac1)-y_num_spikes(trac1)) + ...
                              y_num_spikes(trac1)*(y_num_spikes(trac2)-x_num_spikes(trac2)) + ...
                              y_num_spikes(trac2)*(y_num_spikes(trac1)-x_num_spikes(trac1));
                    
                end
            end
        end
        
        if rsd_tc ~= Inf
            MD = 2/rsd_tc * MD;
        end
        MVRD = sqrt(abs(D+cosalpha*MD));
        
    else
        
        MVRD = sqrt(abs(D));
        
    end
end

function Dxy = f_altcor_exp2(exp_x_spikes,exp_y_spikes,inv_exp_x_spikes,inv_exp_y_spikes,x_markage,y_markage)

x_num_spikes = length(exp_x_spikes);
y_num_spikes = length(exp_y_spikes);

Dxy=0;
for i=1:x_num_spikes
    dummy=find(exp_y_spikes<=exp_x_spikes(i),1,'last');
    if ~isempty(dummy)
        Dxy = Dxy + exp_y_spikes(dummy)*inv_exp_x_spikes(i)*y_markage(dummy);
    end
end

for i=1:y_num_spikes
    dummy=find(exp_x_spikes<exp_y_spikes(i),1,'last');
    if ~isempty(dummy)
        Dxy = Dxy + exp_x_spikes(dummy)*inv_exp_y_spikes(i)*x_markage(dummy);
    end
end
