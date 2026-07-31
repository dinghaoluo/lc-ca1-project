% Calculates pairwise mulit-neuron van Rossum distance using kernels and markage vectors:
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

function MVRDmat = vanRossumMNPW(spikes,rsd_tc,cosalpha)

num_trials = size(spikes,1);
num_trains = size(spikes,2);
num_spikes = zeros(num_trials,num_trains);
for tric=1:num_trials
    for trac=1:num_trains
        num_spikes(tric,trac) = find(spikes(tric,trac,:),1,'last');
    end
end

if cosalpha == 1 && rsd_tc == Inf
    
    MVRDmat=zeros(num_trials);
    for tric1=1:num_trials-1
        for tric2=tric1+1:num_trials
            MVRDmat(tric1,tric2) = sqrt(sum(num_spikes(tric1,:))*(sum(num_spikes(tric1,:))-sum(num_spikes(tric2,:)))+sum(num_spikes(tric2,:))*(sum(num_spikes(tric2,:))-sum(num_spikes(tric1,:))));
        end
    end
    
else
    D=zeros(num_trials);
    if rsd_tc ~= Inf
        
        exp_spikes = exp(spikes/rsd_tc);
        inv_exp_spikes = 1./exp_spikes;
        
        markage=ones(num_trials,num_trains,max(max(num_spikes)));
        Dxx=zeros(num_trials,num_trains);
        for trac=1:num_trains
            
            for tric=1:num_trials
                for spc=2:num_spikes(tric,trac)
                    markage(tric,trac,spc)=1+markage(tric,trac,spc-1)*exp_spikes(tric,trac,spc-1)*inv_exp_spikes(tric,trac,spc);
                end
                mat=bsxfun(@rdivide,shiftdim(exp_spikes(tric,trac,1:num_spikes(tric,trac)),2)',shiftdim(exp_spikes(tric,trac,1:num_spikes(tric,trac)),2));
                Dxx(tric,trac)=num_spikes(tric,trac)+2*sum(sum(tril(mat,-1)));
            end
        
            for tric1=1:num_trials-1
                for tric2=tric1+1:num_trials
                    Dxy=f_altcor_exp2(exp_spikes(tric1,trac,1:num_spikes(tric1,trac)),exp_spikes(tric2,trac,1:num_spikes(tric2,trac)),...
                        inv_exp_spikes(tric1,trac,1:num_spikes(tric1,trac)),inv_exp_spikes(tric2,trac,1:num_spikes(tric2,trac)),...
                        markage(tric1,trac,1:num_spikes(tric1,trac)),markage(tric2,trac,1:num_spikes(tric2,trac)));
                    D(tric1,tric2) = D(tric1,tric2) + (Dxx(tric1,trac)+Dxx(tric2,trac))/2 - Dxy;
                end
            end
            
        end
        D = 2/rsd_tc * D;
        
    else                                                                   % rsd_tc = Inf --- pure rate code
        
        for tric1=1:num_trials-1
            for tric2=tric1+1:num_trials
                D(tric1,tric2) = sum(num_spikes(tric1,:).*(num_spikes(tric1,:)-num_spikes(tric2,:))) + ...
                                 sum(num_spikes(tric2,:).*(num_spikes(tric2,:)-num_spikes(tric1,:)));
            end
        end
        
    end
    
    if cosalpha > 0
        
        MD=zeros(num_trials);
        for trac1=1:num_trains-1
            for trac2=trac1+1:num_trains
                
                if rsd_tc ~= Inf
                    
                    MDxx=zeros(1,num_trials);
                    for tric=1:num_trials
                        MDxx(tric)=f_altcor_exp2(exp_spikes(tric,trac1,1:num_spikes(tric,trac1)),exp_spikes(tric,trac2,1:num_spikes(tric,trac2)),...
                                inv_exp_spikes(tric,trac1,1:num_spikes(tric,trac1)),inv_exp_spikes(tric,trac2,1:num_spikes(tric,trac2)),...
                                markage(tric,trac1,1:num_spikes(tric,trac1)),markage(tric,trac2,1:num_spikes(tric,trac2)));
                    end

                    for tric1=1:num_trials-1
                        for tric2=tric1+1:num_trials
                            MD(tric1,tric2) = MD(tric1,tric2) + MDxx(tric1) + MDxx(tric2);
                            MD(tric1,tric2) = MD(tric1,tric2) - f_altcor_exp2(exp_spikes(tric1,trac1,1:num_spikes(tric1,trac1)),exp_spikes(tric2,trac2,1:num_spikes(tric2,trac2)),...
                                inv_exp_spikes(tric1,trac1,1:num_spikes(tric1,trac1)),inv_exp_spikes(tric2,trac2,1:num_spikes(tric2,trac2)),...
                                markage(tric1,trac1,1:num_spikes(tric1,trac1)),markage(tric2,trac2,1:num_spikes(tric2,trac2)));
                            MD(tric1,tric2) = MD(tric1,tric2) - f_altcor_exp2(exp_spikes(tric2,trac1,1:num_spikes(tric2,trac1)),exp_spikes(tric1,trac2,1:num_spikes(tric1,trac2)),...
                                inv_exp_spikes(tric2,trac1,1:num_spikes(tric2,trac1)),inv_exp_spikes(tric1,trac2,1:num_spikes(tric1,trac2)),...
                                markage(tric2,trac1,1:num_spikes(tric2,trac1)),markage(tric1,trac2,1:num_spikes(tric1,trac2)));
                        end
                    end
                    
                else                                                                   % rsd_tc = Inf --- pure rate code
                    
                    for tric1=1:num_trials-1
                        for tric2=tric1+1:num_trials
                            MD(tric1,tric2) = MD(tric1,tric2) + num_spikes(tric1,trac1)*(num_spikes(tric1,trac2)-num_spikes(tric2,trac2)) + ...
                                                                num_spikes(tric1,trac2)*(num_spikes(tric1,trac1)-num_spikes(tric2,trac1)) + ...
                                                                num_spikes(tric2,trac1)*(num_spikes(tric2,trac2)-num_spikes(tric1,trac2)) + ...
                                                                num_spikes(tric2,trac2)*(num_spikes(tric2,trac1)-num_spikes(tric1,trac1));
                        end
                    end

                end

            end
        end
        
        if rsd_tc ~= Inf
            MD = 2/rsd_tc * MD;
        end
        MVRDmat = sqrt(abs(D+cosalpha*MD));

    else

        MVRDmat = sqrt(abs(D));

    end
end
MVRDmat=MVRDmat+MVRDmat';

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
