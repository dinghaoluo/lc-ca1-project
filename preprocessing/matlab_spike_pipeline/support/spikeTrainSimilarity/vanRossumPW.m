% Calculates pairwise van Rossum distance using kernels and markage vectors:
% x_spikes and y_spikes are matrices with the same number of spike trains (rows)
% rsd_tc: exponential decay (time scale parameter)
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

function Dmat = vanRossumPW(spikes,rsd_tc)

num_trains = size(spikes,1);
num_spikes = zeros(1,num_trains);
for trc=1:num_trains
%     disp(num2str(trc))
    num_spikes(trc) = find(spikes(trc,:),1,'last');
end
Dmat=zeros(num_trains);

if rsd_tc ~=Inf
    
    exp_spikes = exp(spikes/rsd_tc);
    inv_exp_spikes = 1./exp_spikes;

    markage=ones(num_trains,max(num_spikes));
    Dxx=zeros(1,num_trains);
    for trac=1:num_trains
        for spc=2:num_spikes(trac)
            markage(trac,spc)=1+markage(trac,spc-1)*exp_spikes(trac,spc-1)*inv_exp_spikes(trac,spc);
        end
        xmat=bsxfun(@rdivide,exp_spikes(trac,1:num_spikes(trac)),exp_spikes(trac,1:num_spikes(trac))');
        Dxx(trac)=num_spikes(trac)+2*sum(sum(tril(xmat,-1)));
    end
    
    for trac1=1:num_trains-1
        for trac2=trac1+1:num_trains
            Dmat(trac1,trac2)=Dmat(trac1,trac2)+Dxx(trac1)+Dxx(trac2)-2*f_altcor_exp2(exp_spikes(trac1,1:num_spikes(trac1)),exp_spikes(trac2,1:num_spikes(trac2)),inv_exp_spikes(trac1,1:num_spikes(trac1)),inv_exp_spikes(trac2,1:num_spikes(trac2)),markage(trac1,1:num_spikes(trac1)),markage(trac2,1:num_spikes(trac2)));
        end
    end
    Dmat = sqrt(2/rsd_tc*Dmat);
    
else                                                                   % rsd_tc = Inf --- pure rate code
    
    for trac1=1:num_trains-1
        for trac2=trac1+1:num_trains
            Dmat(trac1,trac2) = num_spikes(trac1) * (num_spikes(trac1) - num_spikes(trac2)) + num_spikes(trac2) * (num_spikes(trac2)-num_spikes(trac1));
        end
    end
    
end

Dmat=Dmat+Dmat';

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
