function labdist_faster_qk_mat = labdist_faster_qkpara_opt(sa,la,sb,lb,q,k,max_size)
% LABDIST_FASTER_QKPARA(SA,LA,SB,LB,Q,K), parallel in q and k and avoiding "out of memory"-error
% Calculates the multi-unit metric distance between two spike trains
% Uses a fast version of the algorithm
% SA, SB - spike times on the two spike trains
% LA, LB - spike labels (positive integers)
% Q - vector of timing precision parameters
% K - vector of label reassigning parameters
% Max_Size - Maximum size of matrix m4 (to prevent "out of memory"-error
% (if the input parameter is omitted it is defined below)
%
% Thomas Kreuz, 10/25/08; based on code by Dmitriy Aronov, 6/20/01;

if nargin<7
    max_size=10000000;
end

%Assign labels in the form 1,2,...,L and count spikes of each label
lbs = unique([la lb]);
L = size(lbs,2);
for c = 1:L
    j = find(la==lbs(c));
    la(j) = c;
    numa(c) = size(j,2);
    j = find(lb==lbs(c));
    lb(j) = c;
    numb(c) = size(j,2);
end

%Choose the spike train to separate to subtrains
if prod(numb+1) > prod(numa+1) % corrected
    t = la;
    la = lb;
    lb = t;
    t = sa;
    sa = sb;
    sb = t;
    t = numa;
    numa=numb;
    numb = t;
    % disp('Spike trains exchanged')
end
tb=zeros(L,max(numb));
for c = 1:L
    tb(c,1:numb(c))=sb(logical(lb==c));
end

% subdivide parameter vectors (if necessary)
m2_size=prod(numb+1)*(sum(numa)+1);
if m2_size>max_size
    disp(' ')
    error('Too many spikes! Calculation impossible !!!')
else
    lkqs=[length(k) length(q)];
    m4_size=prod(lkqs)*m2_size;
    lkq=lkqs;
    paras=zeros(2,max(lkqs));
    paras(1,1:lkqs(1))=k;
    paras(2,1:lkqs(2))=q;
    if m4_size>max_size
        while m4_size>max_size
            [ml,mli]=max(lkq);
            lkq(mli)=ceil(lkq(mli)/2);
            m4_size=prod(lkq)*m2_size;
        end
        kq_runs=ceil(lkqs./lkq);
        kq_calls=zeros(2,max(kq_runs)+1);
        for kqc=1:2
            kq_calls(kqc,2:kq_runs(kqc)+1)=fix(lkqs(kqc)/kq_runs(kqc))*ones(1,kq_runs(kqc));
            kq_calls(kqc,2:mod(lkqs(kqc),kq_runs(kqc))+1)=kq_calls(kqc,2:mod(lkqs(kqc),kq_runs(kqc))+1)+1;
        end
        kqm=zeros(2,max(kq_runs),max(max(kq_calls)));
        for kqc=1:2
            for kqc2=1:kq_runs(kqc)
                kqm(kqc,kqc2,1:kq_calls(kqc,kqc2+1))=paras(kqc,sum(kq_calls(kqc,1:kqc2))+(1:kq_calls(kqc,kqc2+1)));
            end
        end
    else
        kq_runs=[1 1];
        kq_calls=[zeros(1,2)' lkqs'];
        kqm=zeros(2,1,max(max(kq_calls)));
        kqm(1,1,1:kq_calls(1,2))=k;
        kqm(2,1,1:kq_calls(2,2))=q;
    end
end

%Set up an indexing system
ind = [];
for c = 1:L
    j = repmat(0:numb(c),prod(numb(c+1:end)+1),1);
    j = repmat(reshape(j,numel(j),1),prod(numb(1:c-1)+1),1);
    ind = [ind j];
end
ind = sortrows([sum(ind,2) ind]);
ind = ind(:,2:end);

%Initialize the array
m2 = zeros(size(ind,1),size(sa,2)+1);
m2(1,:) = 0:size(sa,2);
m2(:,1) = sum(ind,2);

labdist_faster_qk_mat=zeros(length(k),length(q));
for kc=1:kq_runs(1)
    for qc=1:kq_runs(2)
        % kqc=[kc qc]
        kv=shiftdim(kqm(1,kc,1:kq_calls(1,kc+1)),2)';
        qv=shiftdim(kqm(2,qc,1:kq_calls(2,qc+1)),2)';
        clear m4
        m4 = repmat(shiftdim(m2,-2),[length(kv),length(qv),1,1]);

        %Perform the calculation
        for v = 2:size(m4,3)
            fa2=find(shiftdim(m4(1,1,:,1),2)==m4(1,1,v,1)-1);
            fa=fa2(logical(sum(ind(fa2,:)-repmat(ind(v,:),length(fa2),1)==0,2)==L-1));
            fth=find(ind(v,:)>0)';
            bsv=diag(tb(fth,ind(v,fth)));
            for w = 2:size(m4,4)
                m4(:,:,v,w)=min(cat(3,m4(:,:,v,w-1)+1,m4(:,:,fa,w)+1,m4(:,:,fa,w-1)+ ...
                    repmat(shiftdim(qv'*abs(sa(w-1)-bsv'),-1),[length(kv),1,1])+ ...
                    permute(repmat(shiftdim(kv'*not(la(w-1)==fth)',-1),[length(qv),1,1]),[2 1 3])),[],3);
            end
        end
        labdist_faster_qk_mat(sum(kq_calls(1,1:kc))+(1:kq_calls(1,kc+1)),sum(kq_calls(2,1:qc))+(1:kq_calls(2,qc+1))) = m4(:,:,end,end);
    end
end
