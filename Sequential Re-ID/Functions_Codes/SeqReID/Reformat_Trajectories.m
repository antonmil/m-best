function [itlf,fp_idl] = Reformat_Trajectories(Trk,idl0)

for i=1:size(Trk,2)
    itlf(i).t_start=inf; %#ok<*AGROW>
    itlf(i).t_end=0;
    itlf(i).length=0;
    itlf(i).omega=[];
    itlf(i).data=[];
    ixds=find(Trk{i}(1,:)~=0);
    if ~isempty(ixds)
        trzs=size(Trk{i},2);
        Trk{i}(:,ixds(end)+1:trzs)=[];
    end
end
kv=0;
for k=1:size(idl0,2)
    kv=kv+1;
    fp_idl(kv).rect=[];
    fp_idl(kv).xy=[];
    
    for i=1:size(Trk,2)
        kt=find(Trk{i}(3,:)==kv);
        if ~isempty(kt)
            if Trk{i}(1,kt)~=0
                Xx=idl0(Trk{i}(3,kt)).rect(Trk{i}(1,kt),1);
                Yy=idl0(Trk{i}(3,kt)).rect(Trk{i}(1,kt),2);
                Ww=idl0(Trk{i}(3,kt)).rect(Trk{i}(1,kt),3);
                Hh=idl0(Trk{i}(3,kt)).rect(Trk{i}(1,kt),4);
                itlf(i).t_start=min(kv,itlf(i).t_start);
                itlf(i).t_end=kv;
                itlf(i).omega=[itlf(i).omega 1];
                itlf(i).length=size(itlf(i).omega,2);
                itlf(i).data=[itlf(i).data [Xx+Ww/2;Yy+Hh/2]];
            elseif Trk{i}(1,kt)==0
                kln=find(Trk{i}(1,:)~=0);
                [~,ikl]=min(abs(kln-kt));
                kn=kln(ikl);
                Xx=idl0(Trk{i}(3,kn)).rect(Trk{i}(1,kn),1);
                Yy=idl0(Trk{i}(3,kn)).rect(Trk{i}(1,kn),2);
                Ww=idl0(Trk{i}(3,kn)).rect(Trk{i}(1,kn),3);
                Hh=idl0(Trk{i}(3,kn)).rect(Trk{i}(1,kn),4);
                itlf(i).t_start=min(kv,itlf(i).t_start);
                itlf(i).t_end=kv;
                itlf(i).omega=[itlf(i).omega 1];
                itlf(i).length=size(itlf(i).omega,2);
                itlf(i).data=[itlf(i).data [Xx+Ww/2;Yy+Hh/2]];
            end
            
            
        end
    end
end