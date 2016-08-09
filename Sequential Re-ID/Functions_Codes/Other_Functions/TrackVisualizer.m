function TrackVisualizer(I,Trk,idl0,kv,colorord)

imshow(I),
hold on

for i=1:size(Trk,2)
    kt=find(Trk{i}(3,:)==kv);
    if ~isempty(kt)
        if Trk{i}(1,kt)~=0
            Xx=idl0(Trk{i}(3,kt)).rect(Trk{i}(1,kt),1);
            Yy=idl0(Trk{i}(3,kt)).rect(Trk{i}(1,kt),2);
            Ww=idl0(Trk{i}(3,kt)).rect(Trk{i}(1,kt),3);
            Hh=idl0(Trk{i}(3,kt)).rect(Trk{i}(1,kt),4);
        elseif Trk{i}(1,kt)==0
            kln=find(Trk{i}(1,:)~=0);
            [~,ikl]=min(abs(kln-kt));
            kn=kln(ikl);
            Xx=idl0(Trk{i}(3,kn)).rect(Trk{i}(1,kn),1);
            Yy=idl0(Trk{i}(3,kn)).rect(Trk{i}(1,kn),2);
            Ww=idl0(Trk{i}(3,kn)).rect(Trk{i}(1,kn),3);
            Hh=idl0(Trk{i}(3,kn)).rect(Trk{i}(1,kn),4);
        end
        
        rectangle('Position',[Xx,Yy,Ww,Hh],'EdgeColor',colorord(i,:),...
            'LineWidth',2)
        text(Xx+Ww/2,Yy+Hh/3,num2str(i),'Color',colorord(i,:),'FontSize',11)
    end
end

end