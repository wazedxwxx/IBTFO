
data=out(:,3);
AAA=reshape(data,[54 54]);
surf(AAA)
diff_id=AAA(1:end-1,:)-AAA(2:end,:);
XXX=reshape(out(:,4),[54 54]);
YYY=reshape(out(:,5),[54 54]);

XX= XYpoint(:,1);
YY= XYpoint(:,2);
cell_type=reshape(cell_type,[54 54]);
% pcolor(cell_type)
%  sqrt((0.59-XX).)
idx=find(cell_type>0);
% scatter(XXX(idx),YYY(idx),2)
% hold on
% scatter(XX,YY,2)