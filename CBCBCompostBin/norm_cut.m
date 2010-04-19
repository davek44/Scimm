function ncut=norm_cut(W,i,seg)
%
%  ncut=norm_cut(W,i,seg)
%
%  Purpose
%  =======
%  Compute normalized cut value.
%
a=find(seg~=i);
b=find(seg==i);

if( length(a) > 1 & length(b) > 1 )
  %cut_a=sum(sum(W(a,b)))
  %cut_b=sum(sum(W(b,a)))
  cut=sum(sum(W(a,b)));

  asso_a=sum(sum(W(a,:)));
  asso_b=sum(sum(W(b,:)));

  %norm_cut=(cut_a/asso_a)+(cut_b/asso_b);
  ncut=(cut/asso_a)+(cut/asso_b);
else
  cut=sum(sum(W(a,b)));
  asso_a=0;
  asso_b=0;
  if( cut == 0 ) 
    ncut=cut;
  else
    ncut=inf;
  end
end
%disp( sprintf( '%2.5f/(%2.5f+%2.5f) (%d,%d)',cut,asso_a,asso_b,length(a),length(b) ) );
