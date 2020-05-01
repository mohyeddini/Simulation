function kr=rel_perm(sw,index)
nw=3.7;
no=2;
kwstar=0.6;
kostar=0.9;
swi=0.2;
sor=0.15;
s=(sw-swi)/(1-sor-swi);
if index==1
    kr=kwstar*s^nw;
elseif index==2
        kr=kostar*(1-s)^no;
        
 end
end


