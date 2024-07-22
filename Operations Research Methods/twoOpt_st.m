function [opt2,tour] = twoOpt_st(A,stour,sopt)

%% computes the 2-opt heuristic for TSP
% call: [opt,tour] = twoOpt_st(A,stour,sopt);
% input: A  - adjacency matrix of a graph
%        stour - initial tour
%        sopt - length of the initial tour
% output: opt2 - the length of the best tour
%         tour - the best 2-opt tour  
%% remark - in testa you can see the effect of replacing
%%       edges [a,b] and [c,d] with [a,c], [b,d])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tour = stour;

%% initialization
opta=sopt; %start value
testa=0;
tmin = 0;
i=0;
n=length(A);
b=tour(n);
k=0;
while i<n-2 
     a=b;
     i=i+1;
     b=tour(i);
     
     j=i+1;
     d=tour(j);
     while j<= n-1
         c=d;
         j=j+1;
         d=tour(j);
         if d~=a                     
             %to eliminate [a,b] and [c,d] or not?    
             test = A(a,c)+A(b,d)-(A(a,b)+A(c,d));                                     
             
             opta=[opta; sopt+test]; 
             testa=[testa;test];  %evaluates exchanging edges       
             
             if test < tmin             
                 tmin = test;   % the best reduction          
                 imin = i;             
                 jmin = j;         
             end                

         end
     end 
    
end

%Update tour if tmin is smaller than 0. This means that the tour remains
%unchanged if tmin = 0. 
if tmin < 0        
   tour(imin:jmin-1) = tour(jmin-1:-1:imin);       
end

opt2 = sopt + tmin; 
          