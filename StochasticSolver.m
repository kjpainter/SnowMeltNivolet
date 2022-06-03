function [output measure] = StochasticSolver(rho,alpha,beta,gamma,pp,qq,rr)

global mask mask2 incidence elevation n m cells nbors data prop frac means

s = ones(m,n).*mask;
[x y z] = find(s==1);

b = zeros(m,n)./nbors; 
output = zeros(m,n,length(prop));
    
steps = 0; 
datastep = 1; 
snowcells = cells; 
measure = ones(2*length(prop),1); 

carryon = true;

% elevation and incidence angle is fixed in time
scaling = 1.0/((1.0+alpha*(means(1)^pp))*(1.0+beta*(means(2)^qq)));
fixed = scaling*((1+alpha*(incidence.^pp)).*(1+beta*(elevation.^qq)));

while carryon
% randomly order remaining snow-covered patches 
    [xx yy ii] = find(mask>0 & s >0);
    seq = randperm(length(xx));         

% test each snow-covered patch for the melt process    
    for j = 1:length(seq)
        x =xx(seq(j));
        y =yy(seq(j));
        ffunc = 1/(fixed(x,y)*(1+gamma*(b(x,y)^rr)));
        if (rand < exp(-rho*ffunc))
            s(x,y) = 0; 
            snowcells = snowcells-1; 
            if (y<n)
                b(x,y+1) = b(x,y+1)+1/nbors(x,y+1);
            end
            if (x<m)
                b(x+1,y) = b(x+1,y)+1/nbors(x+1,y);
            end
            if (y>1)
                b(x,y-1) = b(x,y-1)+1/nbors(x,y-1);
            end
            if (x>1)
                b(x-1,y) = b(x-1,y)+1/nbors(x-1,y);
            end            
        end

% calculate the measures for the error formula, if melt corresponds to those in snow cover masks        
        if (snowcells<=prop(datastep))
            output(:,:,datastep) = s;
            measure(2*datastep-1) = sum(sum(abs(s-data(:,:,datastep))))/(cells); 
            l1 = find(abs(diff(s.*mask2))==1); 
            l2 = find(abs(diff((s.*mask2)'))==1); 
            measure(2*datastep)=abs(frac(datastep)-(length(l1)+length(l2))/cells);            

            datastep=datastep+1; %          
            if datastep > length(prop) 
                carryon = false;
                break
            end
        end           
    end

% give up simulation if snow melt is too slow for melt in sensible computational timescales    
    steps = steps+1; 
    if (mod(steps,10000)==0) 
        snowcover = snowcells/cells;
        if (snowcover>1-steps/100000) 
            measure = ones(2*length(prop),1);
            carryon = false;
            break
        end
    end
  
end


      