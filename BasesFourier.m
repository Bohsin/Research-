classdef BasesFourier
    properties
        num
        states
        momenta
    end
    
    methods
        function momenta=Momenta(obj,Nsites)
            momenta=cell(1,size(obj.num,2));
            for i=1:size(obj.num,2)
                momen=0:obj.num(i)-1;
                momenta{i}=(Nsites/obj.num(i))*momen;
            end
        end
        
    end
    
end