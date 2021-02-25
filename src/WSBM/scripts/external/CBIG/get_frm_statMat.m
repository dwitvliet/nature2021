function outMat = get_frm_statMat(statMat, fieldName)

sizeR = size(statMat,1);
sizeC = size(statMat,2);

outMat = zeros( [ sizeR sizeC ]);

for idx=1:sizeR
    for jdx=1:sizeC
        
        if jdx < idx
           continue 
        end
        
        outMat(idx,jdx) = getfield(statMat{idx,jdx},fieldName{:});
             
    end
end
    
