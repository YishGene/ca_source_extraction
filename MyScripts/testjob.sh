#!/bin/sh
matlab -nodesktop -nosplash -r "try, disp('boo'); end; quit;"


#for ii=1:10
#          system(['qsub -b y "matlab -nodisplay -r '' disp(' num2str(ii) '); exit ''"']);
#end
#  system(['qsub -b y "matlab -nodisplay -r '' script_class_noJVM; exit ''"']);

