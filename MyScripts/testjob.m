for ii=1:10
          system(['qsub -b y "matlab -nodisplay -r '' disp(' num2str(ii) '); exit ''"']);
end


