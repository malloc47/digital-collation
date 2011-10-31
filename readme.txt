Sapheos Project
Digital Colation Application 

To register pages, call SIFTRegister with the following form from a MATLAB prompt:

registered = SIFTRegister('template_file_name.tif','target_file_name.tif');

To customize the parameters, use:

registered = SIFTRegister('template_file_name.tif','target_file_name.tif',{8,10,500,3,[3,6],0.65,0,false,false});

For more information on the particular parameters involved, examine SIFTRegister.m. 

To generate an overlay so the differences between pages can be observed, use:

composite = constructOverlay(registered,base);

or, optionally:

composite = constructOverlay(registered,base,0.6);

to adjust the alpha value of the resulting composite.

Be sure to write out or display the output returned by SIFTRegister and/or constructOverlay so the results can be examined visually.

Nearest neighbor code in matlab was obtained from: http://www.mathworks.it/matlabcentral/fileexchange/19345-efficient-k-nearest-neighbor-search-using-jit
Required vl_feat package obtained from: http://www.vlfeat.org/


