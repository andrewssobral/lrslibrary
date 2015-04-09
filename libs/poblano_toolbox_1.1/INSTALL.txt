To install the Poblano Toolbox for MATLAB:

1. Unpack the compressed file. In a linux environment, for example, this can
   be done from the command line via:

     unzip poblano_toolbox_1.1.zip

   *or*

     gunzip -c poblano_toolbox_1.1.tgz | tar xvf -

   This should create a directory named *poblano_toolbox_1.1*.

2. Rename the root directory from *poblano_toolbox_1.1* to *poblano_toolbox*.

3. Start MATLAB.

4. Within MATLAB, cd to the *poblano_toolbox* directory and execute the
   following commands:

     addpath(pwd) %<-- Add the Poblano toolbox to the MATLAB path
     savepath %<-- Save for future MATLAB sessions

   OR enter the following command:

     install_poblano
