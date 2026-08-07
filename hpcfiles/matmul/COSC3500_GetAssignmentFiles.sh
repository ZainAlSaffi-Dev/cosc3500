#Running this script will pull in your individual files to whatever directory this is running in.

#Pulls in the PDF spec sheet (password for PDF is your 8-digit student number)
#The -f flag forces any existing file to be overwritten
/home/groups/cosc3500/bin/cosc3500tool get -f COSC3500_Assignment.pdf

#Pulls in your blank matrixMultiply files (by default these will not replace existing files you have)
/home/groups/cosc3500/bin/cosc3500tool get matrixMultiply.cpp
/home/groups/cosc3500/bin/cosc3500tool get matrixMultiplyGPU.cu
/home/groups/cosc3500/bin/cosc3500tool get matrixMultiplyMPI.cpp

#These files are common and shared, but are pulled in as well for convenience.
#The MakeFile
cp /home/groups/cosc3500/shared/matmul/Makefile Makefile
#The directory full of the pre-made slurm scripts
cp -r /home/groups/cosc3500/shared/matmul/slurm ./