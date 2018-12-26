# Cellular-Potts-Model
An implementation of Cell Movement, Cell Migration, and Zebrafish Tailbud Extension Using Cellular Potts Model in Matlab
[![Watch the Video](https://www.youtube.com/watch?v=mms8z7odN0E&feature=youtu.be)](https://youtu.be/mms8z7odN0E)

1. Download the files from https://github.com/mattonics/Cellular-Potts-Model/.
2. The perim function, which takes the most time, was compiled to mex file. By default the cpm_1Cell and cpm_ManyCells files use the mex files for perim and H_J, the adhesion contribution to the Hamiltonian. Replace all instances of the mex versions of perim and H_J to perim and H_J if not using a 64 bit Windows machine.
3. Add a file path to save the generated video.
4. Run and wait a while.
