# Cellular-Potts-Model
An implementation of Cell Movement, Cell Migration, and Zebrafish Tailbud Extension Using Cellular Potts Model in Matlab (click the images)

[![IMAGE ALT TEXT HERE](https://img.youtube.com/vi/t275kJ1ukXs/0.jpg)](https://youtu.be/t275kJ1ukXs)

[![IMAGE ALT TEXT HERE](https://img.youtube.com/vi/v0DvHIlUrzc/0.jpg)](https://youtu.be/v0DvHIlUrzc)

[![IMAGE ALT TEXT HERE](https://img.youtube.com/vi/X9LWbKXnyPM/0.jpg)](https://youtu.be/X9LWbKXnyPM)

1. Download the files from https://github.com/mattonics/Cellular-Potts-Model/.
2. The perim function, which takes the most time, was compiled to mex file. By default the cpm_1Cell and cpm_ManyCells files use the mex files for perim and H_J, the adhesion contribution to the Hamiltonian. Replace all instances of the mex versions of perim and H_J to perim and H_J if not using a 64 bit Windows machine.
3. Add a file path to save the generated video.
4. Run and wait a while.
