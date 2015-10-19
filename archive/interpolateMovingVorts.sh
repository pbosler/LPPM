#!/bin/bash

#dataDir=$HOME/modelData/vtkOut
dataDir=/fscratch/pabosle/vtkOut

cp $HOME/LPPM/interpVTK2NCL.exe $dataDir/.

cd $dataDir

for filename in *mv_conv_direct_triUnif_?__0002.vtk
do
cat <<EOF > $dataDir/InterpVTKtoNCL.namelist
&filenames
inputfile = '${filename}'
outputPrefix = '../${filename%.vtk}'
degreeSpacing = 0.5
/

EOF

./interpVTK2NCL.exe

done