#!/bin/bash

dataDir=$HOME/modelData/vtkOut

cp $HOME/LPPM/interpVTK2NCL.exe $dataDir/.

cd $dataDir

for filename in *__000?.vtk
do
cat <<EOF > $dataDir/InterpVTKtoNCL.namelist
&filenames
inputfile = '${filename}'
outputPrefix = '$HOME/modelData/advectionPaper/${filename%.vtk}'
degreeSpacing = 0.5
/

EOF

./interpVTK2NCL.exe

done