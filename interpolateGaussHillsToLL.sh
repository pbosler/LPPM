#!/bin/bash

dataDir=/ascldap/users/pabosle/modelData/vtkOut

cp interpVTK2NCL.exe $dataDir/.

cd $dataDir
for file in *2__000?.vtk
do
cat <<EOF > $dataDir/InterpVTKtoNCL.namelist
&filenames
inputfile = '${file}'
outputPrefix = '${file%.vtk}'
degreeSpacing = 0.5
/

EOF

./interpVTK2NCL.exe

done