dataDir=$HOME/modelData/vtkOut

#cp interpVTK2NCL.exe $dataDir/.
cp nearestNeighborInterpToLL.exe $dataDir/.

cd $dataDir

for filename in *__000?.vtk
do
cat <<EOF > $dataDir/InterpVTKtoNCLNearestNeighbor.namelist
&filenames
inputfile = '${filename}'
outputPrefix = '${filename%.vtk}'
degreeSpacing = 0.5
/

EOF

#./interpVTK2NCL.exe
./nearestNeighborInterpToLL.exe

done