#include "vtkActor.h"
#include "vtkCamera.h"
#include "vtkPolyDataReader.h"
#include "vtkPolyData.h"
#include "vtkPolyDataMapper.h"
#include "vtkProperty.h"
#include "vtkProperty2D.h"
#include "vtkCaptionActor2D.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkScalarBarActor.h"
#include "vtkTextActor.h"
#include "vtkTextProperty.h"
#include "vtkLookupTable.h"
#include "vtkWindowToImageFilter.h"
#include "vtkPNGWriter.h"
#include "vtkScalarsToColors.h"
#include "vtkMath.h"
#include "vtkPoints.h"
#include "ModelLookupTables.h"

#include <iostream>
#include <cmath>
#include <string>
#include <sstream>

using namespace std;

int main(){
	//
	// model run parameters (edit this section for specific problems)
	//
	double dt = 4*0.01;
	int frameOut = 1;
	int finalFrame = 350;
	double relVortMin = -6.71E-05;
	double relVortMax = 6.71E-05;
	char scalarsDataName[16] = "relVortPanel";
	char scalarsDataTitle[28] = "Rel. Vort.";

	char vtkFileroot[256] = "/Users/pbosler/Desktop/modelData/vtkOut/rhWave_remesh20_noReinit_quadUnif5__";
	char imgFileroot[256] = "/Users/pbosler/Desktop/modelData/imgOut/rhWave_remesh20_noReinit_";
	
	int numColorsInColorbars = 1024;
	//
	// end edit section 
	//
	
	const double EARTH_RADIUS = 0.9 ; //6371.220; // kilometers
	

	//
	// initialize .vtk file readers
	//
	char counterString[4];
	int j;
	j = sprintf(counterString,"%04d",0);
	char vtkFilename[256];
	strcpy(vtkFilename, vtkFileroot);
	strcat(vtkFilename,counterString);
	strcat(vtkFilename,".vtk");
	
	cout << "reading " << finalFrame + 1 << " .vtk data files beginning with : " << endl;
	cout << "        " << vtkFilename << endl;
	
	vtkPolyDataReader * relVortData[finalFrame + 1];
	vtkPolyDataReader * meshData[finalFrame + 1];
	
	for ( int frameJ = 0; frameJ <= finalFrame; frameJ++)
	{
		j = sprintf(counterString,"%04d",frameJ);
		strcpy(vtkFilename, vtkFileroot);
		strcat(vtkFilename,counterString);
		strcat(vtkFilename,".vtk");
		
		relVortData[frameJ] = vtkPolyDataReader::New();
		meshData[frameJ] = vtkPolyDataReader::New();
		
		relVortData[frameJ]->SetFileName(vtkFilename);
		meshData[frameJ]->SetFileName(vtkFilename);
		
		relVortData[frameJ]->SetScalarsName(scalarsDataName);		
	}
	
	//
	// set up scalar colorbars
	//
	vtkLookupTable * relVortScale = vtkLookupTable::New();
		relVortScale->SetRange(relVortMin, relVortMax);
		relVortScale->SetNumberOfTableValues(numColorsInColorbars);
		double r, g, b, val;
		double range = relVortMax - relVortMin;
		for (int i = 0; i < numColorsInColorbars; i++)
		{
			val = relVortMin + ( (double)i/numColorsInColorbars) * range;
			GetColorForValue_BlueWhiteRed( val, r, g, b, relVortMin, relVortMax);
			relVortScale -> SetTableValue(i, r, g, b);
		}
		relVortScale -> Build();
	vtkTextProperty * colorbarLabelProperties = vtkTextProperty::New();
		colorbarLabelProperties -> SetColor(0.0, 0.0, 0.0);
		colorbarLabelProperties -> SetFontSize(12);	
	
	//
	// VTK graphics pipeline
	//
	vtkPolyDataMapper * relVortMapper = vtkPolyDataMapper::New();
		relVortMapper -> SetScalarRange(relVortMin, relVortMax);
		relVortMapper -> SetLookupTable(relVortScale);
	vtkPolyDataMapper * meshMapper = vtkPolyDataMapper::New();
		meshMapper -> ScalarVisibilityOff();
		
	vtkActor * relVortSphere = vtkActor::New();
		relVortSphere -> SetMapper(relVortMapper);
	vtkActor * meshSphere = vtkActor::New();		
		meshSphere -> SetMapper(meshMapper);
		meshSphere -> GetProperty() -> EdgeVisibilityOn();
		meshSphere -> GetProperty() -> SetEdgeColor(0.0, 0.0, 0.0);
		meshSphere -> GetProperty() -> SetColor(1.0, 1.0, 1.0);
	
	vtkScalarBarActor * relVortColorbar = vtkScalarBarActor::New();
		relVortColorbar -> SetLookupTable(relVortScale);
		relVortColorbar -> SetOrientationToVertical();
		relVortColorbar -> SetLabelTextProperty(colorbarLabelProperties);
		relVortColorbar -> SetPosition(0.015, 0.1);
		relVortColorbar -> SetPosition2(0.9, 0.9);
		relVortColorbar -> SetWidth(0.25);
	
	vtkTextActor * annotationText = vtkTextActor::New();
		annotationText -> SetTextScaleModeToNone();
		annotationText -> GetTextProperty() -> SetColor(0.0, 0.0, 0.0);
		annotationText -> GetTextProperty() -> SetFontSize(40);
		annotationText -> SetPosition(0.2, 0.1);
		annotationText -> SetPosition2(0.7, 0.25);
	
	vtkCamera * camera1 = vtkCamera::New();
		camera1 -> SetPosition(5.0 * EARTH_RADIUS, 0.0 * EARTH_RADIUS, 2.0 * EARTH_RADIUS);
		camera1 -> SetViewUp( 0.0, 0.0, 1.0);
		camera1 -> SetFocalPoint( 0.0, -0.25, 0.0);
		
	vtkCamera * camera2 = vtkCamera::New();
		camera2 -> SetPosition(5.0 * EARTH_RADIUS, 0.0 * EARTH_RADIUS, 2.0 * EARTH_RADIUS);
		camera2 -> SetViewUp( 0.0, 0.0, 1.0);
		camera2 -> SetFocalPoint( 0.0, 0.0, 0.0);	
	
	vtkRenderer * relVortRenderer = vtkRenderer::New();
		relVortRenderer -> AddActor(relVortSphere);
		relVortRenderer -> AddActor(relVortColorbar);
		relVortRenderer -> AddActor(annotationText);
		relVortRenderer -> SetBackground(1.0, 1.0, 1.0);
		relVortRenderer -> SetViewport(0.0, 0.0, 0.5, 1.0);
		relVortRenderer -> SetActiveCamera(camera1);
		
	vtkRenderer * meshRenderer = vtkRenderer::New();
		meshRenderer -> AddActor(meshSphere);
		meshRenderer -> SetBackground( 1.0, 1.0, 1.0);
		meshRenderer -> SetViewport( 0.5, 0.0, 1.0, 1.0);
		meshRenderer -> SetActiveCamera(camera2);
	
	vtkRenderWindow * renWin = vtkRenderWindow::New();
		renWin -> AddRenderer(relVortRenderer);
		renWin -> AddRenderer(meshRenderer);
		renWin -> SetSize(1200, 600);
	
	vtkWindowToImageFilter * win2im = vtkWindowToImageFilter::New();
		win2im -> SetInput(renWin);
	
	vtkPNGWriter * writer = vtkPNGWriter::New();
		char imgFilename[256];
	
	cout << "RENDERING ... " << endl;
	
	int nActive = 0;
	int imgCounter = 0;
	for ( int frameJ = 0; frameJ <= finalFrame; frameJ+=frameOut)
	{
		char annotationString[64];		
		j = sprintf(annotationString, "t = %5.2f, N = %6d", frameJ*dt, nActive);
		annotationText -> SetInput(annotationString);
		
		relVortMapper -> SetInput( relVortData[frameJ]->GetOutput());
		meshMapper -> SetInput( meshData[frameJ]->GetOutput());
		
		renWin -> Render();
		win2im -> Modified();	
	
		nActive = relVortData[frameJ]->GetOutput()->GetNumberOfCells();
	
		j = sprintf(counterString,"%04d",imgCounter);
		strcpy(imgFilename, imgFileroot);
		strcat(imgFilename, counterString);
		strcat(imgFilename,".png");
		
		writer -> SetInput(win2im->GetOutput());
		writer -> SetFileName(imgFilename);
		writer -> Write();	
		
		imgCounter++;
		
		if ( frameJ == 0 )
		{
			j = sprintf(annotationString, "t = %5.2f, N = %6d", frameJ*dt, nActive);
			annotationText -> SetInput(annotationString);
		
			relVortMapper -> SetInput( relVortData[frameJ]->GetOutput());
			meshMapper -> SetInput( meshData[frameJ]->GetOutput());
		
			renWin -> Render();
			win2im -> Modified();	
	
			j = sprintf(counterString,"%04d",frameJ);
			strcpy(imgFilename, imgFileroot);
			strcat(imgFilename, counterString);
			strcat(imgFilename,".png");
		
			writer -> SetInput(win2im->GetOutput());
			writer -> SetFileName(imgFilename);
			writer -> Write();	
		}
		if ( frameJ == finalFrame/4){cout << "... 25\% done.\n";}
		if ( frameJ == finalFrame/2){cout << "... 50\% done.\n";}
		if ( frameJ == 3*finalFrame/4){cout << "... 75\% done.\n";}
			
	}

	// 
	// clean up
	//		
	camera1 -> Delete();
	camera2 -> Delete();
	writer -> Delete();
	win2im -> Delete();
	renWin -> Delete();
	meshSphere -> Delete();
	relVortSphere -> Delete();
	annotationText -> Delete();
	relVortColorbar -> Delete();
	relVortRenderer -> Delete();
	meshRenderer -> Delete();
	colorbarLabelProperties -> Delete();
	meshMapper -> Delete();
	relVortMapper -> Delete();
	relVortScale -> Delete();
	for ( int frameJ=0; frameJ <= finalFrame; frameJ++)
	{
		relVortData[frameJ] -> Delete();
		meshData[frameJ] -> Delete();
	}
			
	
	
	
return 0;
};
