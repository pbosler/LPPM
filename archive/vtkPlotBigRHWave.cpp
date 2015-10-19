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
#include <string.h>
#include <sstream>

using namespace std;

int main(){
	//
	// EARTH CONSTANTS 
	//
	double EARTH_RADIUS = 6371220.0;
	double ONE_DAY = 86140.0;
	double GRAVITY = 9.80616;
	double OMEGA = 2.0*3.141592653589793/ONE_DAY;
	double bounds[6];
	
	//
	// USER PARAMETERS 
	// TO DO : Change this section to a user input file (e.g., namelist)
	//
	double dt = 0.01;
	int frameOut = 1;
	int finalFrame = 35;
	double relVortMin = -6.71E-05;
	double relVortMax = 6.71E05;
	bool showEdges = true;
	char vtkFileRoot[128] = "/Users/pbosler/Desktop/modelData/vtkOut/rhWave_noRemesh_quadUnif5__";
	char imgFileRoot[128] = "/Users/pbosler/Desktop/modelData/imgOut/rhWave_noRemesh_2day_";
	/* END USER PARAMETERS */
	
	//
	// initialize input / output
	//
	char outputCounterString[4];
	int outputCounter;

	vtkPolyDataReader* relVortData[finalFrame + 1];
	vtkPolyDataReader* meshData[finalFrame + 1];
	vtkPolyDataReader* flowMapData[finalFrame + 1];
	
	cout << "initializing vtk file readers ... ";
	int j;
	char vtkFileName[128];
	char imgFileName[128];
	for (outputCounter = 0; outputCounter <= finalFrame; outputCounter++ ){
		// set file names
		j = sprintf(outputCounterString,"%04d",outputCounter);
		strcpy(vtkFileName,vtkFileRoot);
		strcat(vtkFileName,outputCounterString);
		strcat(vtkFileName,".vtk");
		
		// initialize vtkPolyData readers
		relVortData[outputCounter] = vtkPolyDataReader::New();
			relVortData[outputCounter] -> SetFileName(vtkFileName);
			relVortData[outputCounter] -> SetScalarsName("RelVort");
		
		meshData[outputCounter] = vtkPolyDataReader::New();
			meshData[outputCounter] -> SetFileName(vtkFileName);
			meshData[outputCounter] -> SetScalarsName("Tracer2");
			
		flowMapData[outputCounter] = vtkPolyDataReader::New();
			flowMapData[outputCounter] -> SetFileName(vtkFileName);
			flowMapData[outputCounter] -> SetScalarsName("Tracer1");
	}
	cout << "... done. \n";
	
	cout << "initializing lookup tables ... ";
	int numColors = 1024;
	double r, g, b, val;
	vtkLookupTable* vortScale = vtkLookupTable::New();
		vortScale->SetRange(relVortMin,relVortMax);
		double vortRange = relVortMax - relVortMin;
		vortScale->SetNumberOfTableValues(numColors);
		for (int i = 0; i < numColors; i++){
			val = relVortMin + ( (double)i/numColors)*vortRange;
			//
			// choose colormap for vorticity
			//
			GetColorForValue_BlueWhiteRed(val,r,g,b,relVortMin,relVortMax);
			
			vortScale -> SetTableValue(i,r,g,b);
		}
	vortScale -> Build();
	
	double latMin = -1.5708;
	double latMax = 1.5708;
	double latRange = latMax - latMin;
	vtkLookupTable* latScale = vtkLookupTable::New();
		latScale -> SetRange(latMin,latMax);
		latScale -> SetNumberOfTableValues(numColors);
		for ( int i = 0; i < numColors; i++){
			val = latMin + ( (double)i/numColors)*latRange;
			//
			// choose colormap for latitude tracers
			//
			GetColorForValue_NCLDetail(val,r,g,b,latMin,latMax);
			
			latScale -> SetTableValue(i,r,g,b);
		}	
	latScale -> Build();
	cout << "... done. \n";	
	
	//
	//	initialize graphics pipeline
	//
	cout << "initializing graphics pipeline ... ";
		//
		//	data mappers
		//
		vtkPolyDataMapper* sphereMapper1 = vtkPolyDataMapper::New();
			sphereMapper1 -> SetScalarRange(relVortMin,relVortMax);
			sphereMapper1 -> SetLookupTable(vortScale);
		vtkPolyDataMapper* sphereMapper2 = vtkPolyDataMapper::New();
			sphereMapper2 -> SetScalarRange(latMin,latMax);
			sphereMapper2 -> SetLookupTable(latScale);
		vtkPolyDataMapper* sphereMapper3 = vtkPolyDataMapper::New();
			sphereMapper3 -> SetScalarRange(latMin,latMax);
			sphereMapper3 -> SetLookupTable(latScale);
		
		//
		// graphics actors
		//
		vtkActor* sphereActor1 = vtkActor::New();
			sphereActor1 -> SetMapper(sphereMapper1);
		vtkActor* sphereActor2 = vtkActor::New();
			sphereActor2 -> SetMapper(sphereMapper2);
		vtkActor* sphereActor3 = vtkActor::New();		
			sphereActor3 -> SetMapper(sphereMapper3);
			if ( showEdges ) {
				sphereActor3 -> GetProperty() -> EdgeVisibilityOn();
				sphereActor3 -> GetProperty() -> SetEdgeColor(0.0,0.0,0.0);
				sphereActor3 -> GetProperty() -> SetColor(1.0,1.0,1.0);
			}
		
		//
		// colorbar actors
		//
		vtkTextProperty* colorbarLabelProperties = vtkTextProperty::New();
			colorbarLabelProperties->SetColor(0.0,0.0,0.0);
			colorbarLabelProperties->SetFontSize(12);
		vtkScalarBarActor* vortColorbar = vtkScalarBarActor::New();
			vortColorbar->SetTitle("Rel. Vort. < s^(-1) >");
			vortColorbar->GetTitleTextProperty()->SetColor(0.0,0.0,0.0);
			vortColorbar->GetTitleTextProperty()->SetFontSize(12);
			vortColorbar->SetLabelTextProperty(colorbarLabelProperties);
			vortColorbar->SetOrientationToVertical();
			vortColorbar->SetPosition(0.015,0.1);
			vortColorbar->SetPosition2(0.9,0.9);
			vortColorbar->SetWidth(0.2);
			vortColorbar->SetLookupTable(vortScale);
		vtkScalarBarActor* latColorbar = vtkScalarBarActor::New();
			latColorbar->SetTitle("Init. Latitude");
			latColorbar->GetTitleTextProperty()->SetColor(0.0,0.0,0.0);
			latColorbar->GetTitleTextProperty()->SetFontSize(12);
			latColorbar->SetLabelTextProperty(colorbarLabelProperties);
			latColorbar->SetOrientationToVertical();
			latColorbar->SetPosition(0.015,0.1);
			latColorbar->SetPosition2(0.9,0.9);
			latColorbar->SetWidth(0.2);
			latColorbar->SetLookupTable(latScale);
		vtkScalarBarActor* alphaColorbar = vtkScalarBarActor::New();
			alphaColorbar->SetTitle("Lat(alpha)");
			alphaColorbar->GetTitleTextProperty()->SetColor(0.0,0.0,0.0);
			alphaColorbar->GetTitleTextProperty()->SetFontSize(12);
			alphaColorbar->SetLabelTextProperty(colorbarLabelProperties);
			alphaColorbar->SetOrientationToVertical();
			alphaColorbar->SetPosition(0.015,0.1);
			alphaColorbar->SetPosition2(0.9,0.9);
			alphaColorbar->SetWidth(0.2);
			alphaColorbar->SetLookupTable(latScale);	
		
		//
		// text annotation
		//
		vtkTextActor* textActor1 = vtkTextActor::New();
			textActor1->SetTextScaleModeToProp();
			textActor1->SetPosition(0.1,0.1);
			textActor1->SetPosition2(0.8,0.3);
			textActor1->GetTextProperty()->SetColor(0.0,0.0,0.0);
			
		
		//
		// cameras
		//
		vtkCamera* camera1 = vtkCamera::New();
			camera1 -> SetPosition(8000,2.0,2.0);
			camera1 -> SetFocalPoint(0.0,0.0,0.0);
		vtkCamera* camera2 = vtkCamera::New();
			camera2 -> SetPosition(2.0*EARTH_RADIUS,2.0*EARTH_RADIUS,2.0*EARTH_RADIUS);
			//camera2 -> SetViewUp(0.0,0.0,1.0);
			camera2 -> SetFocalPoint(0.0,0.0,0.0);
			
		//
		// renderers
		//
		vtkRenderer* renderer1 = vtkRenderer::New();
			renderer1 -> AddActor(sphereActor1);
			renderer1 -> AddActor(vortColorbar);
			renderer1 -> AddActor(textActor1);
			renderer1 -> SetActiveCamera(camera1);
			renderer1 -> SetBackground(1.0,1.0,1.0);
			renderer1 -> SetViewport(0.0,0.0,0.3333,1.0);
		vtkRenderer* renderer2 = vtkRenderer::New();
			renderer2 -> AddActor(sphereActor2);
			renderer2 -> AddActor(latColorbar);
			renderer2 -> SetActiveCamera(camera2);
			renderer2 -> SetBackground(1.0,1.0,1.0);
			renderer2 -> SetViewport(0.3333,0.0,0.6666,1.0);		
		vtkRenderer* renderer3 = vtkRenderer::New();
			renderer3 -> AddActor(sphereActor3);
			renderer3 -> AddActor(alphaColorbar);
			renderer3 -> SetActiveCamera(camera2);
			renderer3 -> SetBackground(1.0,1.0,1.0);
			renderer3 -> SetViewport(0.6666,0.0,1.0,1.0);		
		
		//
		// render window
		//
		vtkRenderWindow* renWin = vtkRenderWindow::New();
			renWin->AddRenderer(renderer1);
			renWin->AddRenderer(renderer2);
			renWin->AddRenderer(renderer3);
			renWin->SetSize(1800,600);
		
		//
		// window to image filter and image writer
		//
		vtkWindowToImageFilter* win2im = vtkWindowToImageFilter::New();
			win2im->SetInput(renWin);
		vtkPNGWriter* writer = vtkPNGWriter::New();
	cout << "... done. \n";	
	
	//
	// RENDERING LOOP
	//
	cout << "... RENDERING ... \n";
	vtkPoints* points1 = vtkPoints::New();
	vtkPoints* points2 = vtkPoints::New();
	vtkPoints* points3 = vtkPoints::New();
	vtkPolyData* polydata1 = vtkPolyData::New();
	vtkPolyData* polydata2 = vtkPolyData::New();
	vtkPolyData* polydata3 = vtkPolyData::New();
	int nActive = 0;
	char titleString[64];
	for ( int frameJ = 0; frameJ <= finalFrame; frameJ++){
		//
		// annotation string
		//
		j = sprintf(titleString,"t = %7.4f, N = %8d",frameJ*dt*frameOut,nActive);
		textActor1 -> SetInput(titleString);
		
		//
		// read data
		//
		polydata1->DeepCopy(relVortData[frameJ]->GetOutput());
		polydata2->DeepCopy(meshData[frameJ]->GetOutput());
		polydata3->DeepCopy(flowMapData[frameJ]->GetOutput());
	
		points1 = polydata1->GetPoints();
		points2 = polydata2->GetPoints();
		points3 = polydata3->GetPoints();
		
		//
		// normalize sphere to fit in single precision OpenGL used by VTK rendering
		//
		for ( vtkIdType i = 0; i < polydata1->GetNumberOfPoints(); i++){
			
			double xyz[3];
			double norm;
			
			points1->GetPoint(i,xyz);
			norm = sqrt( xyz[0]*xyz[0] + xyz[1]*xyz[1] + xyz[2]*xyz[2]);
			
			xyz[0]/=norm;
			xyz[1]/=norm;
			xyz[2]/=norm;
			
			points1->SetPoint(i,xyz);
			points2->SetPoint(i,xyz);
			points3->SetPoint(i,xyz);
		
		}
		polydata1->SetPoints(points1);
		polydata2->SetPoints(points2);
		polydata3->SetPoints(points3);
		sphereMapper1 -> SetInput(polydata1);
		sphereMapper2 -> SetInput(polydata2);
		sphereMapper3 -> SetInput(polydata3);
		
		//
		// render
		//
		renWin -> Render();
		win2im -> Modified();
		
		//
		// output frame to file
		//
		j = sprintf(outputCounterString,"%04d",frameJ);
		strcpy(imgFileName,imgFileRoot);
		strcat(imgFileName,outputCounterString);
		strcat(imgFileName,".png");
		
		writer -> SetInput(win2im->GetOutput());
		writer -> SetFileName(imgFileName);
		writer -> Write();
		
		// TO DO : nActive for AMR
		nActive = relVortData[frameJ]->GetOutput()->GetNumberOfCells();
		nActive = nActive/4;
		
		// 
		// track progress
		//
		if ( frameJ == finalFrame/4){cout << "... 25\% done.\n";}
		if ( frameJ == finalFrame/2){cout << "... 50\% done.\n";}
		if ( frameJ == 3*finalFrame/4){cout << "... 75\% done.\n";}
		
	}		
	cout << "... RENDERING COMPLETE.\n";
				
	//
	// clean up / free memory
	//
	polydata1->Delete();
	polydata2->Delete();
	polydata3->Delete();
	points1->Delete();
	points2->Delete();
	points3->Delete();
	win2im->Delete();
	writer->Delete();
	renWin->Delete();
	renderer1->Delete();
	renderer2->Delete();
	renderer3->Delete();
	camera1->Delete();
	camera2->Delete();
	textActor1->Delete();
	alphaColorbar->Delete();
	vortColorbar->Delete();
	latColorbar->Delete();
	colorbarLabelProperties->Delete();
	sphereActor1->Delete();
	sphereActor2->Delete();
	sphereActor3->Delete();
	sphereMapper1->Delete();
	sphereMapper2->Delete();
	sphereMapper3->Delete();
	vortScale->Delete();
	latScale->Delete();
	for (outputCounter = 0; outputCounter <= finalFrame; outputCounter++ ){
		// delete vtkPolyData readers
		relVortData[outputCounter] -> Delete();
		meshData[outputCounter] -> Delete();
		flowMapData[outputCounter] -> Delete();
	}
};