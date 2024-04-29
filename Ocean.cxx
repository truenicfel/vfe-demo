#include <iostream>

#include <vtkImageData.h>
#include <vtkSmartPointer.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkXMLImageDataReader.h>

#include "vtkCriticalLinesExtractor.h"
#include <vtkXMLPolyDataWriter.h>



// Input arguments:
//   [1] : path to the vti-file in XML format, containing the unsteady 2d vector field (vectors, x partials, y partials and time partials)
//   [2] : path to the vtp-file in XML format, containing the seed points
int main(int argc, char* argv[])
{
	if (argc < 3)
	{
		std::cout << "Please supply both required input arguments!" << std::endl;
		return 0;
	}
	if (argc > 3)
	{
		std::cout << "Note: Only two arguments are required. Further arguments are ignored." << std::endl;
	}

	std::string imageDataFile = argv[1];
	std::string pointDataFile = argv[2];

	// ---------------------------------------------
	// -------- Initialization of test case --------
	// ---------------------------------------------
	std::cout << "Reading input data..." << std::endl;


	vtkSmartPointer<vtkXMLImageDataReader> imageDataReader = vtkSmartPointer<vtkXMLImageDataReader>::New();
	imageDataReader->SetFileName(imageDataFile.c_str());
	imageDataReader->UpdateInformation();
	imageDataReader->Update();
	vtkSmartPointer<vtkImageData> imageData = imageDataReader->GetOutput();


	// read poly data
	vtkNew<vtkXMLPolyDataReader> polyDataReader;
	polyDataReader->SetFileName(pointDataFile.c_str());
	polyDataReader->Update();
	vtkSmartPointer<vtkPolyData> points = polyDataReader->GetOutput();

	// ---------------------------------------------
	// --- Compute vector field in optimal frame ---
	// ---------------------------------------------
	std::cout << "Computing critical lines..." << std::endl;

	vtkNew<vtkCriticalLinesExtractor> filter;
	filter->SetInputData(points);
	filter->SetImageData(imageData);
	filter->SetGrowingIterations(200);
	filter->SetTimeStep(0.05);
	filter->SetVariationalIterations(300);
	filter->SetAdamLearningRate(0.001);
	filter->SetFieldNameV("vectors");
	filter->SetFieldNameVx("xPartials");
	filter->SetFieldNameVy("yPartials");				
	filter->SetFieldNameVt("tPartials");				
	filter->SetUseComputeBoundingBox(true);
	// "safety net": avoid sampling where the grid is not defined
	filter->SetBBMin(vtkVector3d(-3.9, -33.0, 0.00));
	filter->SetBBMax(vtkVector3d(0.5, -28.0, 36.0));
	filter->Update();
	vtkSmartPointer<vtkPolyData> output = filter->GetOutput();


	// write to file
	vtkNew<vtkXMLPolyDataWriter> writer;
	writer->SetFileName("criticalLines.vtp");
	writer->SetInputData(output);
	writer->Update();

	cout << "Done." << endl;
	return 0;
}