#include <iostream>

#include <vtkImageData.h>
#include <vtkPoints.h>
#include <vtkSmartPointer.h>
#include <vtkXMLImageDataReader.h>
#include <vtkVertex.h>

#include "vtkIsoLineExtractor.h"
#include <vtkXMLPolyDataWriter.h>



// Input arguments:
//   [1] : path to the vti-file in XML format, containing the unsteady 2d vector field (vectors, x partials, y partials and time partials)
int main(int argc, char* argv[])
{
	if (argc < 2)
	{
		std::cout << "Please supply the required input argument!" << std::endl;
		return 0;
	}
	if (argc > 2)
	{
		std::cout << "Note: Only one argument is required. Further arguments are ignored." << std::endl;
	}

	std::string imageDataFile = argv[1];

	// ---------------------------------------------
	// -------- Initialization of test case --------
	// ---------------------------------------------
	std::cout << "Reading input data..." << std::endl;


	vtkSmartPointer<vtkXMLImageDataReader> imageDataReader = vtkSmartPointer<vtkXMLImageDataReader>::New();
	imageDataReader->SetFileName(imageDataFile.c_str());
	imageDataReader->UpdateInformation();
	imageDataReader->Update();
	vtkSmartPointer<vtkImageData> imageData = imageDataReader->GetOutput();

	// setup for seed point
	vtkNew<vtkPoints> points;
	points->SetNumberOfPoints(1);
	vtkNew<vtkPolyData> polyData;
	vtkNew<vtkCellArray> vertices;
	points->SetNumberOfPoints(1);

	// add the seed point
	auto point = Eigen::Vector2d(0.56, 0.0);
	points->SetPoint(0, point.x(), point.y(), 0.0);
	vtkNew<vtkVertex> vertex;
	vertex->GetPointIds()->SetId(0, 0);
	vertices->InsertNextCell(vertex);
	points->ComputeBounds();

	// add to polyData
	polyData->SetPoints(points);
	polyData->SetVerts(vertices);


	// ---------------------------------------------
	// --- Compute iso line ---
	// ---------------------------------------------
	std::cout << "Computing iso lines..." << std::endl;

	vtkNew<vtkIsoLineExtractor> filter;
	filter->SetInputData(polyData);
	filter->SetImageData(imageData);
	filter->SetGrowingIterations(1500);
	filter->SetStepSizeGrowing(0.5);
	filter->SetVariationalIterations(250);
	filter->SetStepSizeVariational(0.01);
	filter->SetProximityThreshold(0.01); // Unused
	filter->SetIso(2800);
	filter->SetSmoothnessWeight(500);
	filter->SetBorderRefinement(25);
	filter->SetFieldNameScalars("scalars");
	filter->SetFieldNameGradient("gradient");
	filter->Update();
	vtkSmartPointer<vtkPolyData> output = filter->GetOutput();


	// write to file
	vtkNew<vtkXMLPolyDataWriter> writer;
	writer->SetFileName("isoLine.vtp");
	writer->SetInputData(output);
	writer->Update();

	cout << "Done." << endl;
	return 0;
}