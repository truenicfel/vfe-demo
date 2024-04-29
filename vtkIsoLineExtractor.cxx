#include "vtkIsoLineExtractor.h"

#include <vtkImageData.h>
#include <vtkObjectFactory.h>
#include <vtkInformationVector.h>
#include <vtkInformation.h>
#include <vtkDataObject.h>
#include <vtkSmartPointer.h>
#include <vtkPointData.h>
#include <vtkDataArray.h>
#include <vtkVector.h>
#include <vtkDoubleArray.h>
#include <vtkPolyData.h>
#include <vtkPoints.h>
#include <vtkVertex.h>
#include <vtkPolyLine.h>

#include <iostream>

vtkStandardNewMacro(vtkIsoLineExtractor);

vtkIsoLineExtractor::vtkIsoLineExtractor()
: GrowingIterations(800)
, StepSizeGrowing(0.5)
, VariationalIterations(250)
, StepSizeVariational(0.01)
, ProximityThreshold(0.1)
, Iso(2800)
, SmoothnessWeight(500)
, BorderRefinement(25)
, FieldNameScalars("scalars")
, FieldNameGradient("gradient")
, stepSizeRefinement(0.0)
, stepSizeGrowing(0.0)
{
}

void vtkIsoLineExtractor::SetImageData(vtkSmartPointer<vtkImageData> imageData)
{

	ImageData = imageData;
	
}


int vtkIsoLineExtractor::RequestData(vtkInformation *vtkNotUsed(request),
	vtkInformationVector **inputVector,
	vtkInformationVector *outputVector)
{
	using namespace Eigen;

	// Get the info objects
	vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
	vtkInformation *outInfo = outputVector->GetInformationObject(0);

	// Get the input and output
	vtkPolyData* input = vtkPolyData::SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));
	vtkPolyData* output = vtkPolyData::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));

    // check inputs
	if (ImageData == nullptr)
	{
		std::cout << "Image data not set!" << std::endl;
		return 0;
	}
	if (ImageData->GetPointData()->GetArray(FieldNameScalars) == nullptr || ImageData->GetPointData()->GetArray(FieldNameScalars)->GetNumberOfComponents() != 1)
	{
		std::cout << "Field " << FieldNameScalars << " was not found or does not have 1 components!" << std::endl;
		return 0;
	}
	if (ImageData->GetPointData()->GetArray(FieldNameGradient) == nullptr || ImageData->GetPointData()->GetArray(FieldNameGradient)->GetNumberOfComponents() != 2)
	{
		std::cout << "Field " << FieldNameGradient << " was not found or does not have 2 components!" << std::endl;
		return 0;
	}

	// Get the information on resolution and domain
	int* dimensionsRaw = ImageData->GetDimensions();
	Eigen::Vector2i dimensions(dimensionsRaw[0], dimensionsRaw[1]);
	double* spacingRaw = ImageData->GetSpacing();
	Eigen::Vector2d spacing(spacingRaw[0], spacingRaw[1]);
	double* originRaw = ImageData->GetOrigin();
	Eigen::Vector2d origin(originRaw[0], originRaw[1]);

    // step sizes for growing and refinement based on grid spacing
    stepSizeRefinement = spacing.minCoeff() * StepSizeVariational;
    stepSizeGrowing = spacing.minCoeff() * StepSizeGrowing;

    // extract seed points and extract a line for each of them
    vtkSmartPointer<vtkPoints> seedPoints = input->GetPoints();
    auto numberOfPoints = seedPoints->GetNumberOfPoints();
    std::cout << numberOfPoints << " seed point(s)." << std::endl;
    std::vector<std::deque<Vector2d>> extractedLines;
    for (vtkIdType pointIndex = 0; pointIndex < numberOfPoints; ++pointIndex) {
        double point[3] = {0.0, 0.0, 0.0};
        seedPoints->GetPoint(pointIndex, point);
        Eigen::Vector2d seedPoint(point[0], point[1]);
        std::cout << "Growing from point: " << std::endl << seedPoint << std::endl;

        // extraction
        auto line = extract(seedPoint);

        extractedLines.push_back(line);

        std::cout << "Extracted critical line with length " << line.size() << "." << std::endl;

    }

    // get the total number of points
    vtkIdType totalNumPoints = 0;
    size_t numLines = extractedLines.size();
    for (const auto& line: extractedLines) {
        totalNumPoints += static_cast<vtkIdType>(line.size());
    }

    // convert to vtk data structure
    vtkNew<vtkPolyData> polyData;
    vtkNew<vtkCellArray> cellArray;
    vtkNew<vtkPoints> points;
    points->SetNumberOfPoints(totalNumPoints);
    vtkIdType cpNumber = 0;
    for (size_t iLine = 0; iLine < numLines; ++iLine)
    {
        auto line = extractedLines[iLine];

        auto numPoints = static_cast<vtkIdType>(line.size());
        if (numPoints <= 1)
            continue;

        vtkNew<vtkPolyLine> polyLine;
        polyLine->GetPointIds()->SetNumberOfIds(numPoints);
        for (vtkIdType iPoint = 0; iPoint < numPoints; ++iPoint)
        {
            Vector3d tempPoint(line[iPoint].x(), line[iPoint].y(), 0.0);
            points->SetPoint(cpNumber, tempPoint.data());

            polyLine->GetPointIds()->SetId(iPoint, cpNumber);
            ++cpNumber;
        }
        cellArray->InsertNextCell(polyLine);
    }
    points->ComputeBounds();
    polyData->SetPoints(points);
    polyData->SetLines(cellArray);
	
	// Copy the computed image to the output
	output->ShallowCopy(polyData);

	output->ComputeBounds();

	return 1;
}

std::deque<Eigen::Vector2d> vtkIsoLineExtractor::extract(const Eigen::Vector2d& seed) const {

    // create a new line that only contains the seed
    DoubleBufferedLine2 line(seed);

    int timesGrown = 0;

    // define the sample gradient function
    Gradient<2>::SampleGradientFunction sampleGradientFunction = [&](const DoubleBufferedLine2& line, const int& index) {

        const Eigen::Vector2d& position = line.Read(index);
        // border point?
        bool border = index == 0 || index == (line.GetSize() - 1);

        auto scalar = SampleScalar(position);
        auto gradient = SampleGradient(position);

        // EL gradient for iso
        Eigen::Vector2d result = gradient * (scalar - Iso);

        if (!border)
        {
            // left and right neighbours for second derivative estimation
            Eigen::Vector2d left = line.Read(index - 1);
            Eigen::Vector2d right = line.Read(index + 1);
            // EL gradient for smoothness
            result += SmoothnessWeight * -(left - 2 * position + right);
        }
        return result;
    };

    while (timesGrown < GrowingIterations) {

        // grow
        grow(line, true);
        grow(line, false);
        timesGrown += 1;

        // gradient (one forward one backward)
        Gradient<2> gradientForward(sampleGradientFunction, BorderRefinement, true, false, stepSizeRefinement);
        Gradient<2> gradientBackward(sampleGradientFunction, BorderRefinement, false, false, stepSizeRefinement);

        // optimize once but VariationalIterations times. in between we check if we left the domain
        int timesRefined = 0;
        while (timesRefined < VariationalIterations) {

            // optimize forward and backward
            gradientForward.OptimizeOnce(line);
            gradientBackward.OptimizeOnce(line);

            // setup for next iteration
            timesRefined += 1;
            line.Swap(); // in the end the up-to-date line will be the one we read from
        }

    }

    return line.GetCopy();
}

void vtkIsoLineExtractor::grow(DoubleBufferedLine2& line, const bool& forward) const {

    Eigen::Vector2d direction(0.0, 1.0);
    if (!forward)
    {
        direction = Eigen::Vector2d(0.0, -1.0);
    }

    size_t size = line.GetSize();

    if (size > 1) {
        if (forward) {
            direction = (line.Read(size - 1) - line.Read(size - 2)).stableNormalized();
        } else {
            direction = (line.Read(0) - line.Read(1)).stableNormalized();
        }

    }

    Eigen::Vector2d newPoint = line.Read(size - 1);

    if (!forward)
    {
        newPoint = line.Read(0);
    }

    newPoint += direction * stepSizeGrowing;

    if (forward) {
        line.Append(newPoint);
    } else {
        line.Prepend(newPoint);
    }

}

// sample the image data directly but only scalar
double vtkIsoLineExtractor::SampleScalar(const Eigen::Vector2d& coordinate) const
{

    vtkSmartPointer<vtkDataArray> array = ImageData->GetPointData()->GetArray(FieldNameScalars);

    // Get the information on the domain
    int* dimensionsRaw = ImageData->GetDimensions();
    Eigen::Vector2i dimensions(dimensionsRaw[0], dimensionsRaw[1]);
    double* spacingRaw = ImageData->GetSpacing();
    Eigen::Vector2d spacing(spacingRaw[0], spacingRaw[1]);
    double* originRaw = ImageData->GetOrigin();
    Eigen::Vector2d domainMin(originRaw[0], originRaw[1]);
    Eigen::Vector2d domainMax = domainMin + (dimensions.template cast<double>() - Eigen::Vector2d::Ones()).cwiseProduct(spacing);

    Eigen::Vector2d vf_tex = (coordinate - domainMin).cwiseQuotient(domainMax - domainMin);
    Eigen::Vector2d vf_sample = vf_tex.cwiseProduct(dimensions.template cast<double>() - Eigen::Vector2d::Ones());

    Eigen::Vector2i vi_sample_base0, vi_sample_base1;
    for (int i = 0; i < 2; ++i) {
        vi_sample_base0[i] = std::min(std::max(0, (int)vf_sample[i]), dimensions[i] - 1);
        vi_sample_base1[i] = std::min(std::max(0, vi_sample_base0[i] + 1), dimensions[i] - 1);
    }
    Eigen::Vector2d vf_sample_interpol = vf_sample - vi_sample_base0.template cast<double>();

    int num_corners = (int)std::pow(2, 2);
    double result = 0.0;
    for (int i = 0; i < num_corners; ++i) {
        double weight = 1;
        Eigen::Vector2i grid_index = Eigen::Vector2i::Zero();
        for (int d = 0; d < 2; ++d) {
            if (i & (int(1) << (2 - 1 - d))) {
                grid_index[d] = vi_sample_base1[d];
                weight *= vf_sample_interpol[d];
            }
            else {
                grid_index[d] = vi_sample_base0[d];
                weight *= 1 - vf_sample_interpol[d];
            }
        }
        //result += this->value(grid_index) * weight;
        //std::cout << grid_index.x() + grid_index.y() * dimensions.x() + grid_index.z() * dimensions.x() * dimensions.y() << std::endl;
        double value = array->GetTuple1(grid_index.x()
                + grid_index.y() * dimensions.x());
        result += value * weight;
    }
    //std::cout << "##################################################" << std::endl;
    return result;
}

// sample the image data directly but only the gradient
Eigen::Vector2d vtkIsoLineExtractor::SampleGradient(const Eigen::Vector2d& coordinate) const
{

    vtkSmartPointer<vtkDataArray> array = ImageData->GetPointData()->GetArray(FieldNameGradient);

    // Get the information on the domain
    int* dimensionsRaw = ImageData->GetDimensions();
    Eigen::Vector2i dimensions(dimensionsRaw[0], dimensionsRaw[1]);
    double* spacingRaw = ImageData->GetSpacing();
    Eigen::Vector2d spacing(spacingRaw[0], spacingRaw[1]);
    double* originRaw = ImageData->GetOrigin();
    Eigen::Vector2d domainMin(originRaw[0], originRaw[1]);
    Eigen::Vector2d domainMax = domainMin + (dimensions.template cast<double>() - Eigen::Vector2d::Ones()).cwiseProduct(spacing);

    Eigen::Vector2d vf_tex = (coordinate - domainMin).cwiseQuotient(domainMax - domainMin);
    Eigen::Vector2d vf_sample = vf_tex.cwiseProduct(dimensions.template cast<double>() - Eigen::Vector2d::Ones());

    Eigen::Vector2i vi_sample_base0, vi_sample_base1;
    for (int i = 0; i < 2; ++i) {
        vi_sample_base0[i] = std::min(std::max(0, (int)vf_sample[i]), dimensions[i] - 1);
        vi_sample_base1[i] = std::min(std::max(0, vi_sample_base0[i] + 1), dimensions[i] - 1);
    }
    Eigen::Vector2d vf_sample_interpol = vf_sample - vi_sample_base0.template cast<double>();

    int num_corners = (int)std::pow(2, 2);
    Eigen::Vector2d result = Eigen::Vector2d::Zero();
    for (int i = 0; i < num_corners; ++i) {
        double weight = 1;
        Eigen::Vector2i grid_index = Eigen::Vector2i::Zero();
        for (int d = 0; d < 2; ++d) {
            if (i & (int(1) << (2 - 1 - d))) {
                grid_index[d] = vi_sample_base1[d];
                weight *= vf_sample_interpol[d];
            }
            else {
                grid_index[d] = vi_sample_base0[d];
                weight *= 1 - vf_sample_interpol[d];
            }
        }
        //result += this->value(grid_index) * weight;
        //std::cout << grid_index.x() + grid_index.y() * dimensions.x() + grid_index.z() * dimensions.x() * dimensions.y() << std::endl;
        double* value = array->GetTuple2(grid_index.x()
                                        + grid_index.y() * dimensions.x());
        result += Eigen::Vector2d(value[0], value[1]) * weight;
    }
    //std::cout << "##################################################" << std::endl;
    return result;
}