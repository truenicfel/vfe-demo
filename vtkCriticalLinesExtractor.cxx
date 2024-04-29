#include "vtkCriticalLinesExtractor.h"
#include <vtkImageData.h>
#include <vtkObjectFactory.h>
#include <vtkStreamingDemandDrivenPipeline.h>
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

vtkStandardNewMacro(vtkCriticalLinesExtractor);

vtkCriticalLinesExtractor::vtkCriticalLinesExtractor()
: GrowingIterations(200)
, TimeStep(0.05)
, VariationalIterations(300)
, AdamLearningRate(0.001)
, UseComputeBoundingBox(false)
, BBMin()
, BBMax()
, FieldNameV("vectors")
, FieldNameVx("xPartials")
, FieldNameVy("yPartials")
, FieldNameVt("tPartials")
{
}

void vtkCriticalLinesExtractor::SetImageData(vtkSmartPointer<vtkImageData> imageData)
{

	ImageData = imageData;
	
}


int vtkCriticalLinesExtractor::RequestData(vtkInformation *vtkNotUsed(request),
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

    // check fields...
	if (ImageData == nullptr)
	{
		std::cout << "Image data not set!" << std::endl;
		return 0;
	}
	if (ImageData->GetPointData()->GetArray(FieldNameV) == nullptr || ImageData->GetPointData()->GetArray(FieldNameV)->GetNumberOfComponents() != 3)
	{
		std::cout << "Field " << FieldNameV << " was not found or does not have 2 components!" << std::endl;
		return 0;
	}
	if (ImageData->GetPointData()->GetArray(FieldNameVx) == nullptr || ImageData->GetPointData()->GetArray(FieldNameVx)->GetNumberOfComponents() != 3)
	{
		std::cout << "Field " << FieldNameVx << " was not found or does not have 2 components!" << std::endl;
		return 0;
	}
	if (ImageData->GetPointData()->GetArray(FieldNameVy) == nullptr || ImageData->GetPointData()->GetArray(FieldNameVy)->GetNumberOfComponents() != 3)
	{
		std::cout << "Field " << FieldNameVy << " was not found or does not have 2 components!" << std::endl;
		return 0;
	}
	if (ImageData->GetPointData()->GetArray(FieldNameVt) == nullptr || ImageData->GetPointData()->GetArray(FieldNameVt)->GetNumberOfComponents() != 3)
	{
		std::cout << "Field " << FieldNameVt << " was not found or does not have 2 components!" << std::endl;
		return 0;
	}

	// Get the information on the domain
	int* dimensionsRaw = ImageData->GetDimensions();
	Eigen::Vector3i dimensions(dimensionsRaw[0], dimensionsRaw[1], dimensionsRaw[2]);
	double* spacingRaw = ImageData->GetSpacing();
	Eigen::Vector3d spacing(spacingRaw[0], spacingRaw[1], spacingRaw[2]);
	double* originRaw = ImageData->GetOrigin();
	Eigen::Vector3d origin(originRaw[0], originRaw[1], originRaw[2]);

    // adapt step size based on spacing
    stepSize = spacing.z() * TimeStep;


    // for each seed point extract a critical line
    vtkSmartPointer<vtkPoints> seedPoints = input->GetPoints();
    auto numberOfPoints = seedPoints->GetNumberOfPoints();
    std::cout << numberOfPoints << " seed point(s)." << std::endl;
    std::vector<std::deque<Vector3d>> extractedLines;
    for (vtkIdType pointIndex = 0; pointIndex < numberOfPoints; ++pointIndex) {
        double point[3] = {0.0, 0.0, 0.0};
        seedPoints->GetPoint(pointIndex, point);
        Eigen::Vector3d seedPoint(point[0], point[1],  point[2] + 0.2);
        std::cout << "Growing from point: " << std::endl << seedPoint << std::endl;

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

            points->SetPoint(cpNumber, line[iPoint].data());

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

std::deque<Eigen::Vector3d> vtkCriticalLinesExtractor::extract(const Eigen::Vector3d& seed) const {

    // create a new line
    DoubleBufferedLine3 line(seed);
    // define the bounding box that the computation will happen in
    Eigen::AlignedBox3d boundingBox;
    boundingBox.min() = Eigen::Vector3d(BBMin.GetX(), BBMin.GetY(), BBMin.GetZ());
    boundingBox.max() = Eigen::Vector3d(BBMax.GetX(), BBMax.GetY(), BBMax.GetZ());

    int timesGrown = 0;
    Eigen::Vector3d lastPoint = line.Read(line.GetSize() - 1);

    // define the sample gradient function
    Gradient<3>::SampleGradientFunction sampleGradientFunction = [&](const DoubleBufferedLine3& line, const int& index) {

        // sample line
        const Eigen::Vector3d& position = line.Read(index);

        // sample derivatives and vector
        auto v = Sample(position, "vectors");
        auto vx = Sample(position, "xPartials");
        auto vy = Sample(position, "yPartials");
        auto vt = Sample(position, "tPartials");
        // setup jacobian
        Eigen::Matrix<double, 3, 3> J;
        J.col(0) = vx;
        J.col(1) = vy;
        J.col(2) = vt;

        // Critical EL term
        Eigen::Vector3d result = J.transpose() * v;
        return result;
    };

    // keep going until domain is left or no more growing iterations
    while (timesGrown < GrowingIterations && !leftDomain<3>(lastPoint, boundingBox)) {

        // grow
        grow(line, true);
        timesGrown += 1;
        // remember the last point for bounding box check
        lastPoint = line.Read(line.GetSize() - 1);

        // gradient helper
        Gradient<3> gradient(sampleGradientFunction, 1, true, true, AdamLearningRate);

        // optimize once but VariationalIterations times. in between we check if we left the domain
        int timesRefined = 0;
        // optimize as much as required and keep checking if the domain was left
        while (timesRefined < VariationalIterations && !leftDomain<3>(lastPoint, boundingBox)) {

            gradient.OptimizeOnce(line);

            // setup for next iteration
            timesRefined += 1;
            lastPoint = line.Read(line.GetSize() - 1);
            line.Swap(); // in the end the up-to-date line will be the one we read from
        }

    }

    if (leftDomain<3>(lastPoint, boundingBox)) {
        std::cout << "Stopped growing because the domain was left!" << std::endl;
    } else {
        std::cout << "Stopped growing because the maximum number of growing steps was reached!" << std::endl;
    }

    return line.GetCopy();
}

void vtkCriticalLinesExtractor::grow(DoubleBufferedLine3& line, const bool& forward) const {

    // default growing in "time direction"
    Eigen::Vector3d direction(0.0, 0.0, 1.0);

    size_t size = line.GetSize();

    // compute the tangent instead of fixed direction when there is more than one point
    if (size > 1) {
        direction = (line.Read(size - 1) - line.Read(size - 2)).stableNormalized();
    }

    Eigen::Vector3d newPoint = line.Read(size - 1) + direction * stepSize;
    line.Append(newPoint);
}

// samples the image data using vtk (unused)
Eigen::Vector3d vtkCriticalLinesExtractor::SampleImageData(Eigen::Vector3d coords,
                                                           const char* arrayName) const {
    double tol2 = 0.0;
    vtkSmartPointer<vtkCell> cell;
    vtkSmartPointer<vtkPointData> pd;
    int subId;
    double pcoords[3] = {0, 0, 0};
    double weights[8] = {0, 0, 0, 0, 0, 0, 0, 0};

    pd = ImageData->GetPointData();
    vtkNew<vtkPointData> outPD;
    outPD->InterpolateAllocate(pd, 1, 1);

    cell = ImageData->FindAndGetCell(coords.data(), nullptr, -1, tol2, subId, pcoords, weights);
    //cell->PrintSelf(std::cout, vtkIndent(0));
    outPD->InterpolatePoint(pd, 0, cell->PointIds, weights);
    //cell->PrintSelf(std::cout, vtkIndent(0));
    //outPD->GetArray(arrayName.c_str())->PrintSelf(std::cout, vtkIndent(0));
    double* tuple = outPD->GetArray(arrayName)->GetTuple3(0);

    return { tuple[0], tuple[1], tuple[2] };

}

// samples the image data directly
Eigen::Vector3d vtkCriticalLinesExtractor::Sample(const Eigen::Vector3d& coordinate, const char* arrayName) const
{

    vtkSmartPointer<vtkDataArray> array = ImageData->GetPointData()->GetArray(arrayName);

    // Get the information on the domain
    int* dimensionsRaw = ImageData->GetDimensions();
    Eigen::Vector3i dimensions(dimensionsRaw[0], dimensionsRaw[1], dimensionsRaw[2]);
    double* spacingRaw = ImageData->GetSpacing();
    Eigen::Vector3d spacing(spacingRaw[0], spacingRaw[1], spacingRaw[2]);
    double* originRaw = ImageData->GetOrigin();
    Eigen::Vector3d domainMin(originRaw[0], originRaw[1], originRaw[2]);
    Eigen::Vector3d domainMax = domainMin + (dimensions.template cast<double>() - Eigen::Vector3d::Ones()).cwiseProduct(spacing);

    Eigen::Vector3d vf_tex = (coordinate - domainMin).cwiseQuotient(domainMax - domainMin);
    Eigen::Vector3d vf_sample = vf_tex.cwiseProduct(dimensions.template cast<double>() - Eigen::Vector3d::Ones());

    Eigen::Vector3i vi_sample_base0, vi_sample_base1;
    for (int i = 0; i < 3; ++i) {
        vi_sample_base0[i] = std::min(std::max(0, (int)vf_sample[i]), dimensions[i] - 1);
        vi_sample_base1[i] = std::min(std::max(0, vi_sample_base0[i] + 1), dimensions[i] - 1);
    }
    Eigen::Vector3d vf_sample_interpol = vf_sample - vi_sample_base0.template cast<double>();

    int num_corners = (int)std::pow(2, 3);
    Eigen::Vector3d result = Eigen::Vector3d::Zero();
    for (int i = 0; i < num_corners; ++i) {
        double weight = 1;
        Eigen::Vector3i grid_index = Eigen::Vector3i::Zero();
        for (int d = 0; d < 3; ++d) {
            if (i & (int(1) << (3 - 1 - d))) {
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
        double* value = array->GetTuple3(grid_index.x()
                + grid_index.y() * dimensions.x()
                + grid_index.z() * dimensions.x() * dimensions.y());
        result += Eigen::Vector3d(value[0], value[1], value[2]) * weight;
    }
    //std::cout << "##################################################" << std::endl;
    return result;
}