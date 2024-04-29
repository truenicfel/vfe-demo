#pragma once

#include <vtkImageData.h>
#include <vtkVector.h>
#include <vtkPolyDataAlgorithm.h>
#include "variational/StoppingCriteria.h"
#include "variational/Gradient.h"
#include "variational/DoubleBufferedLine.h"

class vtkCriticalLinesExtractor : public vtkPolyDataAlgorithm
{
public:

	static vtkCriticalLinesExtractor* New();
	vtkTypeMacro(vtkCriticalLinesExtractor, vtkPolyDataAlgorithm)

	void SetImageData(vtkSmartPointer<vtkImageData> imageData);

	vtkGetMacro(GrowingIterations, int)
	vtkSetMacro(GrowingIterations, int)

	vtkGetMacro(TimeStep, double)
	vtkSetMacro(TimeStep, double)

	vtkGetMacro(VariationalIterations, int)
	vtkSetMacro(VariationalIterations, int)

	vtkGetMacro(AdamLearningRate, double)
	vtkSetMacro(AdamLearningRate, double)


	vtkGetMacro(UseComputeBoundingBox, bool)
	vtkSetMacro(UseComputeBoundingBox, bool)

	vtkGetMacro(BBMin, vtkVector3d)
	vtkSetMacro(BBMin, vtkVector3d)

	vtkGetMacro(BBMax, vtkVector3d)
	vtkSetMacro(BBMax, vtkVector3d)

	vtkGetStringMacro(FieldNameV)
	vtkSetStringMacro(FieldNameV)

	vtkGetStringMacro(FieldNameVx)
	vtkSetStringMacro(FieldNameVx)

	vtkGetStringMacro(FieldNameVy)
	vtkSetStringMacro(FieldNameVy)

	vtkGetStringMacro(FieldNameVt)
	vtkSetStringMacro(FieldNameVt)


protected:

	vtkCriticalLinesExtractor();
	~vtkCriticalLinesExtractor() override = default;

	int RequestData(vtkInformation*, vtkInformationVector**, vtkInformationVector*) override;

	int GrowingIterations;
	double TimeStep;
	int VariationalIterations;
	double AdamLearningRate;
	bool UseComputeBoundingBox;
	vtkVector3d BBMin;
	vtkVector3d BBMax;
	char* FieldNameV;
	char* FieldNameVx;
	char* FieldNameVy;
	char* FieldNameVt;

	vtkSmartPointer<vtkImageData> ImageData;
	

private:
	vtkCriticalLinesExtractor(const vtkCriticalLinesExtractor&);  // Not implemented.
	void operator=(const vtkCriticalLinesExtractor&);  // Not implemented.

    std::deque<Eigen::Vector3d> extract(const Eigen::Vector3d& seed) const;

    void grow(DoubleBufferedLine3& line, const bool& forward = true) const;

    Eigen::Vector3d SampleImageData(Eigen::Vector3d coords, const char* arrayName) const;

    Eigen::Vector3d Sample(const Eigen::Vector3d& coordinate, const char* arrayName) const;

    double stepSize;

};
