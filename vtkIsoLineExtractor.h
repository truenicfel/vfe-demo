#pragma once

#include <vtkImageData.h>
#include <vtkVector.h>
#include <vtkPolyDataAlgorithm.h>
#include "variational/StoppingCriteria.h"
#include "variational/Gradient.h"
#include "variational/DoubleBufferedLine.h"

class vtkIsoLineExtractor : public vtkPolyDataAlgorithm
{
public:

	static vtkIsoLineExtractor* New();
	vtkTypeMacro(vtkIsoLineExtractor, vtkPolyDataAlgorithm)

	void SetImageData(vtkSmartPointer<vtkImageData> imageData);

	vtkGetMacro(GrowingIterations, int)
	vtkSetMacro(GrowingIterations, int)

	vtkGetMacro(StepSizeGrowing, double)
	vtkSetMacro(StepSizeGrowing, double)

	vtkGetMacro(VariationalIterations, int)
	vtkSetMacro(VariationalIterations, int)

	vtkGetMacro(StepSizeVariational, double)
	vtkSetMacro(StepSizeVariational, double)

    vtkGetMacro(BorderRefinement, int)
    vtkSetMacro(BorderRefinement, int)

    vtkGetMacro(ProximityThreshold, double)
    vtkSetMacro(ProximityThreshold, double)

    vtkGetMacro(Iso, double)
    vtkSetMacro(Iso, double)

    vtkGetMacro(SmoothnessWeight, double)
    vtkSetMacro(SmoothnessWeight, double)

	vtkGetStringMacro(FieldNameScalars)
	vtkSetStringMacro(FieldNameScalars)

	vtkGetStringMacro(FieldNameGradient)
	vtkSetStringMacro(FieldNameGradient)


protected:

	vtkIsoLineExtractor();
	~vtkIsoLineExtractor() override = default;

	int RequestData(vtkInformation*, vtkInformationVector**, vtkInformationVector*) override;

	int GrowingIterations;
	double StepSizeGrowing;
	int VariationalIterations;
	double StepSizeVariational;
    int BorderRefinement;
    double ProximityThreshold;
    double Iso;
    double SmoothnessWeight;
	char* FieldNameScalars;
	char* FieldNameGradient;

	vtkSmartPointer<vtkImageData> ImageData;
	

private:
	vtkIsoLineExtractor(const vtkIsoLineExtractor&);  // Not implemented.
	void operator=(const vtkIsoLineExtractor&);  // Not implemented.

    std::deque<Eigen::Vector2d> extract(const Eigen::Vector2d& seed) const;

    void grow(DoubleBufferedLine2& line, const bool& forward = true) const;

    double SampleScalar(const Eigen::Vector2d& coordinate) const;

    Eigen::Vector2d SampleGradient(const Eigen::Vector2d& coordinate) const;

    double stepSizeGrowing;
	double stepSizeRefinement;

};
